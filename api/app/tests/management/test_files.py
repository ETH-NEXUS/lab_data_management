"""
The file views of the management page: the folder tree, reading, downloading,
uploading and deleting files of the data folder.
"""

import os
from os.path import join

from django.core.files.uploadedfile import SimpleUploadedFile
from django.urls import reverse

from importer.command_output import add_message, start_command
from tests.management.base import ROOM_NAME, ManagementPageTestCase

# German text as the lab writes it, long enough for the encoding to be recognized
GERMAN_TEXT = (
    "Die Lösung wurde für die Prüfung verdünnt. Größe, Maß und Gefäß "
    "wurden nach der Überprüfung übernommen.\n"
) * 5


def names(children):
    """The names of the entries of a folder, sorted: the order on disk is not fixed."""
    return sorted(child["name"] for child in children)


def child_named(children, name):
    return next(child for child in children if child["name"] == name)


class DirectoryContentTest(ManagementPageTestCase):
    def directory_content(self):
        return self.client.get(reverse("directory_content")).json()["directory_content"]

    def test_the_tree_lists_files_and_folders(self):
        self.write("a.csv", "x")
        os.mkdir(join(self.folder, "echo"))
        self.write(join("echo", "run.xml"), "x")

        tree = self.directory_content()

        self.assertEqual("directory", tree["type"])
        self.assertEqual(self.folder, tree["path"])
        self.assertEqual(["a.csv", "echo"], names(tree["children"]))
        self.assertEqual(
            {"type": "file", "name": "a.csv", "path": join(self.folder, "a.csv")},
            child_named(tree["children"], "a.csv"),
        )
        echo = child_named(tree["children"], "echo")
        self.assertEqual("directory", echo["type"])
        self.assertEqual(["run.xml"], names(echo["children"]))

    def test_a_snapshots_folder_is_listed_without_its_content(self):
        os.mkdir(join(self.folder, ".snapshots"))
        self.write(join(".snapshots", "old.csv"), "x")

        tree = self.directory_content()

        snapshots = child_named(tree["children"], ".snapshots")
        self.assertEqual([], snapshots["children"])

    def test_a_link_to_a_folder_is_not_followed(self):
        # A link back to the data folder would list it again and again
        os.symlink(self.folder, join(self.folder, "loop"))

        tree = self.directory_content()

        self.assertEqual([], tree["children"])

    def test_a_link_to_a_file_is_listed_as_a_file(self):
        target = self.write("a.csv", "x")
        os.symlink(target, join(self.folder, "link.csv"))

        tree = self.directory_content()

        self.assertEqual(["a.csv", "link.csv"], names(tree["children"]))
        self.assertEqual("file", child_named(tree["children"], "link.csv")["type"])


class FileViewsTest(ManagementPageTestCase):
    def post(self, name, data):
        return self.client.post(reverse(name), data, content_type="application/json")

    def test_a_file_is_deleted(self):
        path = self.write("a.csv", "x")

        response = self.post("delete_file", {"path": path})

        self.assertEqual({"status": "ok"}, response.json())
        self.assertFalse(os.path.exists(path))

    def test_deleting_a_missing_file_is_not_an_error(self):
        response = self.post("delete_file", {"path": join(self.folder, "missing.csv")})

        self.assertEqual({"status": "ok"}, response.json())

    def test_a_file_is_downloaded_with_its_name(self):
        path = self.write("a.csv", "1,2,3")

        response = self.post("download_file", {"file_path": path})

        self.assertEqual(200, response.status_code)
        self.assertEqual(b"1,2,3", b"".join(response))
        self.assertEqual(
            'attachment; filename="a.csv"', response["Content-Disposition"]
        )

    def test_downloading_a_missing_file_answers_404(self):
        response = self.post(
            "download_file", {"file_path": join(self.folder, "no.csv")}
        )

        self.assertEqual(404, response.status_code)

    def test_a_file_with_a_german_name_is_downloaded_with_its_name(self):
        path = self.write("Lösung.csv", "1")

        response = self.post("download_file", {"file_path": path})

        self.assertEqual(200, response.status_code)
        self.assertEqual(
            "attachment; filename*=utf-8''L%C3%B6sung.csv",
            response["Content-Disposition"],
        )

    def test_a_folder_is_not_downloaded(self):
        response = self.post("download_file", {"file_path": self.folder})

        self.assertEqual(400, response.status_code)
        self.assertIn("is a folder", response.json()[0])

    def test_a_folder_is_not_read(self):
        response = self.post("get_file_content", {"file_path": self.folder})

        self.assertEqual(400, response.status_code)
        self.assertIn("is a folder", response.json()[0])

    def test_a_file_is_uploaded_into_a_new_folder(self):
        folder = join(self.folder, "new", "folder")

        response = self.client.post(
            reverse("upload_file"),
            {"directory_path": folder, "file": SimpleUploadedFile("a.csv", b"1,2")},
        )

        self.assertEqual(join(folder, "a.csv"), response.json()["file_path"])
        with open(join(folder, "a.csv"), "rb") as file:
            self.assertEqual(b"1,2", file.read())

    def test_an_upload_does_not_replace_a_file_with_the_same_name(self):
        path = self.write("a.csv", "old")

        response = self.client.post(
            reverse("upload_file"),
            {
                "directory_path": self.folder,
                "file": SimpleUploadedFile("a.csv", b"new"),
            },
        )

        self.assertEqual(400, response.status_code)
        self.assertEqual(
            [
                f"A file named a.csv already exists in {self.folder}. "
                "Delete it first or rename the new file."
            ],
            response.json(),
        )
        with open(path) as file:
            self.assertEqual("old", file.read())

    def test_an_upload_does_not_replace_a_folder_with_the_same_name(self):
        os.mkdir(join(self.folder, "echo"))

        response = self.client.post(
            reverse("upload_file"),
            {"directory_path": self.folder, "file": SimpleUploadedFile("echo", b"x")},
        )

        self.assertEqual(400, response.status_code)
        self.assertTrue(os.path.isdir(join(self.folder, "echo")))

    def test_an_upload_without_a_file_says_so(self):
        response = self.client.post(
            reverse("upload_file"), {"directory_path": self.folder}
        )

        self.assertEqual(400, response.status_code)
        self.assertEqual(
            ["The request has no directory path or no file."], response.json()
        )

    def test_a_utf8_file_is_read(self):
        path = join(self.folder, "utf8.txt")
        with open(path, "w", encoding="utf-8") as file:
            file.write(GERMAN_TEXT)

        response = self.post("get_file_content", {"file_path": path})

        self.assertEqual({"content": GERMAN_TEXT}, response.json())

    def test_a_latin1_file_is_read(self):
        # Older instrument software writes Windows/latin-1 files
        path = join(self.folder, "latin1.txt")
        with open(path, "w", encoding="latin-1") as file:
            file.write(GERMAN_TEXT)

        response = self.post("get_file_content", {"file_path": path})

        self.assertEqual({"content": GERMAN_TEXT}, response.json())

    def test_reading_a_missing_file_answers_404(self):
        response = self.post("get_file_content", {"file_path": join(self.folder, "no")})

        self.assertEqual(404, response.status_code)


class LongPollingTest(ManagementPageTestCase):
    def setUp(self):
        super().setUp()
        start_command(ROOM_NAME)
        add_message(ROOM_NAME, "info", "first")
        add_message(ROOM_NAME, "info", "second")

    def test_only_the_messages_from_since_on_are_returned(self):
        output = self.read_output(since=1)

        self.assertEqual([{"level": "info", "text": "second"}], output["messages"])
        self.assertEqual(2, output["next"])
        self.assertEqual("running", output["status"])

    def test_a_since_that_is_not_a_number_starts_from_the_first_message(self):
        url = reverse("long_polling", args=[ROOM_NAME])

        output = self.client.get(f"{url}?since=abc").json()

        self.assertEqual(["first", "second"], [m["text"] for m in output["messages"]])

    def test_a_room_without_a_command_has_no_messages_and_no_status(self):
        url = reverse("long_polling", args=["other_room"])

        output = self.client.get(url).json()

        self.assertEqual({"messages": [], "next": 0, "status": None}, output)

    def test_a_room_name_that_is_not_a_name_is_refused(self):
        url = reverse("long_polling", args=["not-a-room!"])

        response = self.client.get(url)

        self.assertEqual(400, response.status_code)
        self.assertEqual(["This is not a room name: not-a-room!"], response.json())
