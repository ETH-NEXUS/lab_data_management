"""
The views of the management page need a logged in user and, for POST, a CSRF token.
"""

import shutil
import tempfile
from os.path import join

from django.contrib.auth.models import User
from django.test import Client, TestCase, override_settings
from django.urls import reverse

POST_VIEWS = [
    "run_command",
    "delete_file",
    "download_file",
    "upload_file",
    "get_file_content",
]


class LoginRequiredTest(TestCase):
    def test_without_login_every_view_is_refused(self):
        responses = [self.client.post(reverse(name)) for name in POST_VIEWS]
        responses.append(self.client.get(reverse("directory_content")))
        responses.append(self.client.get(reverse("long_polling", args=["room_1"])))

        self.assertEqual([403] * 7, [response.status_code for response in responses])

    def test_with_login_the_output_can_be_read(self):
        self.client.force_login(User.objects.create_user("tester"))

        response = self.client.get(reverse("long_polling", args=["room_1"]))

        self.assertEqual(200, response.status_code)

    def test_a_post_without_csrf_token_is_refused(self):
        client = Client(enforce_csrf_checks=True)
        client.force_login(User.objects.create_user("tester"))

        response = client.post(
            reverse("run_command"),
            {"form_data": {"command": "map", "room_name": "room_1"}},
            content_type="application/json",
        )

        self.assertEqual(403, response.status_code)
        self.assertIn("CSRF", response.json()["detail"])

    def test_a_post_with_login_and_csrf_token_works_like_the_ui(self):
        # The UI gets the CSRF cookie first and sends the token in a header
        client = Client(enforce_csrf_checks=True)
        client.force_login(User.objects.create_user("tester"))
        client.get(reverse("auth-cookie"))
        token = client.cookies["csrftoken"].value
        folder = tempfile.mkdtemp()
        self.addCleanup(shutil.rmtree, folder)
        path = join(folder, "notes.txt")
        with open(path, "w") as file:
            file.write("hello")

        with override_settings(MANAGEMENT_DATA_ROOT=folder):
            response = client.post(
                reverse("get_file_content"),
                {"file_path": path},
                content_type="application/json",
                HTTP_X_CSRFTOKEN=token,
            )

        self.assertEqual(200, response.status_code)
        self.assertEqual({"content": "hello"}, response.json())
