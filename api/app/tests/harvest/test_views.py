"""
The Harvest views: the list of Harvest projects and taking over the name and
notes of a Harvest project into a project of this app.

The Harvest client is always replaced by a mock, so no test calls Harvest.
"""

from unittest import mock

import requests
from django.contrib.auth.models import User
from django.test import TestCase, override_settings
from django.urls import reverse

from core.models import Project

# Shortened projects as GET https://api.harvestapp.com/v2/projects lists them
SCREENING = {"id": 7, "name": "SNL Screening 2026", "notes": "Notes from Harvest"}
OTHER = {"id": 8, "name": "Other project", "notes": None}


def harvest_error(status_code):
    """The error HarvestClient.get raises when Harvest answers with this status."""
    response = requests.Response()
    response.status_code = status_code
    return requests.HTTPError(f"{status_code} from Harvest", response=response)


def harvest_answer(endpoint, params=None):
    """Answers like Harvest: "projects" lists them, "projects/7" is one project."""
    if endpoint == "projects":
        return {"projects": [SCREENING, OTHER], "total_pages": 1}
    for project in (SCREENING, OTHER):
        if endpoint == f"projects/{project['id']}":
            return project
    raise harvest_error(404)


@override_settings(HARVEST_ACCESS_TOKEN="token", HARVEST_PROJECT_FILTER=None)
class HarvestViewsTest(TestCase):
    def setUp(self):
        self.client.force_login(User.objects.create_user("tester"))
        patcher = mock.patch("harvest.views.client.get", side_effect=harvest_answer)
        self.harvest_get = patcher.start()
        self.addCleanup(patcher.stop)

    def projects(self):
        return self.client.get(reverse("harvest_projects"))

    def update(self, project):
        return self.client.get(reverse("update_harvest_info", args=[project.id]))

    def test_the_harvest_projects_are_listed(self):
        response = self.projects()

        self.assertEqual({"projects": [SCREENING, OTHER]}, response.json())

    @override_settings(HARVEST_PROJECT_FILTER="SNL")
    def test_only_projects_with_the_filter_text_in_their_name_are_listed(self):
        response = self.projects()

        self.assertEqual({"projects": [SCREENING]}, response.json())

    @override_settings(HARVEST_ACCESS_TOKEN=None)
    def test_without_a_token_the_list_is_empty_and_harvest_is_not_asked(self):
        response = self.projects()

        self.assertEqual({"projects": []}, response.json())
        self.harvest_get.assert_not_called()

    def test_the_list_says_so_when_harvest_cannot_be_reached(self):
        self.harvest_get.side_effect = requests.ConnectionError("no network")

        response = self.projects()

        self.assertEqual(502, response.status_code)
        self.assertIn("Harvest could not be asked", response.json()["error"])

    def test_the_name_and_notes_are_taken_from_harvest(self):
        project = Project.objects.create(name="Old name", harvest_id=7)

        response = self.update(project)

        self.assertEqual({"success": True}, response.json())
        project.refresh_from_db()
        self.assertEqual("SNL Screening 2026", project.name)
        self.assertEqual("Notes from Harvest", project.harvest_notes)
        self.harvest_get.assert_called_once_with("projects/7")

    def test_a_project_without_harvest_id_stays_as_it_is(self):
        project = Project.objects.create(name="Manual project")

        response = self.update(project)

        self.assertEqual({"success": True}, response.json())
        project.refresh_from_db()
        self.assertEqual("Manual project", project.name)
        self.assertIsNone(project.harvest_notes)
        self.harvest_get.assert_not_called()

    def test_a_project_that_is_no_longer_in_harvest_answers_404(self):
        project = Project.objects.create(name="Old name", harvest_id=99)

        response = self.update(project)

        self.assertEqual(404, response.status_code)
        self.assertEqual(
            {"error": "Project 99 is no longer in Harvest."}, response.json()
        )
        project.refresh_from_db()
        self.assertEqual("Old name", project.name)

    def test_an_update_says_so_when_harvest_refuses_the_token(self):
        self.harvest_get.side_effect = harvest_error(401)
        project = Project.objects.create(name="Old name", harvest_id=7)

        response = self.update(project)

        self.assertEqual(502, response.status_code)
        project.refresh_from_db()
        self.assertEqual("Old name", project.name)

    def test_an_update_says_so_when_harvest_cannot_be_reached(self):
        self.harvest_get.side_effect = requests.Timeout("too slow")
        project = Project.objects.create(name="Old name", harvest_id=7)

        response = self.update(project)

        self.assertEqual(502, response.status_code)

    @override_settings(HARVEST_ACCESS_TOKEN=None)
    def test_an_update_without_a_token_says_so_and_harvest_is_not_asked(self):
        project = Project.objects.create(name="Old name", harvest_id=7)

        response = self.update(project)

        self.assertEqual(400, response.status_code)
        self.assertEqual(
            {"error": "Harvest is not set up on this server."}, response.json()
        )
        self.harvest_get.assert_not_called()

    def test_an_unknown_project_answers_404(self):
        response = self.client.get(reverse("update_harvest_info", args=[999999]))

        self.assertEqual(404, response.status_code)

    def test_without_login_both_views_are_refused_and_harvest_is_not_asked(self):
        self.client.logout()
        project = Project.objects.create(name="Old name", harvest_id=7)

        responses = [self.projects(), self.update(project)]

        self.assertEqual([403, 403], [response.status_code for response in responses])
        self.harvest_get.assert_not_called()
        project.refresh_from_db()
        self.assertEqual("Old name", project.name)
