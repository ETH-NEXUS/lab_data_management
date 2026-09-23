"""
The Harvest views: the list of Harvest projects and taking over the name and
notes of a Harvest project into a project of this app.

The Harvest client is always replaced by a mock, so no test calls Harvest.
"""

from unittest import mock

from django.contrib.auth.models import User
from django.test import TestCase, override_settings
from django.urls import reverse

from core.models import Project

# A shortened answer of GET https://api.harvestapp.com/v2/projects
HARVEST_PROJECTS = {
    "projects": [
        {"id": 7, "name": "Screening 2026", "notes": "Notes from Harvest"},
        {"id": 8, "name": "Other project", "notes": None},
    ]
}


class HarvestViewsTest(TestCase):
    def setUp(self):
        self.client.force_login(User.objects.create_user("tester"))
        patcher = mock.patch("harvest.views.client.get", return_value=HARVEST_PROJECTS)
        self.harvest_get = patcher.start()
        self.addCleanup(patcher.stop)

    def update(self, project):
        return self.client.get(reverse("update_harvest_info", args=[project.id]))

    @override_settings(HARVEST_ACCESS_TOKEN=None)
    def test_without_a_token_the_list_is_empty_and_harvest_is_not_asked(self):
        response = self.client.get(reverse("harvest_projects"))

        self.assertEqual({"projects": []}, response.json())
        self.harvest_get.assert_not_called()

    def test_the_name_and_notes_are_taken_from_harvest(self):
        project = Project.objects.create(name="Old name", harvest_id=7)

        response = self.update(project)

        self.assertEqual({"success": True}, response.json())
        project.refresh_from_db()
        self.assertEqual("Screening 2026", project.name)
        self.assertEqual("Notes from Harvest", project.harvest_notes)
        self.harvest_get.assert_called_once_with("projects")

    def test_a_project_without_harvest_id_stays_as_it_is(self):
        project = Project.objects.create(name="Manual project")

        response = self.update(project)

        self.assertEqual({"success": True}, response.json())
        project.refresh_from_db()
        self.assertEqual("Manual project", project.name)
        self.assertIsNone(project.harvest_notes)

    def test_an_unknown_project_answers_404(self):
        response = self.client.get(reverse("update_harvest_info", args=[999999]))

        self.assertEqual(404, response.status_code)
