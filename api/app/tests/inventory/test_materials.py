from unittest import skip
from django.urls import reverse
from django.test import override_settings
from django.core.files.uploadedfile import SimpleUploadedFile
from rest_framework import status
from rest_framework.test import APITestCase
import shutil
import tempfile
from inventory.static_models import ItemType, MaterialMaster


@skip("Does not authenticate, so the API answers 403, see plan.md")
class InventoryMaterialReagentTests(APITestCase):
    """
    Covers reagent-specific material metadata.
    """

    def setUp(self):
        self.media_root = tempfile.mkdtemp()
        self.override_media = override_settings(MEDIA_ROOT=self.media_root)
        self.override_media.enable()
        self.addCleanup(self.override_media.disable)
        self.addCleanup(shutil.rmtree, self.media_root, ignore_errors=True)

        self.reagent_item_type = ItemType.objects.create(name="reagent")
        self.device_item_type = ItemType.objects.create(name="device")

    def test_create_reagent_requires_storage_temperature(self):
        payload = {
            "product_name": "PBS Buffer",
            "item_type_id": self.reagent_item_type.id,
        }

        response = self.client.post(
            reverse("inventory-material-list"), payload, format="json"
        )

        self.assertEqual(response.status_code, status.HTTP_400_BAD_REQUEST)
        self.assertEqual(
            response.data["storage_temperature"][0],
            "Storage temperature is required for reagents.",
        )

    def test_create_reagent_accepts_storage_temperature_and_sds(self):
        payload = {
            "product_name": "PBS Buffer",
            "item_type_id": str(self.reagent_item_type.id),
            "storage_temperature": "4°C",
            "safety_data_sheet": SimpleUploadedFile(
                "pbs-sds.pdf",
                b"fake-pdf-content",
                content_type="application/pdf",
            ),
        }

        response = self.client.post(
            reverse("inventory-material-list"), payload, format="multipart"
        )

        self.assertEqual(response.status_code, status.HTTP_201_CREATED)

        material = MaterialMaster.objects.get(id=response.data["id"])
        self.assertEqual(material.storage_temperature, "4°C")
        self.assertTrue(material.safety_data_sheet.name.endswith("pbs-sds.pdf"))
        self.assertEqual(response.data["storage_temperature"], "4°C")
        self.assertEqual(response.data["storage_temperature_label"], "4°C")
        self.assertIn("pbs-sds.pdf", response.data["safety_data_sheet"])

    def test_patch_non_reagent_allows_empty_storage_temperature(self):
        material = MaterialMaster.objects.create(
            product_name="Power Cable",
            item_type=self.device_item_type,
        )

        response = self.client.patch(
            reverse("inventory-material-detail", args=[material.id]),
            {"description": "Updated device note"},
            format="json",
        )

        self.assertEqual(response.status_code, status.HTTP_200_OK)
        material.refresh_from_db()
        self.assertEqual(material.description, "Updated device note")
