from datetime import timedelta

from django.contrib.auth import get_user_model
from django.urls import reverse
from django.test import override_settings
from django.core.files.uploadedfile import SimpleUploadedFile
from django.utils import timezone
from rest_framework import status
from rest_framework.test import APITestCase
import shutil
import tempfile

from core.models import Experiment, Project
from inventory.dynamic_models import InventoryDashboardTilePreference, InventoryStock, MaterialUsage, Order, Room, Sector
from inventory.history_models import InventoryChangeRecord
from inventory.static_models import ItemType, MaterialMaster, MaterialUnit, UnitOfMeasure

User = get_user_model()


class InventoryStockMultiSectorTests(APITestCase):
    """
    Covers the backend stock flow for one item stored in multiple sectors.
    """

    def setUp(self):
        self.first_user = User.objects.create_user(
            username="first-inventory-user",
            password="test-password",
        )
        self.second_user = User.objects.create_user(
            username="second-inventory-user",
            password="test-password",
        )
        self.client.force_authenticate(user=self.first_user)
        self.item_type = ItemType.objects.create(name="Consumable")
        self.unit_of_measure = UnitOfMeasure.objects.create(name="box")
        self.material = MaterialMaster.objects.create(
            product_name="Reservoir tips",
            item_type=self.item_type,
        )
        self.stock_unit = MaterialUnit.objects.create(
            material=self.material,
            unit=self.unit_of_measure,
            base_units_per_unit="1",
            is_stock_unit=True,
        )
        self.room = Room.objects.create(name="C75")
        self.primary_sector = Sector.objects.create(room=self.room, name="3.1")
        self.secondary_sector = Sector.objects.create(room=self.room, name="3.2")
        self.other_room = Room.objects.create(name="C41")
        self.other_room_sector = Sector.objects.create(room=self.other_room, name="1.1")
        self.project = Project.objects.create(name="Inventory test project")
        self.order = Order.objects.create(
            material=self.material,
            order_unit=self.stock_unit,
            amount="2",
            order_date="2026-07-05T10:00:00Z",
            status=Order.STATUS_ORDERED,
            project=self.project,
        )

    def test_create_stock_accepts_multiple_sector_ids(self):
        payload = {
            "material_id": self.material.id,
            "sector_ids": [self.primary_sector.id, self.secondary_sector.id],
            "stock_unit_id": self.stock_unit.id,
            "quantity": "2",
            "minimum_quantity": "1",
            "notes": "Split across shelves",
        }

        response = self.client.post(reverse("inventory-stock-list"), payload, format="json")

        self.assertEqual(response.status_code, status.HTTP_201_CREATED)

        stock = InventoryStock.objects.get(id=response.data["id"])
        self.assertEqual(stock.sector_id, self.primary_sector.id)
        self.assertEqual(
            list(stock.additional_sectors.order_by("id").values_list("id", flat=True)),
            [self.secondary_sector.id],
        )
        self.assertEqual(response.data["location_label"], "C75 / 3.1, 3.2")
        self.assertEqual(
            [sector["id"] for sector in response.data["sectors"]],
            [self.primary_sector.id, self.secondary_sector.id],
        )

    def test_favorites_only_include_the_current_users_entries(self):
        stock = InventoryStock.objects.create(
            material=self.material,
            sector=self.primary_sector,
            stock_unit=self.stock_unit,
            quantity="2",
            minimum_quantity="1",
        )

        self.client.force_authenticate(user=self.first_user)
        mark_response = self.client.post(reverse("inventory-stock-mark-favorite", args=[stock.id]))
        own_favorites_response = self.client.get(reverse("inventory-stock-favorites"))

        self.client.force_authenticate(user=self.second_user)
        other_favorites_response = self.client.get(reverse("inventory-stock-favorites"))
        other_detail_response = self.client.get(reverse("inventory-stock-detail", args=[stock.id]))

        self.assertEqual(mark_response.status_code, status.HTTP_200_OK)
        self.assertTrue(mark_response.data["is_favorite"])
        self.assertEqual(own_favorites_response.data["count"], 1)
        self.assertEqual(own_favorites_response.data["results"][0]["id"], stock.id)
        self.assertEqual(other_favorites_response.data["count"], 0)
        self.assertFalse(other_detail_response.data["is_favorite"])

    def test_stock_list_sorts_favorites_for_current_user(self):
        non_favorite_stock = InventoryStock.objects.create(
            material=self.material,
            sector=self.primary_sector,
            stock_unit=self.stock_unit,
            quantity="2",
            minimum_quantity="1",
        )
        favorite_stock = InventoryStock.objects.create(
            material=self.material,
            sector=self.primary_sector,
            stock_unit=self.stock_unit,
            quantity="2",
            minimum_quantity="1",
        )
        favorite_stock.favorite_users.add(self.first_user)

        response = self.client.get(
            reverse("inventory-stock-list"),
            {"ordering": "-favorite"},
        )

        self.assertEqual(response.status_code, status.HTTP_200_OK)
        self.assertEqual(
            [stock_data["id"] for stock_data in response.data["results"]],
            [favorite_stock.id, non_favorite_stock.id],
        )

    def test_create_stock_rejects_sectors_from_multiple_rooms(self):
        payload = {
            "material_id": self.material.id,
            "sector_ids": [self.primary_sector.id, self.other_room_sector.id],
            "stock_unit_id": self.stock_unit.id,
            "quantity": "2",
            "minimum_quantity": "1",
        }

        response = self.client.post(reverse("inventory-stock-list"), payload, format="json")

        self.assertEqual(response.status_code, status.HTTP_400_BAD_REQUEST)
        self.assertEqual(
            response.data["sector_ids"][0],
            "All selected sectors must belong to the same room.",
        )

    def test_sector_filter_matches_primary_and_additional_sectors(self):
        stock = InventoryStock.objects.create(
            material=self.material,
            sector=self.primary_sector,
            stock_unit=self.stock_unit,
            quantity="2",
            minimum_quantity="1",
        )
        stock.additional_sectors.add(self.secondary_sector)

        response = self.client.get(
            reverse("inventory-stock-list"),
            {"sector": self.secondary_sector.id},
        )

        self.assertEqual(response.status_code, status.HTTP_200_OK)
        self.assertEqual(response.data["count"], 1)
        self.assertEqual(response.data["results"][0]["id"], stock.id)

    def test_patch_notes_keeps_existing_additional_sectors(self):
        stock = InventoryStock.objects.create(
            material=self.material,
            sector=self.primary_sector,
            stock_unit=self.stock_unit,
            quantity="2",
            minimum_quantity="1",
            notes="Initial note",
        )
        stock.additional_sectors.add(self.secondary_sector)

        response = self.client.patch(
            reverse("inventory-stock-detail", args=[stock.id]),
            {"notes": "Updated note"},
            format="json",
        )

        self.assertEqual(response.status_code, status.HTTP_200_OK)

        stock.refresh_from_db()
        self.assertEqual(stock.notes, "Updated note")
        self.assertEqual(
            list(stock.additional_sectors.order_by("id").values_list("id", flat=True)),
            [self.secondary_sector.id],
        )

    def test_patch_sector_ids_replaces_multi_sector_assignment(self):
        stock = InventoryStock.objects.create(
            material=self.material,
            sector=self.primary_sector,
            stock_unit=self.stock_unit,
            quantity="2",
            minimum_quantity="1",
        )
        stock.additional_sectors.add(self.secondary_sector)

        replacement_sector = Sector.objects.create(room=self.room, name="3.3")

        response = self.client.patch(
            reverse("inventory-stock-detail", args=[stock.id]),
            {"sector_ids": [replacement_sector.id]},
            format="json",
        )

        self.assertEqual(response.status_code, status.HTTP_200_OK)

        stock.refresh_from_db()
        self.assertEqual(stock.sector_id, replacement_sector.id)
        self.assertEqual(stock.additional_sectors.count(), 0)
        self.assertEqual(response.data["location_label"], "C75 / 3.3")

    def test_create_stock_accepts_source_order_id(self):
        payload = {
            "material_id": self.material.id,
            "sector_ids": [self.primary_sector.id],
            "stock_unit_id": self.stock_unit.id,
            "quantity": "2",
            "minimum_quantity": "1",
            "source_order_id": self.order.id,
        }

        response = self.client.post(reverse("inventory-stock-list"), payload, format="json")

        self.assertEqual(response.status_code, status.HTTP_201_CREATED)

        stock = InventoryStock.objects.get(id=response.data["id"])
        history_record = stock.change_records.get(
            performed_action=InventoryChangeRecord.ACTION_STOCK_CREATED,
        )
        self.assertEqual(stock.source_order_id, self.order.id)
        self.assertEqual(response.data["source_order"]["id"], self.order.id)
        self.assertEqual(history_record.performed_by, self.first_user)
        self.assertEqual(history_record.order_id, self.order.id)
        self.assertEqual(history_record.project_id, self.project.id)

    def test_archive_sets_timestamp_and_creates_history_record(self):
        stock = InventoryStock.objects.create(
            material=self.material,
            sector=self.primary_sector,
            stock_unit=self.stock_unit,
            quantity="2",
            minimum_quantity="1",
        )

        response = self.client.post(reverse("inventory-stock-archive", args=[stock.id]))

        self.assertEqual(response.status_code, status.HTTP_200_OK)
        self.assertTrue(response.data["is_archived"])
        self.assertIsNotNone(response.data["archived_at"])

        stock.refresh_from_db()
        self.assertIsNotNone(stock.archived_at)
        self.assertTrue(
            stock.change_records.filter(performed_action="stock_archived").exists()
        )

    def test_archived_endpoint_orders_by_archive_timestamp_and_restore_clears_it(self):
        earlier_stock = InventoryStock.objects.create(
            material=self.material,
            sector=self.primary_sector,
            stock_unit=self.stock_unit,
            quantity="2",
            minimum_quantity="1",
        )
        later_stock = InventoryStock.objects.create(
            material=self.material,
            sector=self.primary_sector,
            stock_unit=self.stock_unit,
            quantity="2",
            minimum_quantity="1",
        )

        self.client.post(reverse("inventory-stock-archive", args=[earlier_stock.id]))
        self.client.post(reverse("inventory-stock-archive", args=[later_stock.id]))
        earlier_stock.archived_at = timezone.now() - timedelta(days=1)
        earlier_stock.save(update_fields=["archived_at"])

        archived_response = self.client.get(reverse("inventory-stock-archived"))
        restore_response = self.client.post(reverse("inventory-stock-restore", args=[later_stock.id]))

        self.assertEqual(archived_response.status_code, status.HTTP_200_OK)
        self.assertEqual(
            [stock_data["id"] for stock_data in archived_response.data["results"]],
            [later_stock.id, earlier_stock.id],
        )
        self.assertEqual(restore_response.status_code, status.HTTP_200_OK)
        self.assertIsNone(restore_response.data["archived_at"])

        later_stock.refresh_from_db()
        self.assertFalse(later_stock.is_archived)
        self.assertIsNone(later_stock.archived_at)

    def test_awaiting_check_in_endpoint_excludes_orders_with_linked_stock(self):
        awaiting_check_in_order = Order.objects.create(
            material=self.material,
            order_unit=self.stock_unit,
            amount="2",
            order_date="2026-07-05T10:00:00Z",
            status=Order.STATUS_PRODUCT_ARRIVED,
        )
        checked_in_order = Order.objects.create(
            material=self.material,
            order_unit=self.stock_unit,
            amount="2",
            order_date="2026-07-04T10:00:00Z",
            status=Order.STATUS_PRODUCT_ARRIVED,
        )
        InventoryStock.objects.create(
            material=self.material,
            sector=self.primary_sector,
            stock_unit=self.stock_unit,
            quantity="2",
            minimum_quantity="1",
            source_order=checked_in_order,
        )

        response = self.client.get(reverse("inventory-order-awaiting-check-in"))

        self.assertEqual(response.status_code, status.HTTP_200_OK)
        self.assertEqual([order_data["id"] for order_data in response.data], [awaiting_check_in_order.id])

    def test_awaiting_check_in_endpoint_returns_only_the_latest_five_orders(self):
        for day in range(1, 7):
            Order.objects.create(
                material=self.material,
                order_unit=self.stock_unit,
                amount="2",
                order_date=f"2026-07-0{day}T10:00:00Z",
                status=Order.STATUS_PRODUCT_ARRIVED,
            )

        response = self.client.get(reverse("inventory-order-recent-awaiting-check-in"))

        self.assertEqual(response.status_code, status.HTTP_200_OK)
        self.assertEqual(len(response.data), 5)
        self.assertEqual(response.data[0]["order_date"][:10], "2026-07-06")

    def test_recent_project_usages_endpoint_returns_only_the_latest_five_usages(self):
        self.project.harvest_id = 42
        self.project.save(update_fields=["harvest_id"])
        stock = InventoryStock.objects.create(
            material=self.material,
            sector=self.primary_sector,
            stock_unit=self.stock_unit,
            quantity="12",
            minimum_quantity="1",
        )

        for day in range(1, 7):
            usage = MaterialUsage.objects.create(
                project=self.project,
                inventory_stock=stock,
                usage_unit=self.stock_unit,
                quantity_used="1",
            )
            usage.used_at = timezone.now() - timedelta(days=day)
            usage.save(update_fields=["used_at"])

        non_harvest_project = Project.objects.create(name="Non-Harvest project")
        MaterialUsage.objects.create(
            project=non_harvest_project,
            inventory_stock=stock,
            usage_unit=self.stock_unit,
            quantity_used="1",
        )

        response = self.client.get(reverse("inventory-material-usage-recent-project"))

        self.assertEqual(response.status_code, status.HTTP_200_OK)
        self.assertEqual(len(response.data), 5)
        self.assertEqual(response.data[0]["used_at"][:10], str((timezone.now() - timedelta(days=1)).date()))
        self.assertTrue(all(usage["project"]["id"] == self.project.id for usage in response.data))

    def test_recent_experiment_usages_endpoint_returns_only_the_latest_five_usages(self):
        stock = InventoryStock.objects.create(
            material=self.material,
            sector=self.primary_sector,
            stock_unit=self.stock_unit,
            quantity="12",
            minimum_quantity="1",
        )
        experiment = Experiment.objects.create(
            name="Inventory test experiment",
            project=self.project,
        )

        for day in range(1, 7):
            usage = MaterialUsage.objects.create(
                project=self.project,
                experiment=experiment,
                inventory_stock=stock,
                usage_unit=self.stock_unit,
                quantity_used="1",
            )
            usage.used_at = timezone.now() - timedelta(days=day)
            usage.save(update_fields=["used_at"])

        response = self.client.get(reverse("inventory-material-usage-recent-experiment"))

        self.assertEqual(response.status_code, status.HTTP_200_OK)
        self.assertEqual(len(response.data), 5)
        self.assertEqual(response.data[0]["used_at"][:10], str((timezone.now() - timedelta(days=1)).date()))
        self.assertTrue(all(usage["experiment"]["id"] == experiment.id for usage in response.data))

    def test_dashboard_tile_preferences_are_created_for_the_current_user(self):
        response = self.client.get(reverse("inventory-dashboard-tile-preference-current"))

        self.assertEqual(response.status_code, status.HTTP_200_OK)
        self.assertEqual(len(response.data), 10)
        self.assertEqual(sum(preference["is_visible"] for preference in response.data), 6)
        self.assertEqual(
            InventoryDashboardTilePreference.objects.filter(user=self.first_user).count(),
            10,
        )

    def test_dashboard_tile_preferences_are_saved_per_user(self):
        selected_tile_keys = [
            "low_stock_items",
            "awaiting_check_in",
            "favorite_items",
            "device_items",
        ]

        response = self.client.put(
            reverse("inventory-dashboard-tile-preference-current"),
            {"tile_keys": selected_tile_keys},
            format="json",
        )

        self.assertEqual(response.status_code, status.HTTP_200_OK)
        visible_tile_keys = [preference["key"] for preference in response.data if preference["is_visible"]]
        self.assertEqual(visible_tile_keys, selected_tile_keys)

        self.client.force_authenticate(user=self.second_user)
        second_user_response = self.client.get(reverse("inventory-dashboard-tile-preference-current"))

        self.assertEqual(second_user_response.status_code, status.HTTP_200_OK)
        self.assertEqual(sum(preference["is_visible"] for preference in second_user_response.data), 6)

    def test_dashboard_tile_preferences_allow_an_empty_selection(self):
        response = self.client.put(
            reverse("inventory-dashboard-tile-preference-current"),
            {"tile_keys": []},
            format="json",
        )

        self.assertEqual(response.status_code, status.HTTP_200_OK)
        self.assertEqual(sum(preference["is_visible"] for preference in response.data), 0)

    def test_dashboard_tile_preferences_allow_all_available_tiles(self):
        current_preferences_response = self.client.get(reverse("inventory-dashboard-tile-preference-current"))
        tile_keys = [preference["key"] for preference in current_preferences_response.data]

        response = self.client.put(
            reverse("inventory-dashboard-tile-preference-current"),
            {"tile_keys": tile_keys},
            format="json",
        )

        self.assertEqual(response.status_code, status.HTTP_200_OK)
        self.assertEqual(sum(preference["is_visible"] for preference in response.data), len(tile_keys))

    def test_dashboard_tile_preferences_reject_duplicate_tile_keys(self):
        response = self.client.put(
            reverse("inventory-dashboard-tile-preference-current"),
            {"tile_keys": ["low_stock_items", "low_stock_items"]},
            format="json",
        )

        self.assertEqual(response.status_code, status.HTTP_400_BAD_REQUEST)

    def test_dashboard_tile_preferences_reject_unknown_tile_keys(self):
        response = self.client.put(
            reverse("inventory-dashboard-tile-preference-current"),
            {"tile_keys": ["not_an_inventory_dashboard_tile"]},
            format="json",
        )

        self.assertEqual(response.status_code, status.HTTP_400_BAD_REQUEST)

    def test_history_endpoint_paginates_results(self):
        for _index in range(6):
            InventoryChangeRecord.objects.create(
                performed_action=InventoryChangeRecord.ACTION_STOCK_CREATED,
                performed_by=self.first_user,
            )

        response = self.client.get(
            reverse("inventory-history-list"),
            {"page": 1, "page_size": 5},
        )

        self.assertEqual(response.status_code, status.HTTP_200_OK)
        self.assertEqual(response.data["count"], 6)
        self.assertEqual(len(response.data["results"]), 5)
        self.assertIsNotNone(response.data["next"])

    def test_history_endpoint_filters_check_in_and_check_out_records(self):
        check_in_record = InventoryChangeRecord.objects.create(
            performed_action=InventoryChangeRecord.ACTION_STOCK_CREATED,
            performed_by=self.first_user,
            quantity_delta="2",
        )
        stock_increase_record = InventoryChangeRecord.objects.create(
            performed_action=InventoryChangeRecord.ACTION_STOCK_UPDATED,
            performed_by=self.first_user,
            quantity_delta="1",
        )
        check_out_record = InventoryChangeRecord.objects.create(
            performed_action=InventoryChangeRecord.ACTION_USAGE_CREATED,
            performed_by=self.first_user,
            quantity_delta="1",
        )
        InventoryChangeRecord.objects.create(
            performed_action=InventoryChangeRecord.ACTION_STOCK_UPDATED,
            performed_by=self.first_user,
            quantity_delta="-1",
        )
        InventoryChangeRecord.objects.create(
            performed_action=InventoryChangeRecord.ACTION_STOCK_ARCHIVED,
            performed_by=self.first_user,
        )

        response = self.client.get(
            reverse("inventory-history-list"),
            {"activity_group": "check_in_out", "page": 1, "page_size": 5},
        )

        self.assertEqual(response.status_code, status.HTTP_200_OK)
        self.assertEqual(response.data["count"], 3)
        self.assertEqual(
            {record["id"] for record in response.data["results"]},
            {check_in_record.id, stock_increase_record.id, check_out_record.id},
        )

    def test_history_returns_the_current_users_favorite_state(self):
        stock = InventoryStock.objects.create(
            material=self.material,
            sector=self.primary_sector,
            stock_unit=self.stock_unit,
            quantity="2",
            minimum_quantity="1",
        )
        stock.favorite_users.add(self.first_user)
        InventoryChangeRecord.objects.create(
            performed_action=InventoryChangeRecord.ACTION_STOCK_CREATED,
            inventory_stock=stock,
        )

        response = self.client.get(reverse("inventory-history-list"))

        self.assertEqual(response.status_code, status.HTTP_200_OK)
        self.assertTrue(response.data["results"][0]["inventory_stock"]["is_favorite"])


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

        response = self.client.post(reverse("inventory-material-list"), payload, format="json")

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

        response = self.client.post(reverse("inventory-material-list"), payload, format="multipart")

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
