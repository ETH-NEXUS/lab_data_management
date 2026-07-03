from django.urls import reverse
from rest_framework import status
from rest_framework.test import APITestCase

from inventory.dynamic_models import InventoryStock, Room, Sector
from inventory.static_models import ItemType, MaterialMaster, MaterialUnit, UnitOfMeasure


class InventoryStockMultiSectorTests(APITestCase):
    """
    Covers the backend stock flow for one item stored in multiple sectors.
    """

    def setUp(self):
        self.item_type = ItemType.objects.create(name="Consumable")
        self.unit_of_measure = UnitOfMeasure.objects.create(name="box")
        self.material = MaterialMaster.objects.create(
            product_name="Reservoir tips",
            item_type=self.item_type,
        )
        self.stock_unit = MaterialUnit.objects.create(
            material=self.material,
            unit=self.unit_of_measure,
            display_name="box",
            base_units_per_unit="1",
            is_stock_unit=True,
        )
        self.room = Room.objects.create(name="C75")
        self.primary_sector = Sector.objects.create(room=self.room, name="3.1")
        self.secondary_sector = Sector.objects.create(room=self.room, name="3.2")
        self.other_room = Room.objects.create(name="C41")
        self.other_room_sector = Sector.objects.create(room=self.other_room, name="1.1")

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
