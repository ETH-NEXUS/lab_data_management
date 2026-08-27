from django.conf import settings
from django.db import migrations, models
import django.db.models.deletion


DASHBOARD_TILES = [
    ("low_stock_items", "Items empty or low in stock", 1, True),
    ("awaiting_check_in", "Ordered items not checked-in yet", 2, True),
    ("device_items", "Items specific for a certain device", 3, True),
    ("favorite_items", "List of favorite items", 4, True),
    ("expired_items", "Recently expired items", 5, True),
    ("archived_items", "Recently archived items", 6, True),
    ("harvest_project_usages", "Items recently linked to a Harvest project", 7, False),
    ("ldm_experiment_usages", "Items recently linked to an LDM experiment", 8, False),
    ("recent_activities", "Recent activities", 9, False),
    ("check_in_out", "Recent check-ins / check-outs", 10, False),
]


def seed_dashboard_tiles(apps, schema_editor):
    InventoryDashboardTile = apps.get_model("inventory", "InventoryDashboardTile")

    for key, name, default_position, is_default_visible in DASHBOARD_TILES:
        InventoryDashboardTile.objects.update_or_create(
            key=key,
            defaults={
                "name": name,
                "default_position": default_position,
                "is_default_visible": is_default_visible,
            },
        )


def remove_dashboard_tiles(apps, schema_editor):
    InventoryDashboardTile = apps.get_model("inventory", "InventoryDashboardTile")
    InventoryDashboardTile.objects.filter(key__in=[tile[0] for tile in DASHBOARD_TILES]).delete()


class Migration(migrations.Migration):
    dependencies = [
        ("inventory", "0010_inventorystock_favorite_users"),
        migrations.swappable_dependency(settings.AUTH_USER_MODEL),
    ]

    operations = [
        migrations.CreateModel(
            name="InventoryDashboardTile",
            fields=[
                ("id", models.BigAutoField(auto_created=True, primary_key=True, serialize=False, verbose_name="ID")),
                ("key", models.CharField(help_text="Stable frontend key, e.g. low_stock_items.", max_length=100, unique=True)),
                ("name", models.CharField(help_text="Readable dashboard tile name, e.g. Items empty or low in stock.", max_length=255)),
                ("default_position", models.PositiveSmallIntegerField(help_text="Position used when this tile is selected for a new user.")),
                ("is_default_visible", models.BooleanField(default=False, help_text="Whether this tile is selected for a new user by default.")),
            ],
            options={
                "verbose_name": "Inventory dashboard tile",
                "verbose_name_plural": "Inventory dashboard tiles",
                "ordering": ["default_position", "id"],
            },
        ),
        migrations.CreateModel(
            name="InventoryDashboardTilePreference",
            fields=[
                ("id", models.BigAutoField(auto_created=True, primary_key=True, serialize=False, verbose_name="ID")),
                ("is_visible", models.BooleanField(default=False, help_text="Whether the tile appears on the user's dashboard.")),
                ("position", models.PositiveSmallIntegerField(help_text="Display order when the tile is visible.")),
                ("created_at", models.DateTimeField(auto_now_add=True)),
                ("updated_at", models.DateTimeField(auto_now=True)),
                ("tile", models.ForeignKey(help_text="Dashboard tile configured by this user.", on_delete=django.db.models.deletion.CASCADE, related_name="user_preferences", to="inventory.inventorydashboardtile")),
                ("user", models.ForeignKey(help_text="User who owns this dashboard configuration.", on_delete=django.db.models.deletion.CASCADE, related_name="inventory_dashboard_tile_preferences", to=settings.AUTH_USER_MODEL)),
            ],
            options={
                "verbose_name": "Inventory dashboard tile preference",
                "verbose_name_plural": "Inventory dashboard tile preferences",
                "ordering": ["user_id", "position", "tile_id"],
                "unique_together": {("user", "tile")},
            },
        ),
        migrations.RunPython(seed_dashboard_tiles, remove_dashboard_tiles),
    ]
