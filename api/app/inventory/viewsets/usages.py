from django.db.models import Q
from django.db import transaction
from rest_framework import viewsets
from rest_framework.decorators import action
from rest_framework.response import Response

from ..dynamic_models import MaterialUsage
from ..history_models import InventoryChangeRecord
from ..history_utils import record_inventory_action
from ..serializers.dynamic_models_serializers import MaterialUsageDetailSerializer, MaterialUsageListSerializer


class MaterialUsageViewSet(viewsets.ModelViewSet):
    queryset = (
        MaterialUsage.objects.all()
        .select_related(
            "project",
            "experiment",
            "inventory_stock",
            "inventory_stock__material",
            "inventory_stock__material__brand",
            "inventory_stock__material__manufacturer",
            "inventory_stock__material__vendor",
            "inventory_stock__material__item_type",
            "inventory_stock__material__device_type",
            "inventory_stock__sector",
            "inventory_stock__sector__room",
            "inventory_stock__stock_unit",
            "inventory_stock__stock_unit__unit",
            "usage_unit",
            "usage_unit__unit",
        )
        .prefetch_related(
            "inventory_stock__material__attributes",
            "inventory_stock__material__units__unit",
        )
        .order_by("-used_at")
    )

    def get_serializer_class(self):
        if self.action == "list":
            return MaterialUsageListSerializer
        return MaterialUsageDetailSerializer

    def _build_usage_history_note(self, material_usage):
        """
        Builds one small readable usage context string.

        Returned data examples:
        - "PCR Plates used for Project A"
        - "Tips used for Project A / Experiment B"
        """
        material_name = "Unknown material"
        if material_usage.inventory_stock_id and material_usage.inventory_stock.material_id:
            material_name = material_usage.inventory_stock.material.product_name

        project_name = material_usage.project.name if material_usage.project_id else "Unknown project"

        if material_usage.experiment_id:
            return f"{material_name} used for {project_name} / {material_usage.experiment.name}"

        return f"{material_name} used for {project_name}"

    def perform_create(self, serializer):
        with transaction.atomic():
            material_usage = serializer.save()
            performed_by = self.request.user if self.request.user.is_authenticated else None

            record_inventory_action(
                performed_action=InventoryChangeRecord.ACTION_USAGE_CREATED,
                performed_by=performed_by,
                inventory_stock=material_usage.inventory_stock,
                material_usage=material_usage,
                project=material_usage.project,
                experiment=material_usage.experiment,
                quantity_delta=material_usage.quantity_used,
                quantity_unit=material_usage.usage_unit,
                notes=f"Created usage entry for {self._build_usage_history_note(material_usage)}.",
            )

    def perform_update(self, serializer):
        previous_usage = self.get_object()
        previous_quantity = previous_usage.quantity_used
        previous_note = self._build_usage_history_note(previous_usage)

        with transaction.atomic():
            material_usage = serializer.save()
            performed_by = self.request.user if self.request.user.is_authenticated else None

            record_inventory_action(
                performed_action=InventoryChangeRecord.ACTION_USAGE_UPDATED,
                performed_by=performed_by,
                inventory_stock=material_usage.inventory_stock,
                material_usage=material_usage,
                project=material_usage.project,
                experiment=material_usage.experiment,
                quantity_delta=material_usage.quantity_used - previous_quantity,
                quantity_unit=material_usage.usage_unit,
                notes=(
                    f"Updated usage entry from {previous_note} "
                    f"to {self._build_usage_history_note(material_usage)}."
                ),
            )

    def perform_destroy(self, instance):
        performed_by = self.request.user if self.request.user.is_authenticated else None
        usage_note = self._build_usage_history_note(instance)

        with transaction.atomic():
            record_inventory_action(
                performed_action=InventoryChangeRecord.ACTION_USAGE_DELETED,
                performed_by=performed_by,
                inventory_stock=instance.inventory_stock,
                material_usage=instance,
                project=instance.project,
                experiment=instance.experiment,
                quantity_delta=instance.quantity_used,
                quantity_unit=instance.usage_unit,
                notes=f"Deleted usage entry for {usage_note}.",
            )
            instance.delete()

    def get_queryset(self):
        queryset = super().get_queryset()

        project_id = self.request.query_params.get("project")
        experiment_id = self.request.query_params.get("experiment")
        material_id = self.request.query_params.get("material")
        search = self.request.query_params.get("search")

        if project_id:
            queryset = queryset.filter(project_id=project_id)

        if experiment_id:
            queryset = queryset.filter(experiment_id=experiment_id)

        if material_id:
            queryset = queryset.filter(inventory_stock__material_id=material_id)

        if search:
            queryset = queryset.filter(
                Q(project__name__icontains=search)
                | Q(inventory_stock__material__product_name__icontains=search)
                | Q(notes__icontains=search)
            )

        return queryset

    @action(detail=False, methods=["get"])
    def by_project(self, request):
        project_id = request.query_params.get("project")
        queryset = self.get_queryset()

        if project_id:
            queryset = queryset.filter(project_id=project_id)

        serializer = MaterialUsageListSerializer(queryset, many=True, context={"request": request})
        return Response(serializer.data)

    @action(detail=False, methods=["get"])
    def recent_project(self, request):
        """
        Returns the five latest material usages linked to a Harvest project.

        Returned item example:
        - {"id": 42, "project": {"name": "Project A"}, "used_at": "2026-08-27T10:00:00Z"}
        """
        queryset = self.get_queryset().filter(project__harvest_id__isnull=False)[:5]
        serializer = MaterialUsageListSerializer(queryset, many=True, context={"request": request})
        return Response(serializer.data)

    @action(detail=False, methods=["get"])
    def recent_experiment(self, request):
        """
        Returns the five latest material usages linked to an LDM experiment.

        Returned item example:
        - {"id": 42, "experiment": {"name": "Experiment A"}, "used_at": "2026-08-27T10:00:00Z"}
        """
        queryset = self.get_queryset().filter(experiment__isnull=False)[:5]
        serializer = MaterialUsageListSerializer(queryset, many=True, context={"request": request})
        return Response(serializer.data)
