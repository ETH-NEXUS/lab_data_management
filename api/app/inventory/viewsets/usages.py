from django.db.models import Q
from rest_framework import viewsets
from rest_framework.decorators import action
from rest_framework.response import Response

from ..dynamic_models import MaterialUsage
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
