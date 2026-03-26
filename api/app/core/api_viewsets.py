from django_filters.rest_framework import DjangoFilterBackend
from rest_framework import serializers, viewsets
from rest_framework.filters import OrderingFilter

from .models import (
    PlateDimension,
    Project,
    WellCompound,
    WellWithdrawal,
    MeasurementFeature,
    Measurement,
    Threshold,
    BarcodeSpecification,
)
from .serializers import (
    ProjectSerializer,
    ThresholdSerializer,
    BarcodeSpecificationSerializer,
)
from .views import (
    PlateViewSet,
    WellViewSet,
    ExperimentViewSet,
    PlateMappingViewSet,
)


class PlateEndpointViewSet(PlateViewSet):
    filter_backends = (DjangoFilterBackend, OrderingFilter)
    filterset_fields = (
        "barcode",
        "library",
        "experiment",
        "is_control_plate",
        "use_as_template_to_select",
    )
    ordering_fields = ("barcode",)


class PlateDimensionSerializer(serializers.ModelSerializer):
    class Meta:
        model = PlateDimension
        fields = "__all__"


class PlateDimensionEndpointViewSet(viewsets.ModelViewSet):
    serializer_class = PlateDimensionSerializer
    queryset = PlateDimension.objects.all()


class WellEndpointViewSet(WellViewSet):
    filter_backends = (DjangoFilterBackend,)
    filterset_fields = ("plate__barcode",)


class ProjectEndpointViewSet(viewsets.ModelViewSet):
    serializer_class = ProjectSerializer
    queryset = Project.objects.all()


class ExperimentEndpointViewSet(ExperimentViewSet):
    pass


class BarcodeSpecificationEndpointViewSet(viewsets.ModelViewSet):
    serializer_class = BarcodeSpecificationSerializer
    queryset = BarcodeSpecification.objects.all()


class WellCompoundSerializer(serializers.ModelSerializer):
    class Meta:
        model = WellCompound
        fields = "__all__"


class WellCompoundEndpointViewSet(viewsets.ModelViewSet):
    serializer_class = WellCompoundSerializer
    queryset = WellCompound.objects.all()
    filter_backends = (DjangoFilterBackend,)
    filterset_fields = ("well__plate__barcode",)


class WellWithdrawalSerializer(serializers.ModelSerializer):
    class Meta:
        model = WellWithdrawal
        fields = "__all__"


class WellWithdrawalEndpointViewSet(viewsets.ModelViewSet):
    serializer_class = WellWithdrawalSerializer
    queryset = WellWithdrawal.objects.all()


class MeasurementFeatureSerializer(serializers.ModelSerializer):
    class Meta:
        model = MeasurementFeature
        fields = "__all__"


class MeasurementFeatureEndpointViewSet(viewsets.ModelViewSet):
    serializer_class = MeasurementFeatureSerializer
    queryset = MeasurementFeature.objects.all()


class MeasurementSerializer(serializers.ModelSerializer):
    class Meta:
        model = Measurement
        fields = "__all__"


class MeasurementEndpointViewSet(viewsets.ModelViewSet):
    serializer_class = MeasurementSerializer
    queryset = Measurement.objects.all()


class PlateMappingEndpointViewSet(PlateMappingViewSet):
    pass


class ThresholdEndpointViewSet(viewsets.ModelViewSet):
    serializer_class = ThresholdSerializer
    queryset = Threshold.objects.all()
