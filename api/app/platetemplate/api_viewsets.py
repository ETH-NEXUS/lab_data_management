from django_filters.rest_framework import DjangoFilterBackend
from rest_framework import serializers, viewsets

from .models import PlateTemplateCategory, PlateTemplate
from .serializers import PlateTemplateCategorySerializer


class PlateTemplateCategoryEndpointViewSet(viewsets.ModelViewSet):
    serializer_class = PlateTemplateCategorySerializer
    queryset = PlateTemplateCategory.objects.all()
    filter_backends = (DjangoFilterBackend,)
    filterset_fields = ("name",)


class PlateTemplateSerializer(serializers.ModelSerializer):
    class Meta:
        model = PlateTemplate
        fields = "__all__"


class PlateTemplateEndpointViewSet(viewsets.ModelViewSet):
    serializer_class = PlateTemplateSerializer
    queryset = PlateTemplate.objects.all()
    filter_backends = (DjangoFilterBackend,)
    filterset_fields = ("name",)
