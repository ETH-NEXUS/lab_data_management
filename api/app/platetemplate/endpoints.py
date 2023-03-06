from drf_auto_endpoint.router import register

from .models import PlateTemplate, PlateTemplateCategory
from .serializers import PlateTemplateCategorySerializer
from core.endpoints import DefaultEndpoint


@register
class PlateTemplateCategoryEndpoint(DefaultEndpoint):
    url = "templatecategories"
    model = PlateTemplateCategory
    base_serializer = PlateTemplateCategorySerializer
    filter_fields = ("name",)


@register
class PlateTemplateEndpoint(DefaultEndpoint):
    url = "templates"
    model = PlateTemplate
    filter_fields = ("name",)
