from rest_framework import serializers
from .models import PlateTemplateCategory, PlateTemplate
from core.serializers import SimplePlateSerializer


class PlateTemplateSerializer(serializers.ModelSerializer):
  plate = SimplePlateSerializer()

  class Meta:
    model = PlateTemplate
    fields = ('id', 'name', 'category', 'plate')


class PlateTemplateCategorySerializer(serializers.ModelSerializer):
  templates = PlateTemplateSerializer(many=True, read_only=True)

  class Meta:
    model = PlateTemplateCategory
    fields = '__all__'
