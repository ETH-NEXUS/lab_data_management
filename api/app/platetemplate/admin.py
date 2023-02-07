from django.contrib import admin

from .models import PlateTemplate, PlateTemplateCategory


@admin.register(PlateTemplateCategory)
class PlateTemplateCategoryAdmin(admin.ModelAdmin):
  pass


@admin.register(PlateTemplate)
class PlateTemplateAdmin(admin.ModelAdmin):
  pass
