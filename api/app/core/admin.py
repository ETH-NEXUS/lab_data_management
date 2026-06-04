from django import forms
from django.contrib import admin, messages
from django.contrib.admin import helpers
from django.template.response import TemplateResponse
from django.utils.html import mark_safe

from .models import (
    Plate,
    PlateDimension,
    Measurement,
    MeasurementFeature,
    Sample,
    Well,
    Location,
    Project,
    Experiment,
    BarcodeSpecification,
    PlateMapping,
    MeasurementAssignment,
    WellType,
    Threshold,
    PlateInfo,
)
from .utils import copy_library_plates


class CopySelectedPlatesForm(forms.Form):
    target_volume = forms.FloatField(label="Target volume (nL)")

    def clean_target_volume(self):
        target_volume = self.cleaned_data["target_volume"]
        if target_volume <= 0:
            raise forms.ValidationError("Target volume must be greater than 0.")
        return target_volume


@admin.action(description="Copy selected plates")
def copy_selected_plates(modeladmin, request, queryset):
    if "apply" in request.POST:
        form = CopySelectedPlatesForm(request.POST)
        if form.is_valid():
            target_volume = form.cleaned_data["target_volume"]
            try:
                copied_plates = copy_library_plates(queryset, target_volume)
            except ValueError as error:
                form.add_error(None, str(error))
            else:
                modeladmin.message_user(
                    request,
                    f"Copied {len(copied_plates)} plate(s) with {target_volume} nL per well.",
                    level=messages.SUCCESS,
                )
                return None
    else:
        form = CopySelectedPlatesForm()

    context = {
        **modeladmin.admin_site.each_context(request),
        "action_checkbox_name": helpers.ACTION_CHECKBOX_NAME,
        "form": form,
        "opts": modeladmin.model._meta,
        "plates": queryset,
        "select_across": request.POST.get("select_across", "0"),
        "title": "Copy selected plates",
    }
    return TemplateResponse(request, "admin/core/copy_selected_plates.html", context)


@admin.register(Plate)
class PlateAdmin(admin.ModelAdmin):
    list_display = (
        "barcode",
        "is_control_plate",
        "dimension",
        "library",
        "experiment",
        "use_as_template_to_select",
    )
    search_fields = ("barcode",)

    list_filter = (
        "dimension",
        "is_control_plate",
        "library",
        "experiment",
        "project",
        "template",
        "use_as_template_to_select",
    )
    actions = (copy_selected_plates,)

    def get_wells(self, plate: Plate):
        out = "<table><tr>"
        cols = plate.dimension.cols
        for pos, well in enumerate(plate.wells.order_by("position").all()):
            if pos % cols == 0:
                out += "</tr><tr>"
            if well.compounds:
                out += f"<td><a href='#' title='{', '.join([c.identifier for c in well.compounds])}'>{well.hr_position}</a></td>"
            else:
                out += f"<td><a href='#' title='EMPTY'>{well.hr_position}<br/></a></td>"
        out += "</tr></table>"
        return mark_safe(out)

    get_wells.short_description = "wells"


@admin.register(PlateDimension)
class PlateDimensionAdmin(admin.ModelAdmin):
    list_display = ("name", "cols", "rows")
    search_fields = ("name",)


@admin.register(Well)
class WellAdmin(admin.ModelAdmin):
    list_display = (
        "__str__",
        "plate",
        "position",
        "hr_position",
        "sample",
        "amount",
        "type",
        "get_compounds",
        "get_measurements",
    )
    search_fields = ("plate__barcode", "compounds__identifier")
    list_filter = ("plate__barcode",)
    autocomplete_fields = ("plate", "compounds", "sample")

    def get_measurements(self, well: Well):
        return ", ".join([str(m) for m in well.measurements.all()])

    get_measurements.short_description = "measurements"

    def get_compounds(self, well: Well):
        return ", ".join([str(c) for c in well.compounds.all()])

    get_compounds.short_description = "compounds"


@admin.register(Measurement)
class MeasurementAdmin(admin.ModelAdmin):
    raw_id_fields = ("well",)
    list_display = (
        "value",
        "label",
        "identifier",
        "measured_at",
    )
    list_filter = (
        "measured_at",
        "label",
        ("measurement_assignment__plate", admin.RelatedOnlyFieldListFilter),
    )
    list_per_page = 1000


@admin.register(MeasurementFeature)
class MeasurementFeatureAdmin(admin.ModelAdmin):
    list_display = ("name", "abbrev", "unit")


@admin.register(Sample)
class SampleAdmin(admin.ModelAdmin):
    list_display = ("name",)
    search_fields = ("name",)


@admin.register(Location)
class LocationAdmin(admin.ModelAdmin):
    list_display = ("name",)


@admin.register(Project)
class ProjectAdmin(admin.ModelAdmin):
    list_display = ("name",)
    search_fields = ("name",)



@admin.register(Experiment)
class ExperimentAdmin(admin.ModelAdmin):
    list_display = ("name",)
    search_fields = ("name",)



@admin.register(BarcodeSpecification)
class BarcodeSpecificationAdmin(admin.ModelAdmin):
    list_display = (
        "prefix",
        "number_of_plates",
    )


@admin.register(PlateMapping)
class PlateMappingAdmin(admin.ModelAdmin):
    list_display = (
        "source_plate",
        "target_plate",
        "mapping_file",
    )
    search_fields = (
        "source_plate",
        "target_plate",
        "mapping_file",
    )


@admin.register(MeasurementAssignment)
class MeasurementAssignmentAdmin(admin.ModelAdmin):
    list_display = (
        "filename",
        "status",
        "created_at",
    )


@admin.register(WellType)
class WellTypeAdmin(admin.ModelAdmin):
    list_display = (
        "name",
        "description",
    )


@admin.register(Threshold)
class ThresholdAdmin(admin.ModelAdmin):
    list_display = (
        "amount",
        "dmso",
    )
    search_fields = (
        "amount",
        "dmso",
    )


@admin.register(PlateInfo)
class PlateInfoAdmin(admin.ModelAdmin):
    list_display = (
        "plate",
        "lib_plate_barcode",
        "label",
        "replicate",
        "measurement_time",
        "cell_type",
        "condition",
    )
    search_fields = (
        "plate",
        "lib_plate_barcode",
        "label",
        "replicate",
        "measurement_time",
        "cell_type",
        "condition",
    )
    list_filter = (
        "plate__barcode",
        "lib_plate_barcode",
        "label",
        "replicate",
        "measurement_time",
        "cell_type",
        "condition",
    )
