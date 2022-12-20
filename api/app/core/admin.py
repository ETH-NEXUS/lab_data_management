from django.contrib import admin
from django.utils.html import mark_safe

from .models import (Plate, PlateDimension, Measurement, Sample, Well, Location, Project)


@admin.register(Plate)
class PlateAdmin(admin.ModelAdmin):
    list_display = ('barcode', 'dimension', 'library', 'experiment')
    search_fields = ('barcode',)
    readonly_fields = ('get_wells',)

    def get_wells(self, plate: Plate):
        out = '<table><tr>'
        cols = plate.dimension.cols
        for pos, well in enumerate(plate.wells.order_by('position').all()):
            if pos % cols == 0:
                out += '</tr><tr>'
            if well.compounds:
                out += f"<td><a href='#' title='{', '.join([c.identifier for c in well.compounds])}'>{well.hr_position}</a></td>"
            else:
                out += f"<td><a href='#' title='EMPTY'>{well.hr_position}<br/></a></td>"
        out += '</tr></table>'
        return mark_safe(out)

    get_wells.short_description = 'wells'


@admin.register(PlateDimension)
class PlateDimensionAdmin(admin.ModelAdmin):
    list_display = ('name', 'cols', 'rows')
    search_fields = ('name',)


@admin.register(Well)
class WellAdmin(admin.ModelAdmin):
    list_display = ('__str__', 'plate', 'position', 'hr_position', 'sample', 'amount', 'get_measurements')
    search_fields = ('plate__barcode', 'compounds__identifier')
    list_filter = ('plate__barcode', )
    autocomplete_fields = ('plate', 'compounds', 'sample')

    def get_measurements(self, well: Well):
        return ', '.join([str(m) for m in well.measurements.all()])
    get_measurements.short_description = 'measurements'


@admin.register(Measurement)
class MeasurementAdmin(admin.ModelAdmin):
    list_display = ('name', 'abbrev', 'unit', 'value')


@admin.register(Sample)
class SampleAdmin(admin.ModelAdmin):
    list_display = ('name',)
    search_fields = ('name',)


@admin.register(Location)
class LocationAdmin(admin.ModelAdmin):
    list_display = ('name',)


@admin.register(Project)
class ProjectAdmin(admin.ModelAdmin):
    list_display = ('name',)
