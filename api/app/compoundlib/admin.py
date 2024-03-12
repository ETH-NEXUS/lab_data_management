from django.contrib import admin
from django.utils.html import mark_safe
from .models import Compound, CompoundLibrary


@admin.register(CompoundLibrary)
class CompoundLibraryAdmin(admin.ModelAdmin):
    list_display = ("name", "file_name")
    search_fields = ("name", "file_name")


@admin.register(Compound)
class CompoundAdmin(admin.ModelAdmin):
    list_display = (
        "get_structure",
        "structure",
        "get_data",
    )
    search_fields = ("name", "structure")
    list_filter = ("name",)

    def get_structure(self, compound: Compound):
        if not compound.structure:
            return None
        return mark_safe(
            '<img src="{}" width="150" height="150" />'.format(compound.structure_image)
        )

    get_structure.short_description = "structure"

    def get_data(self, compound: Compound):
        out = "<table>"
        n_items = 0
        if compound.data and compound.data.items():
            for k, v in compound.data.items():
                out += f"<tr><td>{k}</td><td>{v}</td></tr>"
                if n_items >= 5:
                    break
                else:
                    n_items += 1
            out += "</table>"
            return mark_safe(out)

    get_data.short_description = "data"
