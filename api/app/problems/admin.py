from django.contrib import admin

from problems.models import Problem

#
#
@admin.register(Problem)
class ProblemAdmin(admin.ModelAdmin):
    list_display = (
        "well",
        "plate",
        "library",
        "type",
        "status",
        "details",
        "show",
    )

    list_filter = (
        "type",
        "show",
    )
