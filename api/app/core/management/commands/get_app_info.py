from django.core.management.base import BaseCommand
from core.models import (
    Project,
    Experiment,
    Plate,
    Well,
    CompoundLibrary,
    Compound,
)


class Command(BaseCommand):
    help = "Draws some basic insights from the database."

    def handle(self, *args, **options):
        num_projects = Project.objects.count()
        num_experiments = Experiment.objects.count()
        num_plates = Plate.objects.count()
        num_wells = Well.objects.count()
        num_libraries = CompoundLibrary.objects.count()
        num_compounds = Compound.objects.count()

        self.stdout.write(self.style.SUCCESS("===== DATA INSIGHTS ====="))
        self.stdout.write(self.style.SUCCESS(f"Number of Projects: {num_projects}"))
        self.stdout.write(
            self.style.SUCCESS(f"Number of Experiments: {num_experiments}")
        )
        self.stdout.write(self.style.SUCCESS(f"Number of Plates: {num_plates}"))
        self.stdout.write(self.style.SUCCESS(f"Number of Wells: {num_wells}"))
        self.stdout.write(
            self.style.SUCCESS(f"Number of Compound Libraries: {num_libraries}")
        )
        self.stdout.write(self.style.SUCCESS(f"Number of Compounds: {num_compounds}"))
        self.stdout.write(self.style.SUCCESS("===== END OF INSIGHTS ====="))
