"""
Experiments.
"""

from django.core.exceptions import ValidationError
from django.db import IntegrityError
from django.http import Http404
from django.utils.translation import gettext as _
from rest_framework import viewsets, status
from rest_framework.decorators import action
from rest_framework.response import Response
from helpers.logger import logger
from ..models import (
    Plate,
    Experiment,
    BarcodeSpecification,
    PlateDimension,
    PlateDetail,
    WellDetail,
)
from ..serializers import (
    ExperimentSerializer,
    ExperimentDetail,
)


# custom pagination class with 1000 items per page


class ExperimentViewSet(viewsets.ModelViewSet):
    serializer_class = ExperimentSerializer
    queryset = Experiment.objects.all()
    pagination_class = None

    @action(detail=False, methods=["post"])
    def move_plates(self, request):
        """Moves selected plates from their original experiment to another experiment"""
        data = request.data
        try:
            experiment = Experiment.objects.get(name=data["experiment"])
            plate_barcodes = data["plate_barcodes"]
            for plate_barcode in plate_barcodes:
                plate = Plate.objects.get(barcode=plate_barcode)
                plate.experiment = experiment
                plate.save()
                logger.info(
                    f"Moved plate {plate_barcode} to experiment {experiment.id}"
                )
            PlateDetail.refresh(concurrently=True)
            WellDetail.refresh(concurrently=True)
            ExperimentDetail.refresh(concurrently=True)
            return Response(status=status.HTTP_200_OK)
        except Experiment.DoesNotExist:
            logger.error(f"Experiment {data['experiment']} does not exist")
            return Response(status=status.HTTP_404_NOT_FOUND)

    @action(detail=False, methods=["post"])
    def barcodes(self, request):
        """Saves a barcode specifications for an experiment"""
        data = request.data
        experiment = Experiment.objects.get(pk=data["experiment_id"])
        barcode_specification = BarcodeSpecification(
            prefix=data["prefix"],
            number_of_plates=data["number_of_plates"],
            sides=data["sides"],
            experiment=experiment,
        )
        barcode_specification.save()
        return Response(status=status.HTTP_200_OK)

    @action(detail=False, methods=["post"])
    def bulk_apply_template(self, request, pk=None):
        """Applies a template plate to all the plate of the experiment"""
        experiment = Experiment.objects.get(pk=request.data.get("experiment_id"))
        template_plate_id = request.data.get("template")
        if template_plate_id:
            template_plate = Plate.objects.get(pk=template_plate_id)
            plates = Plate.objects.filter(experiment=experiment)
            for plate in plates:
                logger.debug(f"Applying template to plate {plate.barcode}")
                plate.apply_template(template_plate)
            return Response(status.HTTP_200_OK)
        else:
            raise Http404(
                _("Parameters 'template' and 'experiment_id' are " "required.")
            )

    @action(detail=False, methods=["post"])
    def bulk_add_plates(self, request):
        try:
            experiment_id = int(request.data["experiment_id"])
            barcode_specification_id = int(request.data["barcode_specification_id"])
            plate_dimension_id = int(request.data["plate_dimension_id"])
            experiment = Experiment.objects.get(pk=experiment_id)
            barcode_specification = BarcodeSpecification.objects.get(
                pk=barcode_specification_id
            )
            plate_dimension = PlateDimension.objects.get(pk=plate_dimension_id)
            number_of_plates = barcode_specification.number_of_plates
        except Experiment.DoesNotExist:
            return Response(
                {"error": _("Experiment not found")},
                status=status.HTTP_404_NOT_FOUND,
            )
        except BarcodeSpecification.DoesNotExist:
            return Response(
                {"error": _("Barcode specification not found")},
                status=status.HTTP_404_NOT_FOUND,
            )
        except PlateDimension.DoesNotExist:
            return Response(
                {"error": _("Plate dimension not found")},
                status=status.HTTP_404_NOT_FOUND,
            )

        plates = []
        for i in range(number_of_plates):
            plate = Plate(
                barcode=barcode_specification.get_barcode_by_number(i + 1),
                experiment=experiment,
                dimension=plate_dimension,
            )
            plates.append(plate)

        try:
            Plate.objects.bulk_create(plates)
        except IntegrityError:
            return Response(
                {
                    "error": _(
                        "Could not save plates to database. Probably you have already added the plates "
                        "with these barcode prefix to the current experiment"
                    )
                },
                status=status.HTTP_500_INTERNAL_SERVER_ERROR,
            )
        except ValidationError:
            return Response(
                {"error": _("Invalid plate data")},
                status=status.HTTP_400_BAD_REQUEST,
            )

        return Response(
            {"success": _("Plates added successfully")}, status=status.HTTP_200_OK
        )
