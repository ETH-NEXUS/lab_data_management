"""
Projects and adding control layouts to them.
"""

import traceback
from rest_framework import status
from rest_framework.decorators import api_view
from rest_framework.response import Response
from helpers.logger import logger
from ..models import (
    Plate,
    Project,
    PlateDetail,
    WellDetail,
)


@api_view(["POST"])
def add_control_layout(request):
    """Add a selected control layout to a project"""

    data = request.data
    project_id = data.get("project_id")
    barcode_new = data.get("barcode_new")
    barcode_old = data.get("barcode_old")

    try:
        project = Project.objects.get(id=project_id)
        plate_old = Plate.objects.get(barcode=barcode_old)

        plate_new = Plate(
            barcode=barcode_new,
            dimension=plate_old.dimension,
            library=plate_old.library,
            template=plate_old.template,
            is_control_plate=True,
            project=project,
        )
        plate_new.save()
        plate_old.copy(plate_new, map_type=True)
        PlateDetail.refresh(concurrently=True)
        WellDetail.refresh(concurrently=True)
        return Response(status=status.HTTP_200_OK)

    except Project.DoesNotExist:
        logger.error(f"Project not found: {project_id}")
        return Response(
            {"error": "Project not found"}, status=status.HTTP_404_NOT_FOUND
        )
    except Plate.DoesNotExist:
        logger.error(f"Old plate not found: {barcode_old}")
        return Response(
            {"error": "Old plate not found"}, status=status.HTTP_404_NOT_FOUND
        )
    except Exception as e:
        traceback.print_exc()
        logger.error(f"Error adding control layout: {e}")
        return Response({"error": str(e)}, status=status.HTTP_500_INTERNAL_SERVER_ERROR)
