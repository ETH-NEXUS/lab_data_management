"""
Models of the core app, split by topic. Import them from `core.models` as before.

Every model has to be imported here: Django only registers the models it finds
when it loads `core.models`.
"""

from .mappings import PlateMapping
from .materialized_views import (
    DictField,
    ExperimentDetail,
    MaterializedViewModel,
    PlateDetail,
    WellDetail,
)
from .measurements import Measurement, MeasurementAssignment, MeasurementFeature
from .plate_info import PlateInfo
from .plates import Location, MappingError, Plate, PlateDimension
from .projects import BarcodeSpecification, Experiment, Project
from .thresholds import Threshold
from .wells import Sample, Well, WellCompound, WellType, WellWithdrawal

__all__ = [
    "BarcodeSpecification",
    "DictField",
    "Experiment",
    "ExperimentDetail",
    "Location",
    "MappingError",
    "MaterializedViewModel",
    "Measurement",
    "MeasurementAssignment",
    "MeasurementFeature",
    "Plate",
    "PlateDetail",
    "PlateDimension",
    "PlateInfo",
    "PlateMapping",
    "Project",
    "Sample",
    "Threshold",
    "Well",
    "WellCompound",
    "WellDetail",
    "WellType",
    "WellWithdrawal",
]
