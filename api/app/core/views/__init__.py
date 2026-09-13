"""
Views of the core app, split by topic. Import them from `core.views` as before.
"""

from .auth import CsrfCookieView, LoginView, LogoutView
from .experiments import ExperimentViewSet
from .mappings import MappingPreviewView, PlateMappingViewSet
from .plate_info import prefill_plate_info, save_plate_info
from .plates import PlateViewSet
from .projects import ProjectViewSet, add_control_layout
from .reports import (
    download_csv_data,
    download_pdf_report,
    generate_pdf_report,
    list_files,
)
from .system import DocsView, VersionView, refresh
from .thresholds import ThresholdViewSet
from .wells import WellViewSet

__all__ = [
    "CsrfCookieView",
    "DocsView",
    "ExperimentViewSet",
    "LoginView",
    "LogoutView",
    "MappingPreviewView",
    "PlateMappingViewSet",
    "PlateViewSet",
    "ProjectViewSet",
    "ThresholdViewSet",
    "VersionView",
    "WellViewSet",
    "add_control_layout",
    "download_csv_data",
    "download_pdf_report",
    "generate_pdf_report",
    "list_files",
    "prefill_plate_info",
    "refresh",
    "save_plate_info",
]
