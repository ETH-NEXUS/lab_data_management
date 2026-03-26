from rest_framework.routers import DefaultRouter

from compoundlib.api_viewsets import (
    CompoundLibraryEndpointViewSet,
    CompoundEndpointViewSet,
)
from core.api_viewsets import (
    PlateEndpointViewSet,
    PlateDimensionEndpointViewSet,
    WellEndpointViewSet,
    ProjectEndpointViewSet,
    ExperimentEndpointViewSet,
    BarcodeSpecificationEndpointViewSet,
    WellCompoundEndpointViewSet,
    WellWithdrawalEndpointViewSet,
    MeasurementFeatureEndpointViewSet,
    MeasurementEndpointViewSet,
    PlateMappingEndpointViewSet,
    ThresholdEndpointViewSet,
)
from platetemplate.api_viewsets import (
    PlateTemplateCategoryEndpointViewSet,
    PlateTemplateEndpointViewSet,
)

router = DefaultRouter()

router.register("plates", PlateEndpointViewSet, basename="plate")
router.register("platedimensions", PlateDimensionEndpointViewSet, basename="platedimension")
router.register("wells", WellEndpointViewSet, basename="well")
router.register("projects", ProjectEndpointViewSet, basename="project")
router.register("experiments", ExperimentEndpointViewSet, basename="experiment")
router.register(
    "barcodespecifications",
    BarcodeSpecificationEndpointViewSet,
    basename="barcodespecification",
)
router.register("wellcompounds", WellCompoundEndpointViewSet, basename="wellcompound")
router.register("wellwithdrawals", WellWithdrawalEndpointViewSet, basename="wellwithdrawal")
router.register(
    "measurementfeatures",
    MeasurementFeatureEndpointViewSet,
    basename="measurementfeature",
)
router.register("measurements", MeasurementEndpointViewSet, basename="measurement")
router.register("platemappings", PlateMappingEndpointViewSet, basename="platemapping")
router.register("thresholds", ThresholdEndpointViewSet, basename="threshold")
router.register(
    "compoundlibraries",
    CompoundLibraryEndpointViewSet,
    basename="compoundlibrary",
)
router.register("compounds", CompoundEndpointViewSet, basename="compound")
router.register(
    "templatecategories",
    PlateTemplateCategoryEndpointViewSet,
    basename="templatecategory",
)
router.register("templates", PlateTemplateEndpointViewSet, basename="template")
