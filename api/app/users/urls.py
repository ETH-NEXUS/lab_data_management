from rest_framework.routers import DefaultRouter
from users.views import UserViewSet

router = DefaultRouter()
router.register("auth/users", UserViewSet, basename="users")
