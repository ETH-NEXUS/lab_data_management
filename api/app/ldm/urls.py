"""app URL Configuration

The `urlpatterns` list routes URLs to views. For more information please see:
    https://docs.djangoproject.com/en/4.1/topics/http/urls/
Examples:
Function views
    1. Add an import:  from my_app import views
    2. Add a URL to urlpatterns:  path('', views.home, name='home')
Class-based views
    1. Add an import:  from other_app.views import Home
    2. Add a URL to urlpatterns:  path('', Home.as_view(), name='home')
Including another URLconf
    1. Import the include() function: from django.urls import include, path
    2. Add a URL to urlpatterns:  path('blog/', include('blog.urls'))
"""
from django.contrib import admin
from django.urls import path, include, re_path
from drf_auto_endpoint.router import router
from django.conf import settings
from core.views import MappingPreviewView, VersionView, DocsView
from rest_framework_simplejwt.views import (
    TokenObtainPairView,
    TokenRefreshView,
)
from users.urls import router as user_router

urlpatterns = [
    path("admin/", admin.site.urls),
    path(
        "api/auth/token/",
        TokenObtainPairView.as_view(),
        name="token_obtain_pair",
    ),
    path(
        "api/auth/token/refresh/",
        TokenRefreshView.as_view(),
        name="token_refresh",
    ),
    path("api/", include(router.urls)),
    path("api/", include(user_router.urls)),
    path("api/mapping_preview/", MappingPreviewView.as_view()),
    path("api/version/", VersionView.as_view()),
    re_path(r'^docs/(?P<uri>.*)$', DocsView.as_view(), name='docs'),
]

if not settings.DISABLE_BROWSABLE_API and not settings.DISABLE_AUTH:
    urlpatterns += [path("api-auth/", include("rest_framework.urls"))]

# if settings.DEBUG:
#     urlpatterns += static(settings.STATIC_URL, document_root=settings.STATIC_ROOT)
#     urlpatterns += static(settings.MEDIA_URL, document_root=settings.MEDIA_ROOT)
