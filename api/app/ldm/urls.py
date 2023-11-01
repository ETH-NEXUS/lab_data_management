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
from core.views import (
    MappingPreviewView,
    VersionView,
    DocsView,
    CsrfCookieView,
    LoginView,
    LogoutView,
    generate_pdf_report,
    list_files,
    download_pdf_report,
)
from django.conf import settings
from jupyter.views import JupyterProxyView
from users.urls import router as user_router
from management.views import (
    directory_content,
    run_command,
    long_polling,
    delete_file,
    download_file,
    upload_file,
    get_file_content,
)

urlpatterns = [
    path("admin/", admin.site.urls),
    path("api/auth/cookie/", CsrfCookieView.as_view(), name="auth-cookie"),
    path("api/auth/login/", LoginView.as_view(), name="login"),
    path("api/auth/logout/", LogoutView.as_view(), name="logout"),
    path("api/", include(router.urls)),
    path("api/", include(user_router.urls)),
    path("api/mapping_preview/", MappingPreviewView.as_view()),
    path("api/version/", VersionView.as_view()),
    re_path(r"^docs/(?P<uri>.*)$", DocsView.as_view(), name="docs"),
    path("api/harvest/", include("harvest.urls")),
    re_path(
        "(?P<path>notebook/.*)$",
        JupyterProxyView.as_view(),
        name="notebook",
    ),
    path("api/directory_content/", directory_content, name="directory_content"),
    path("api/run_command/", run_command, name="run_command"),
    path("api/long_polling/<str:room_name>/", long_polling, name="long_polling"),
    path("api/delete_file/", delete_file, name="delete_file"),
    path("api/download_file/", download_file, name="download_file"),
    path("api/upload_file/", upload_file, name="upload_file"),
    path("api/get_file_content/", get_file_content, name="get_file_content"),
    path("api/generate_pdf_report/", generate_pdf_report, name="generate_pdf_report"),
    path("api/list_files/", list_files, name="list_files"),
    path("api/download_pdf_report/", download_pdf_report, name="download_pdf_report"),
]

if not settings.DISABLE_BROWSABLE_API and not settings.DISABLE_AUTH:
    urlpatterns += [path("api-auth/", include("rest_framework.urls"))]

# if settings.DEBUG:
#     urlpatterns += static(settings.STATIC_URL, document_root=settings.STATIC_ROOT)
#     urlpatterns += static(settings.MEDIA_URL, document_root=settings.MEDIA_ROOT)
