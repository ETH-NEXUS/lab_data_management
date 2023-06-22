"""
Django settings for app project.

Generated by 'django-admin startproject' using Django 4.1.3.

For more information on this file, see
https://docs.djangoproject.com/en/4.1/topics/settings/

For the full list of settings and their values, see
https://docs.djangoproject.com/en/4.1/ref/settings/
"""

import ldap
from pathlib import Path
from os import environ
from django_auth_ldap.config import LDAPSearch

from corsheaders.defaults import default_headers
import logging

# Build paths inside the project like this: BASE_DIR / 'subdir'.
BASE_DIR = Path(__file__).resolve().parent.parent


# Quick-start development settings - unsuitable for production
# See https://docs.djangoproject.com/en/4.1/howto/deployment/checklist/

# SECURITY WARNING: keep the secret key used in production secret!
SECRET_KEY = "django-insecure-ba*bz(ouxsm(gwf_^r=+2*-0h0yh^5d4dz(t6$k+x2@drh(kv_"

# SECURITY WARNING: don't run with debug turned on in production!
DEBUG = (environ.get("DJANGO_DEBUG", "False")) == "True"
LOG_LEVEL = environ.get("DJANGO_LOG_LEVEL", "INFO")
LOG_SQL = False
LOG_LDAP = False
LOG_WEBSOCKET = False
LOG_REVPROXY = False
LOG_DAPHNE = False
LOG_URLLIB = False

# Application definition

INSTALLED_APPS = [
    "daphne",
    "jupyter",
    "django.contrib.admin",
    "django.contrib.auth",
    "django.contrib.contenttypes",
    "django.contrib.sessions",
    "django.contrib.messages",
    "django.contrib.staticfiles",
    "django_filters",
    "django_extensions",
    "corsheaders",
    "rest_framework",
    "revproxy",
    "drf_auto_endpoint",
    "pg_ext",
    "core",
    "importer",
    "compoundlib",
    "platetemplate",
    "harvest",
]

MIDDLEWARE = [
    "django.middleware.security.SecurityMiddleware",
    "django.contrib.sessions.middleware.SessionMiddleware",
    "corsheaders.middleware.CorsMiddleware",
    "django.middleware.common.CommonMiddleware",
    "django.middleware.csrf.CsrfViewMiddleware",
    "django.contrib.auth.middleware.AuthenticationMiddleware",
    "django.contrib.messages.middleware.MessageMiddleware",
    "django.middleware.clickjacking.XFrameOptionsMiddleware",
]

ROOT_URLCONF = "ldm.urls"

TEMPLATES = [
    {
        "BACKEND": "django.template.backends.django.DjangoTemplates",
        "DIRS": [],
        "APP_DIRS": True,
        "OPTIONS": {
            "context_processors": [
                "django.template.context_processors.debug",
                "django.template.context_processors.request",
                "django.contrib.auth.context_processors.auth",
                "django.contrib.messages.context_processors.messages",
            ],
        },
    },
]

ASGI_APPLICATION = "ldm.asgi.application"


# Database
# https://docs.djangoproject.com/en/4.1/ref/settings/#databases

DATABASES = {
    "default": {
        "ENGINE": "django.db.backends.postgresql_psycopg2",
        "HOST": environ.get("POSTGRES_HOST"),
        "PORT": environ.get("POSTGRES_PORT"),
        "NAME": environ.get("POSTGRES_DB"),
        "USER": environ.get("POSTGRES_USER"),
        "PASSWORD": environ.get("POSTGRES_PASSWORD"),
    }
}

# LDAP Authentication
# https://django-auth-ldap.readthedocs.io/en/latest/authentication.html

AUTHENTICATION_BACKENDS = [
    "django_auth_ldap.backend.LDAPBackend",
    "django.contrib.auth.backends.ModelBackend",
]

AUTH_LDAP_GLOBAL_OPTIONS = {ldap.OPT_X_TLS_REQUIRE_CERT: ldap.OPT_X_TLS_NEVER}
AUTH_LDAP_SERVER_URI = environ.get(
    "LDAP_SERVER_URI",
    "ldaps://ldaps-rz-1.ethz.ch,ldaps://ldaps-rz-2.ethz.ch,ldaps://ldaps-hit-2.ethz.ch,ldaps://ldaps-hit-1.ethz.ch",
)
AUTH_LDAP_BIND_DN = environ.get(
    "LDAP_BIND_DN",
    "cn=nexus-tpreports_proxy,ou=admins,ou=nethz,ou=id,ou=auth,o=ethz,c=ch",
)
AUTH_LDAP_BIND_PASSWORD = environ.get("LDAP_PASSWORD")
AUTH_LDAP_USER_ATTR_MAP = {
    "username": "cn",
    "first_name": "givenName",
    "last_name": "sn",
    "email": "mail",
}
AUTH_LDAP_USER_SEARCH = LDAPSearch(
    environ.get("LDAP_SEARCH_BASE", "ou=users,ou=nethz,ou=id,ou=auth,o=ethz,c=ch"),
    ldap.SCOPE_SUBTREE,
    environ.get(
        "LDAP_FILTER",
        "(&(ou=@nexus.ethz.ch)(eduPersonAffiliation=staff)(cn=%(user)s))",
    ),
)
# since we've manually updated some users' names, we need to ensure they don't get overwritten
# AUTH_LDAP_ALWAYS_UPDATE_USER = False
# disable auto-creation of users since we're not taking calls for new users anymore
# AUTH_LDAP_NO_NEW_USERS = False


# Password validation
# https://docs.djangoproject.com/en/4.1/ref/settings/#auth-password-validators

AUTH_PASSWORD_VALIDATORS = [
    {
        "NAME": "django.contrib.auth.password_validation.UserAttributeSimilarityValidator",
    },
    {
        "NAME": "django.contrib.auth.password_validation.MinimumLengthValidator",
    },
    {
        "NAME": "django.contrib.auth.password_validation.CommonPasswordValidator",
    },
    {
        "NAME": "django.contrib.auth.password_validation.NumericPasswordValidator",
    },
]


# Internationalization
# https://docs.djangoproject.com/en/4.1/topics/i18n/

LANGUAGE_CODE = "en-US"
TIME_ZONE = "Europe/Zurich"
USE_I18N = True
USE_L10N = True
USE_TZ = True


# Static files (CSS, JavaScript, Images)
# https://docs.djangoproject.com/en/4.1/howto/static-files/

STATIC_URL = "/static/"
STATIC_ROOT = "/vol/web/static"

MEDIA_URL = "/media/"
MEDIA_ROOT = "/vol/web/media"

# Default primary key field type
# https://docs.djangoproject.com/en/4.1/ref/settings/#default-auto-field

DEFAULT_AUTO_FIELD = "django.db.models.BigAutoField"

DISABLE_BROWSABLE_API = False
DISABLE_AUTH = False

REST_FRAMEWORK = {
    "DEFAULT_AUTHENTICATION_CLASSES": (
        "rest_framework.authentication.SessionAuthentication",
    ),
    "DEFAULT_PERMISSION_CLASSES": [
        "rest_framework.permissions.IsAuthenticated",
    ],
    "DEFAULT_METADATA_CLASS": "meta.serializers.APIMetadata",
    "DEFAULT_RENDERER_CLASSES": (
        "rest_framework.renderers.JSONRenderer",
        "rest_framework.renderers.BrowsableAPIRenderer",
        "rest_framework_csv.renderers.CSVRenderer",
    ),
    "DEFAULT_PAGINATION_CLASS": "rest_framework.pagination.PageNumberPagination",
    "PAGE_SIZE": 10,
}

if DISABLE_BROWSABLE_API:
    REST_FRAMEWORK["DEFAULT_RENDERER_CLASSES"] = [
        "rest_framework.renderers.JSONRenderer",
        "rest_framework_csv.renderers.CSVRenderer",
    ]

if DISABLE_AUTH:
    REST_FRAMEWORK["DEFAULT_AUTHENTICATION_CLASSES"] = []
    REST_FRAMEWORK["DEFAULT_PERMISSION_CLASSES"] = []

FIXTURE_DIRS = []

###
# SECURITY
###
SITE_URL = environ.get("DJANGO_SITE_URL")
ALLOWED_HOSTS = environ.get("DJANGO_ALLOWED_HOSTS").split(",")
DISABLE_2FA = environ.get("DJANGO_DISABLE_2FA", "false").lower() == "true"

# CORS configuration
CORS_ALLOW_ALL_ORIGINS = False
CORS_ALLOWED_ORIGINS = environ.get("DJANGO_CORS_ALLOWED_ORIGINS").split(",")
CORS_ALLOW_HEADERS = default_headers + (
    "cache-control",
    "pragma",
    "expires",
    "X-CSRFTOKEN",
)
CORS_EXPOSE_HEADERS = ["Content-Type", "X-CSRFToken"]
CORS_ALLOW_CREDENTIALS = True

# CSRF configuration
CSRF_TRUSTED_ORIGINS = environ.get("DJANGO_CSRF_TRUSTED_ORIGINS").split(",")
CSRF_USE_SESSIONS = False
CSRF_COOKIE_HTTPONLY = False
CSRF_COOKIE_SAMESITE = "Strict"
SESSION_COOKIE_SAMESITE = "Strict"
SESSION_COOKIE_AGE = 1209600  # (1209600) default: 2 weeks in seconds

# PROD ONLY
# CSRF_COOKIE_SECURE = True
# SESSION_COOKIE_SECURE = True


###
# LOGGING
###
LOGGING = {
    "version": 1,
    "disable_existing_loggers": False,
    "handlers": {
        "console": {
            "class": "logging.StreamHandler",
        },
    },
    "root": {
        "level": LOG_LEVEL,
        "handlers": ["console"],
    },
    "loggers": {},
}
if LOG_SQL:
    LOGGING["loggers"]["django.db.backends"] = {
        "level": "DEBUG",
        "handlers": ["console"],
        "propagate": False,
    }
if LOG_LDAP:
    LOGGING["loggers"]["django_auth_ldap"] = {
        "level": "DEBUG",
        "handlers": ["console"],
        "propagate": False,
    }
if not LOG_DAPHNE:
    logging.getLogger("asyncio").setLevel(logging.WARNING)
    logging.getLogger("daphne").setLevel(logging.WARNING)
    logging.getLogger("daphne.http_protocol").setLevel(logging.WARNING)
    logging.getLogger("daphne.ws_protocol").setLevel(logging.WARNING)
    logging.getLogger("daphne.server").setLevel(logging.WARNING)
if not LOG_WEBSOCKET:
    logging.getLogger("websockets").setLevel(logging.WARNING)
    logging.getLogger("django.channels").setLevel(logging.WARNING)
    logging.getLogger("django.channels.server").setLevel(logging.WARNING)
    logging.getLogger("jupyter.consumers").setLevel(logging.WARNING)
    logging.getLogger("jupyter").setLevel(logging.WARNING)
if not LOG_REVPROXY:
    logging.getLogger("revproxy").setLevel(logging.WARNING)
    logging.getLogger("revproxy.cookies").setLevel(logging.WARNING)
    logging.getLogger("revproxy.response").setLevel(logging.WARNING)
if not LOG_URLLIB:
    logging.getLogger("urllib3").setLevel(logging.WARNING)
    logging.getLogger("urllib3.util").setLevel(logging.WARNING)
    logging.getLogger("urllib3.util.retry").setLevel(logging.WARNING)
    logging.getLogger("urllib3.connection").setLevel(logging.WARNING)
    logging.getLogger("urllib3.response").setLevel(logging.WARNING)
    logging.getLogger("urllib3.connectionpool").setLevel(logging.WARNING)
    logging.getLogger("urllib3.poolmanager").setLevel(logging.WARNING)

FLOAT_PRECISION = 6


"""
https://docs.djangoproject.com/en/4.1/ref/settings/
The maximum number of parameters that may be received via GET or POST before a SuspiciousOperation (TooManyFields) is raised.
You can set this to None to disable the check. 
Applications that are expected to receive an unusually large number of form fields should tune this setting.
"""
DATA_UPLOAD_MAX_NUMBER_FIELDS = 30000
USE_TZ = False

###
# Jupyter Notebook
###

NOTEBOOK_SECURE_TOKEN = environ.get("NOTEBOOK_SECURE_TOKEN", "")

NOTEBOOK_ARGUMENTS = [
    "--ip",
    "*",
    "--allow-root",
    "--no-browser",
    "--notebook-dir",
    "/notebooks",
    "--NotebookApp.token",
    NOTEBOOK_SECURE_TOKEN,
    "--NotebookApp.password",
    "",
    "--NotebookApp.base_url",
    "/notebook",
]

JUPYTER_URL = environ.get("JUPYTER_URL", "http://api:8888")

###
# Harvest
###
HARVEST_ACCESS_TOKEN = environ.get("HARVEST_ACCESS_TOKEN")
HARVEST_ACCOUNT_ID = environ.get("HARVEST_ACCOUNT_ID")
