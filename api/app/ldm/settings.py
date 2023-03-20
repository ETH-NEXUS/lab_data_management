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
from datetime import timedelta as td
from django_auth_ldap.config import LDAPSearch

from corsheaders.defaults import default_headers

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

# Application definition

INSTALLED_APPS = [
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
    "rest_framework_simplejwt",
    "drf_auto_endpoint",
    "core",
    "importer",
    "compoundlib",
    "platetemplate",
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

WSGI_APPLICATION = "ldm.wsgi.application"


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
        # 'rest_framework.authentication.TokenAuthentication',
        # 'rest_framework_jwt.authentication.JSONWebTokenAuthentication',
        "rest_framework_simplejwt.authentication.JWTAuthentication",
        "rest_framework.authentication.SessionAuthentication",
    ),
    "DEFAULT_PERMISSION_CLASSES": [
        # 'rest_framework.permissions.DjangoModelPermissionsOrAnonReadOnly',
        "rest_framework.permissions.IsAuthenticatedOrReadOnly",
        "rest_framework.permissions.IsAuthenticated",
    ],
    # We should not activate this if we are using drf-schema-adapter
    # 'DEFAULT_FILTER_BACKENDS': [
    #     'django_filters.rest_framework.DjangoFilterBackend'
    # ],
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

# JWT
SIMPLE_JWT = {
    "ACCESS_TOKEN_LIFETIME": td(minutes=30),
    "REFRESH_TOKEN_LIFETIME": td(days=1),
    "ROTATE_REFRESH_TOKENS": False,
    "BLACKLIST_AFTER_ROTATION": False,
    "UPDATE_LAST_LOGIN": False,
    "ALGORITHM": "HS256",
    "SIGNING_KEY": SECRET_KEY,
    "VERIFYING_KEY": None,
    "AUDIENCE": None,
    "ISSUER": None,
    "JWK_URL": None,
    "LEEWAY": 0,
    "AUTH_HEADER_TYPES": ("Bearer", "JWT"),
    "AUTH_HEADER_NAME": "HTTP_AUTHORIZATION",
    "USER_ID_FIELD": "id",
    "USER_ID_CLAIM": "user_id",
    "USER_AUTHENTICATION_RULE": "rest_framework_simplejwt.authentication.default_user_authentication_rule",
    "AUTH_TOKEN_CLASSES": ("rest_framework_simplejwt.tokens.AccessToken",),
    "TOKEN_TYPE_CLAIM": "token_type",
    "TOKEN_USER_CLASS": "rest_framework_simplejwt.models.TokenUser",
    # We overwrite the token obtain serializer to support OTP
    # 'TOKEN_OBTAIN_SERIALIZER': 'users.serializers.TokenObtainPair2FASerializer',
    "JTI_CLAIM": "jti",
    # 'SLIDING_TOKEN_REFRESH_EXP_CLAIM': 'refresh_exp',
    # 'SLIDING_TOKEN_LIFETIME': td(minutes=30),
    # 'SLIDING_TOKEN_REFRESH_LIFETIME': td(days=1),
}

###
# SECURITY
###
SITE_URL = environ.get("DJANGO_SITE_URL")
ALLOWED_HOSTS = environ.get("DJANGO_ALLOWED_HOSTS").split(",")
DISABLE_2FA = environ.get("DJANGO_DISABLE_2FA", "false").lower() == "true"

# CORS configuration
CORS_ALLOW_ALL_ORIGINS = False
CORS_ALLOWED_ORIGINS = environ.get("DJANGO_CORS_ALLOWED_ORIGINS").split(",")
CORS_ALLOW_HEADERS = default_headers + ("cache-control", "pragma", "expires")
CORS_ALLOW_CREDENTIALS = True

# CSRF configuration
CSRF_TRUSTED_ORIGINS = environ.get("DJANGO_CSRF_TRUSTED_ORIGINS").split(",")

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
    }
if LOG_LDAP:
    LOGGING["loggers"]["django_auth_ldap"] = {
        "level": "DEBUG",
        "handlers": ["console"],
    }

FLOAT_PRECISION = 6
DATA_UPLOAD_MAX_NUMBER_FIELDS = 20000

NOTEBOOK_ARGUMENTS = [
    "--ip",
    "*",
    "--allow-root",
    "--no-browser",
    "--notebook-dir",
    "/notebooks",
    "--NotebookApp.token",
    "",
    "--NotebookApp.password",
    "",
    "--NotebookApp.base_url",
    "/notebook",
]
