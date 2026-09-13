"""
Session login, logout and the CSRF cookie for the web frontend.
"""

import json
from django.contrib.auth import authenticate, login, logout
from django.core.handlers.wsgi import WSGIRequest
from django.utils.decorators import method_decorator
from django.utils.translation import gettext as _
from django.views.decorators.csrf import ensure_csrf_cookie
from django.views.generic import View
from django.http import JsonResponse


class CsrfCookieView(View):
    @method_decorator(ensure_csrf_cookie)
    def get(self, request: WSGIRequest, *args, **kwargs):
        return JsonResponse({"details": _("CSRF cookie set")})


class LoginView(View):
    def post(self, request: WSGIRequest, *args, **kwargs):
        data = json.loads(request.body)
        username = data.get("username")
        password = data.get("password")

        if username is None or password is None:
            return JsonResponse(
                {"detail": _("Please provide username and password.")}, status=400
            )

        user = authenticate(username=username, password=password)

        if user is None:
            return JsonResponse({"detail": _("Invalid credentials.")}, status=400)

        login(request, user)
        return JsonResponse({"detail": _("Successfully logged in.")})


class LogoutView(View):
    def get(self, request: WSGIRequest, *args, **kwargs):
        if not request.user.is_authenticated:
            return JsonResponse({"detail": _("You're not logged in.")}, status=400)

        logout(request)
        return JsonResponse({"detail": _("Successfully logged out.")})
