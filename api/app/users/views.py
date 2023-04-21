from django.http import JsonResponse
from rest_framework import viewsets
from rest_framework.decorators import action
from rest_framework.response import Response
from django.contrib.auth import get_user_model
from .permissions import IsOwnUser
from .serializers import UserSerializer
from rest_framework.exceptions import AuthenticationFailed
from rest_framework_simplejwt.authentication import JWTAuthentication


class UserViewSet(viewsets.ModelViewSet):
    permission_classes = (IsOwnUser,)
    queryset = get_user_model().objects.all()
    serializer_class = UserSerializer

    @action(detail=False, methods=["get"])
    def me(self, request):
        user = get_user_model().objects.get(id=request.user.id)
        serializer = self.get_serializer(user)
        return Response(serializer.data)


from django.contrib.auth.decorators import user_passes_test


def check_auth(request):
    if not request.user.is_authenticated:
        response = JsonResponse({"authenticated": False})
        response.status_code = 401
        return response
    else:
        return JsonResponse({"authenticated": True})


# def authenticate_token(request):
#     auth = JWTAuthentication()
#     try:
#         user_and_token = auth.authenticate(request)
#         if user_and_token is not None:
#             user, token = user_and_token
#             return user
#     except AuthenticationFailed:
#         pass
#
#     return None
#
#
# def check_auth(request):
#     user = request.user
#     if not user.is_authenticated:
#         user = authenticate_token(request)
#
#     if user and user.is_authenticated:
#         return JsonResponse({"authenticated": True})
#     else:
#         response = JsonResponse({"authenticated": False})
#         response.status_code = 401
#         return response
