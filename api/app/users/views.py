from rest_framework import viewsets
from rest_framework.decorators import action
from rest_framework.response import Response
from django.contrib.auth import get_user_model
from .permissions import IsOwnUser
from .serializers import UserSerializer


class UserViewSet(viewsets.ModelViewSet):
    permission_classes = (IsOwnUser,)
    queryset = get_user_model().objects.all()
    serializer_class = UserSerializer

    @action(detail=False, methods=["get"])
    def me(self, request):
        user = get_user_model().objects.get(id=request.user.id)
        serializer = self.get_serializer(user)
        return Response(serializer.data)
