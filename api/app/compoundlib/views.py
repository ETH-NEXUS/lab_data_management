from rest_framework import viewsets
from .serializers import CompoundLibrarySerializer
from .models import CompoundLibrary


# class CompoundLibraryViewSet(viewsets.ModelViewSet):
#     serializer_class = CompoundLibrarySerializer
#     filterset_fields = ('name',)
