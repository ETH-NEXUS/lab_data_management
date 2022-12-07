from rest_framework import serializers
from .models import CompoundLibrary, Compound
from core.serializers import SimplePlateSerializer


class CompoundLibrarySerializer(serializers.ModelSerializer):
    plates = SimplePlateSerializer(many=True)

    class Meta:
        model = CompoundLibrary
        fields = '__all__'
