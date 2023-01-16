from rest_framework import serializers
from .models import CompoundLibrary, Compound
from core.serializers import SimplePlateSerializer, WellPlateSerializer


class CompoundLibrarySerializer(serializers.ModelSerializer):
    plates = SimplePlateSerializer(many=True)

    class Meta:
        model = CompoundLibrary
        fields = '__all__'


class SimpleCompoundLibrarySerializer(serializers.ModelSerializer):
    class Meta:
        model = CompoundLibrary
        fields = ('id', 'name')


class CompoundSerializer(serializers.ModelSerializer):
    wells = serializers.SerializerMethodField()

    def get_wells(self, compound: Compound):
        wells = []
        for well_compound in compound.well_compounds.all():
            wells.append(WellPlateSerializer(well_compound.well).data)
        return wells

    class Meta:
        model = Compound
        fields = '__all__'
