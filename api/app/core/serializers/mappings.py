"""
Plate mappings: applying a mapping (csv file or plate copy) when it is saved.
"""

from rest_framework import serializers

from ..models import MappingError, PlateDetail, PlateMapping, WellDetail
from ..utils.plates.mapping import MappingList
from .base import UndefinedAffineModelSerializer


class PlateMappingSerializer(UndefinedAffineModelSerializer):
    def save(self, **kwargs):
        try:
            # Since it's a new object creation, 'id' won't be in validated_data
            mapping_file = self.validated_data.get("mapping_file")
            from_column = self.validated_data.get("from_column")
            to_column = self.validated_data.get("to_column")
            amount_column = self.validated_data.get("amount_column")
            delimiter = self.validated_data.get("delimiter")
            quotechar = self.validated_data.get("quotechar")
            sourcePlate = self.validated_data[
                "source_plate"
            ]  # it needs the id and not a plate object
            targetPlate = self.validated_data[
                "target_plate"
            ]  # it needs the id and not a plate object

            if mapping_file and from_column and to_column:
                sourcePlate.map(
                    MappingList.from_csv(
                        mapping_file["name"],
                        from_column,
                        to_column,
                        amount_column,
                        delimiter,
                        quotechar,
                    ),
                    targetPlate,
                )
            else:
                sourcePlate.copy(targetPlate, self.validated_data["amount"])
            # TODO: Additional logic for update or delete, if necessary
        except MappingError as ex:
            raise serializers.ValidationError({"detail": ex})

        PlateDetail.refresh(concurrently=True)
        WellDetail.refresh(concurrently=True)
        return super().save(**kwargs)

    class Meta:
        model = PlateMapping
        fields = "__all__"
