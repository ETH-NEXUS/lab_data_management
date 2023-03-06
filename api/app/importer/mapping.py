import json
import yaml
from .helper import sameSchema
from os.path import isfile


class MappingFileSchemaError:
    def __init__(
        self, mapping={}, message="Error in mapping file schema. Should be {}"
    ):
        self.message = message.format(mapping)
        super().__init__(self.message)


class SdfMapping:
    DEFAULT_MAPPING = {
        "compound": {
            "identifier": "IDNUMBER",
            "name": "NAME",
            "structure": "STRUCTURE",
        },
        "plate": {
            "barcode": "PLATE_NUMBER1",
            "position": "POS_IN_PLATE",
            "amount": "PLATE_AMOUNT1",
        },
    }

    def __init__(self, mappingFile: str = None):
        self.mapping = self.DEFAULT_MAPPING
        if mappingFile:
            self.load(mappingFile)

    def __str__(self):
        return json.dumps(self.mapping)

    def __repr__(self):
        return self.__str__()

    def load(self, mappingFile: str):
        if not isfile(mappingFile):
            raise FileNotFoundError(f"Cannot find mapping file {mappingFile}.")

        with open(mappingFile, "r") as mf:
            mapping = yaml.load(mf, yaml.SafeLoader)
        if not sameSchema(self.mapping, mapping):
            raise MappingFileSchemaError(self.DEFAULT_MAPPING)

        self.mapping.update(mapping)

    @property
    def name(self) -> str:
        return self.mapping["compound"]["name"]

    @property
    def identifier(self) -> str:
        return self.mapping["compound"]["identifier"]

    @property
    def structure(self) -> str:
        return self.mapping["compound"]["structure"]

    @property
    def barcodes(self) -> tuple[str]:
        if isinstance(self.mapping["plate"]["barcode"], str):
            return (self.mapping["plate"]["barcode"],)
        return tuple(self.mapping["plate"]["barcode"])

    @property
    def position(self) -> str:
        return self.mapping["plate"]["position"]

    @property
    def amounts(self) -> tuple[float]:
        if isinstance(self.mapping["plate"]["amount"], str):
            return (self.mapping["plate"]["amount"],)
        return tuple(self.mapping["plate"]["amount"])
