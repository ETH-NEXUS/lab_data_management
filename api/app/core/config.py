import yaml
from easydict import EasyDict
from django.utils.functional import classproperty
from schema import Schema, SchemaError
from helpers.logger import logger


class Config:
    """
    The configuration of LDM.

    USAGE:

    from core.config import Config

    Config.current. ...
    """

    config_schema = Schema(
        {
            "importer": {
                "echo": {
                    "default": {
                        "file_blob": str,
                        "xml_blob": str,
                        "columns": {
                            "source_plate_barcode": str,
                            "source_plate_type": str,
                            "source_plate_name": str,
                            "source_well": str,
                            "destination_plate_name": str,
                            "destination_plate_barcode": str,
                            "destination_plate_type": str,
                            "destination_well": str,
                            "actual_volume": str,
                            "transfer_status": str,
                            "current_fluid_volume": str,
                            "DMSO": str,
                        },
                    }
                },
                "m1000": {"default": {"file_blob": str}},
                "microscope": {"default": {"file_blob": str, "txt_blob": str}},
                "dat": {"default": {"file_blob": str}},
            }
        }
    )

    @classproperty
    def current(cls):
        """
        This makes sure the config is read every time it's accessed,
        with the advantage to always have the latest configuration and
        the disadvantage that it takes a bit longer to return the config.
        """
        with open("./ldm.yaml", "r") as cf:
            cfg = yaml.safe_load(cf)
            try:
                cls.config_schema.validate(cfg)
                return EasyDict(cfg)
            except SchemaError as se:
                logger.critical(se)
                raise se
