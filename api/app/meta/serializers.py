from rest_framework.metadata import SimpleMetadata
from rest_framework.schemas.openapi import AutoSchema


def get_label_translation(label):
    label_list = label.split(">>>")
    if len(label_list) > 1 and label_list[1] != "None":
        return {"en": label_list[0], "de": label_list[1]}
    else:
        return {"en": label_list[0], "de": label_list[0]}


class APIMetadata(SimpleMetadata):
    # TODO Dynamic depth?
    def get_field_info(self, field):
        field_info = super().get_field_info(field)

        # Add extra validators using the OpenAPI schema generator.
        validators = {}
        AutoSchema().map_field_validators(field, validators)
        extra_validators = ["format", "pattern"]
        for validator in extra_validators:
            if validators.get(validator, None):
                field_info[validator] = validators[validator]

        # Add additional data from serializer.
        field_info["initial"] = field.initial
        field_info["field"] = field.field_name
        field_info["write_only"] = field.write_only
        if hasattr(field, "choices") and field.choices:
            field_info["type"] = "single_select"
            field_info["choices"] = []

            for value, label in field.choices.items():
                label_translation = get_label_translation(label)
                choice = {
                    "label_de": label_translation["de"],
                    "label": label_translation["en"],
                    "value": value,
                }
                # You can flag a foreign key field value (__str__) with '__DEFAULT__'
                if label.startswith("__DEFAULT__"):
                    choice["label"] = choice["label"].replace("__DEFAULT__", "")
                    choice["label_de"] = choice["label_de"].replace("__DEFAULT__", "")

                    choice["default"] = True
                field_info["choices"].append(choice)

        return field_info
