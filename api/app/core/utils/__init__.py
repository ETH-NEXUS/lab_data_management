"""
Helpers of the core app that are not models, views or serializers.

- `config.py`  reading and checking the LDM configuration file (ldm.yaml)
- `plates/`    helpers for plates: well positions, plate mappings, copying, archiving
- `wells/`     helpers for the content of wells: volume units, threshold checks

Import from the module itself, e.g. `from core.utils.plates.mapping import MappingList`.
Nothing is re-exported here on purpose, so the import shows where the code lives.
"""
