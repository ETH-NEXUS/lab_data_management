"""
Tests for reading Echo transfer reports (CSV) into mapping entries.
"""

from io import StringIO
from unittest import mock

from django.test import SimpleTestCase

from importer.mappers import EchoMapper

HEADER = (
    "Source Plate Name,Source Plate Barcode,Source Plate Type,Source Well,"
    "Source Concentration,Source Concentration Units,Destination Plate Name,"
    "Destination Plate Barcode,Destination Well,Destination Concentration,"
    "Destination Concentration Units,Compound Name,Transfer Volume,Actual Volume,"
    "Transfer Status,Current Fluid Height,Current Fluid Volume,% DMSO"
)

# A shortened report with both sections the Echo writes: failed transfers
# under [EXCEPTIONS], then all transfers under [DETAILS].
REPORT_WITH_EXCEPTIONS = f"""Run ID,446,,,,,,,,,,,,,,,,
User Name,cellario,,,,,,,,,,,,,,,,
,,,,,,,,,,,,,,,,,
[EXCEPTIONS],,,,,,,,,,,,,,,,,
{HEADER}
384LDV_DMSO,Drug08_J,384LDV_DMSO,B4,0,N/A,Greiner_384PS_781904,2026Wagner12,B4,0,N/A,N/A,10,0,MM0202007: Problem calc. well fluid volume fc: 0 ft: 0,0,0,0
[DETAILS],,,,,,,,,,,,,,,,,
{HEADER}
384LDV_DMSO,Drug08_J,384LDV_DMSO,A3,0,N/A,Greiner_384PS_781904,2026Wagner12,A3,0,N/A,N/A,10,10,,1.996,9.957,99.497
384LDV_DMSO,Drug08_J,384LDV_DMSO,A4,0,N/A,Greiner_384PS_781904,2026Wagner12,A4,0,N/A,N/A,10,10,,1.931,9.647,99.597
"""


# A real transfer (A5) where the Echo left "Actual Volume" empty.
REPORT_WITH_EMPTY_VOLUME = f"""Run ID,446,,,,,,,,,,,,,,,,
[DETAILS],,,,,,,,,,,,,,,,,
{HEADER}
384LDV_DMSO,Drug08_J,384LDV_DMSO,A3,0,N/A,Greiner_384PS_781904,2026Wagner12,A3,0,N/A,N/A,10,10,,1.996,9.957,99.497
384LDV_DMSO,Drug08_J,384LDV_DMSO,A5,0,N/A,Greiner_384PS_781904,2026Wagner12,A5,0,N/A,N/A,10,,,1.984,9.886,99.086
"""


# A shortened Echo XML report: two transfers in <printmap>, and a plate
# measurement under <platemap> that also uses <w> elements.
XML_REPORT = """<?xml version="1.0"?>
<transfer style="startstop" date="2024-06-21 09:28:06.000"><plateInfo><plate type="source" name="384LDV_DMSO" barcode="24MP_Test1_1"/><plate type="destination" name="Corning_384_4516" barcode="240618MP-1_2"/></plateInfo><printmap total="2"><w n="A1" r="0" c="0" dn="A1" dr="0" dc="0" cvl="5.543" vl="5.583" fld="DMSO" fc="99.701" reason=""/><w n="B2" r="1" c="1" dn="C3" dr="2" dc="2" cvl="5.262" vl="5.302" fld="DMSO" fc="98.443" reason=""/></printmap><platemap><measbottof><w n="P24" r="15" c="23" mbtof="32.414"/></measbottof></platemap></transfer>
"""


class EchoParseTest(SimpleTestCase):
    def test_section_lines_are_not_read_as_transfers(self):
        # The "[DETAILS],,,," line used to become an entry with empty wells,
        # which crashed the mapping with a PositionMappingError.
        entries = EchoMapper().parse(StringIO(REPORT_WITH_EXCEPTIONS))

        source_wells = [entry["source_well"] for entry in entries]
        self.assertEqual(["B4", "A3", "A4"], source_wells)

    @mock.patch("importer.mappers.echo.message")
    def test_section_lines_are_skipped_without_a_warning(self, message):
        EchoMapper().parse(StringIO(REPORT_WITH_EXCEPTIONS), room_name="room_1")

        message.assert_not_called()

    @mock.patch("importer.mappers.echo.message")
    def test_a_transfer_with_an_empty_field_is_skipped_with_one_warning(self, message):
        entries = EchoMapper().parse(
            StringIO(REPORT_WITH_EMPTY_VOLUME), room_name="room_1"
        )

        self.assertEqual(["A3"], [entry["source_well"] for entry in entries])
        message.assert_called_once_with(
            "Skipped 1 Echo transfers with empty required fields: "
            "Drug08_J A5 -> 2026Wagner12 A5 (Actual Volume empty)",
            "warning",
            "room_1",
        )

    def test_an_xml_report_is_read_from_the_printmap(self):
        # Only the <w> elements of <printmap> are transfers; the <w> elements
        # under <platemap> are plate measurements and must be ignored.
        entries = EchoMapper().parse(StringIO(XML_REPORT), xml_file=True)

        plates = {
            "source_plate_name": "384LDV_DMSO",
            "source_plate_barcode": "24MP_Test1_1",
            "destination_plate_name": "Corning_384_4516",
            "destination_plate_barcode": "240618MP-1_2",
            "transfer_status": "",
        }
        self.assertEqual(
            [
                {
                    **plates,
                    "source_well": "A1",
                    "destination_well": "A1",
                    "actual_volume": "5.583",
                    "current_fluid_volume": "5.543",
                    "DMSO": "99.701",
                },
                {
                    **plates,
                    "source_well": "B2",
                    "destination_well": "C3",
                    "actual_volume": "5.302",
                    "current_fluid_volume": "5.262",
                    "DMSO": "98.443",
                },
            ],
            entries,
        )
