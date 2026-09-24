"""
Tests for reading Echo transfer reports (CSV) into mapping entries.
"""

from io import StringIO
from unittest import mock

from django.core.management.base import CommandError
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

    def test_a_report_without_the_optional_columns_is_read(self):
        # Some Echo exports have no "Current Fluid Volume" and no "% DMSO"
        report = (
            "Run ID,14618\n"
            "\n"
            "Source Plate Name,Source Plate Barcode,Source Plate Type,Source Well,"
            "Destination Plate Name,Destination Plate Barcode,Destination Well,"
            "Actual Volume,Transfer Status\n"
            "384LDV_DMSO,Drug08_J,384LDV_DMSO,A3,Greiner_384PS_781904,2026Wagner12,"
            "A3,10,\n"
        )

        entries = EchoMapper().parse(StringIO(report))

        self.assertEqual(["A3"], [entry["source_well"] for entry in entries])
        self.assertNotIn("DMSO", entries[0])
        self.assertNotIn("current_fluid_volume", entries[0])

    def test_a_report_without_a_needed_column_is_refused(self):
        # The report has no "Source Well", so no transfer could be read from it
        report = (
            "Source Plate Name,Source Plate Barcode,Source Plate Type,"
            "Destination Plate Name,Destination Plate Barcode,Destination Well,"
            "Actual Volume,Transfer Status\n"
            "384LDV_DMSO,Drug08_J,384LDV_DMSO,Greiner_384PS_781904,2026Wagner12,"
            "A3,10,\n"
        )

        with self.assertRaisesMessage(
            CommandError, "The report has no column Source Well."
        ):
            EchoMapper().parse(StringIO(report))

    def test_a_file_without_a_header_row_is_refused(self):
        with self.assertRaisesMessage(CommandError, "The report has no column"):
            EchoMapper().parse(StringIO("just some text\nand another line\n"))

    def test_a_report_without_a_single_transfer_is_refused(self):
        report = (
            "Source Plate Name,Source Plate Barcode,Source Plate Type,Source Well,"
            "Destination Plate Name,Destination Plate Barcode,Destination Well,"
            "Actual Volume,Transfer Status\n"
        )

        with self.assertRaisesMessage(
            CommandError, "does not contain a single transfer"
        ):
            EchoMapper().parse(StringIO(report))
