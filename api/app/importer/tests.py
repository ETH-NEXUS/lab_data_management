from copy import deepcopy
from os import makedirs
from os.path import join
from pathlib import Path

from django.test import TestCase

from core.helper import posToAlphaChar
from core.models import PlateDimension
from .helper import sameSchema
from .mappers import BaseMapper, EchoMapper


class HelperTests(TestCase):

    def test_sameSchemaTrue(self):
        a = {'a': {'b': {'c': 'xxx'}, 'x': {'y': 111}}}
        b = deepcopy(a)
        self.assertTrue(sameSchema(a, b))

    def test_sameSchemaFalse(self):
        a = {'a': {'b': {'c': 'xxx'}, 'x': {'y': 111}}}
        b = deepcopy(a)
        del b['a']['x']['y']
        b['a']['x']['z'] = 'changed'
        self.assertFalse(sameSchema(a, b))


class ConvertPositionToIndexTests(TestCase):
    fixtures = ['plate_dimensions']

    def setUp(self):
        self.dimension_96 = PlateDimension.objects.get(name='dim_96_8x12')
        self.dimension_384 = PlateDimension.objects.get(name='dim_384_16x24')
        self.dimension_1536 = PlateDimension.objects.get(name='dim_1536_32x48')

    def test_single_letter(self):
        '''
        Test a position with a single-letter row label
        '''

        position = "B3"
        expected_index = 26
        self.assertEqual(self.dimension_384.position(position), expected_index)

    def test_first_column(self):
        '''
        Test a position in the first column
        '''

        position = "A1"
        expected_index = 0
        self.assertEqual(self.dimension_96.position(position), expected_index)

    def test_random_column(self):
        '''
        Test a position in the last column
        '''

        position = "D3"
        expected_index = 74
        self.assertEqual(self.dimension_384.position(position), expected_index)

    def test_last_column(self):
        '''
        Test a position in the last column
        '''

        position = "AA11"
        expected_index = 1258
        self.assertEqual(self.dimension_1536.position(position),
                         expected_index)

    def test_index_to_letter(self):
        '''
        Test a position in the last column
        '''

        index = 27
        expected_letter = 'AA'
        self.assertEqual(posToAlphaChar(index), expected_letter)
        self.assertEqual(posToAlphaChar(28), "AB")
        self.assertEqual(posToAlphaChar(52), "AZ")
        self.assertEqual(posToAlphaChar(703), "AAA")
        self.assertEqual(posToAlphaChar(73116), "DDDD")


class MapperTests(TestCase):
    fixtures = ['plate_dimensions', 'well_types', 'test/compound_library']

    TEST_DATA_FOLDER = './temp'
    ECHO_DIR = join(TEST_DATA_FOLDER, 'echo')
    BARCODES = ['P1', 'P2', 'P3']

    def create_echo_test_data(self):
        for barcode in self.BARCODES:
            _dir = join(self.ECHO_DIR, barcode)
            Path(join(_dir, 'something.csv')).touch()
            Path(join(_dir, 'dfg-transfer-file.tsv')).touch()
            makedirs(_dir, exist_ok=True)
            with open(join(_dir, 'ID-123-transfer-Echo_01_123.csv'),
                      'w') as ef:
                ef.write(f"""
                    Run ID,14618
                    Run Date/Time,02/09/2021 10:30:32
                    Application Name,
                    Application Version,1.0.0.0
                    Protocol Name,
                    User Name,cellario
                    
                    
                    [DETAILS]
                    Source Plate Name,Source Plate Barcode,Source Plate Type,Source Well,Source Concentration,Source Concentration Units,Destination Plate Name,Destination Plate Barcode,Destination Well,Destination Concentration,Destination Concentration Units,Compound Name,Transfer Volume,Actual Volume,Transfer Status,Current Fluid Height,Current Fluid Volume,% DMSO
                    384LDV_DMSO,LLD_4541_C,384LDV_DMSO,A10,0,N/A,Corning_384_3577,{barcode},A10,0,N/A,N/A,30,30,,2.003,9.959,99.203
                    384LDV_DMSO,LLD_4541_C,384LDV_DMSO,A11,0,N/A,Corning_384_3577,{barcode},A11,0,N/A,N/A,30,30,,1.968,9.786,99.008
                    384LDV_DMSO,LLD_4541_C,384LDV_DMSO,A12,0,N/A,Corning_384_3577,{barcode},A12,0,N/A,N/A,30,30,,2.001,9.959,99.348
                    384LDV_DMSO,LLD_4541_C,384LDV_DMSO,A13,0,N/A,Corning_384_3577,{barcode},A13,0,N/A,N/A,30,30,,2.013,9.983,98.519
                    384LDV_DMSO,LLD_4541_C,384LDV_DMSO,A14,0,N/A,Corning_384_3577,{barcode},A14,0,N/A,N/A,30,30,,2.001,9.973,99.701
                    384LDV_DMSO,LLD_4541_C,384LDV_DMSO,A15,0,N/A,Corning_384_3577,{barcode},A15,0,N/A,N/A,30,30,,1.964,9.793,99.643
                    """)

    def setUp(self):
        self.create_echo_test_data()

    def tearDown(self):
        # TODO: Delete files
        pass

    def test_get_files(self):
        mapper = BaseMapper()
        self.assertListEqual(
            [join(self.ECHO_DIR, barcode, 'ID-123-transfer-Echo_01_123.csv')
                for barcode in self.BARCODES],
            mapper.get_files(join(self.ECHO_DIR, '**', '*-transfer-*.csv')))

    def test_parse_echo_file(self):
        mapper = EchoMapper()
        with open(join(self.ECHO_DIR, self.BARCODES[0],
            'ID-123-transfer-Echo_01_123.csv'), 'r') as file:
            data = mapper.parse('test.csv', file)
            for item in data:
                print('\n', item)
        self.assertEqual(len(data), 6)
        self.assertDictEqual(data[0], {
            'source_plate_barcode': 'LLD_4541_C',
            'destination_plate_barcode': self.BARCODES[0]
        })
