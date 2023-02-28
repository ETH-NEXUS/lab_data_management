from django.test import TestCase
from copy import deepcopy
from .helper import sameSchema
from .mappers import BaseMapper
from os import makedirs
from os.path import join
from core.models import Plate


class HelperTests(TestCase):
  def test_sameSchemaTrue(self):
    a = {
        'a': {
            'b': {
                'c': 'xxx'
            },
            'x': {
                'y': 111
            }
        }
    }
    b = deepcopy(a)
    self.assertTrue(sameSchema(a, b))

  def test_sameSchemaFalse(self):
    a = {
        'a': {
            'b': {
                'c': 'xxx'
            },
            'x': {
                'y': 111
            }
        }
    }
    b = deepcopy(a)
    del b['a']['x']['y']
    b['a']['x']['z'] = 'changed'
    self.assertFalse(sameSchema(a, b))



class ConvertPositionToIndexTests(TestCase):
    def test_single_letter(self):
        '''
        Test a position with a single-letter row label
        '''

        position = "B3"
        number_of_columns = 4
        expected_index = 7
        self.assertEqual(
            Plate.convert_position_to_index(position, number_of_columns),
            expected_index)

    def test_first_column(self):
        '''
        Test a position in the first column
        '''

        position = "A1"
        number_of_columns = 4
        expected_index = 0
        self.assertEqual(
            Plate.convert_position_to_index(position, number_of_columns),
            expected_index)

    def test_last_column(self):
        '''
        Test a position in the last column
        '''

        position = "D3"
        number_of_columns = 4
        expected_index = 11
        self.assertEqual(
            Plate.convert_position_to_index(position, number_of_columns),
            expected_index)


    def test_multi_letter(self):
        '''
        Test a position with a multi-letter row label
        '''

        position = "AA11"
        number_of_columns = 4
        expected_index = 105
        self.assertEqual(
            Plate.convert_position_to_index(position, number_of_columns),
            expected_index)

class MapperTests(TestCase):
    TEST_DATA_FOLDER = './temp'
    BARCODES = [
        'P1',
        'P2',
        'P3'
    ]
    ECHO_FILES = [

    ]
    def create_library_plates(self):
        pass


    def create_echo_test_data(self):
        ECHO_DIR = join(self.TEST_DATA_FOLDER, 'echo')

        for barcode in self.BARCODES:
            _dir = join(ECHO_DIR, barcode)
            makedirs(_dir, exist_ok=True)
            with open(join(_dir, 'ID-123-transfer-Echo_01_123.csv'), 'w') as\
                    ef:
                ef.write(
                    f"""
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
            
                


    def startUp(self):
        self.create_echo_test_data()


    def test_get_files(self):
        mapper = BaseMapper()
        self.assertEqual('',mapper.get_files('/data/echo/'))
