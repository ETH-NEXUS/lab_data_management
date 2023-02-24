from django.test import TestCase
from copy import deepcopy
from .helper import sameSchema
from importer.mapping import EchoMapping


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



class ConvertPositionToIndexTestCase(TestCase):
    def test_single_letter(self):
        '''
        Test a position with a single-letter row label
        '''

        position = "B3"
        number_of_columns = 4
        expected_index = 7
        self.assertEqual(
            EchoMapping.convert_position_to_index(position, number_of_columns),
            expected_index)

    def test_first_column(self):
        '''
        Test a position in the first column
        '''

        position = "A1"
        number_of_columns = 4
        expected_index = 0
        self.assertEqual(
            EchoMapping.convert_position_to_index(position, number_of_columns),
            expected_index)

    def test_last_column(self):
        '''
        Test a position in the last column
        '''

        position = "D3"
        number_of_columns = 4
        expected_index = 11
        self.assertEqual(
            EchoMapping.convert_position_to_index(position, number_of_columns),
            expected_index)


    def test_multi_letter(self):
        '''
        Test a position with a multi-letter row label
        '''

        position = "AA11"
        number_of_columns = 4
        expected_index = 105
        self.assertEqual(
            EchoMapping.convert_position_to_index(position, number_of_columns),
            expected_index)

