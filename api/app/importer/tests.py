from django.test import TestCase
from copy import deepcopy
from .helper import sameSchema


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
