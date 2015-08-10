import unittest
import os
import unittest.mock
import json
from unittest.mock import patch

from app.metamap import MetaMap
from app import citations

def fake_input():
    obj = []
    with open(os.path.join('./tests/resources', 'eFetch_sample.json'), 'r') as datafile:
        obj = json.load(datafile)
        datafile.close()
    return obj

def fake_output():
    file_txt = open(os.path.join('./tests/resources', 'metamap_output.txt'), 'r')
    output = file_txt.read()
    file_txt.close()
    return output

class TestMetaMap(unittest.TestCase):
    def setUp(self):
        self.patcher = patch('app.metamap.write_file', fake_output())
        self.patcher.start()
        self.mm = MetaMap()

    def tearDown(self):
        self.patcher.stop()

    def test_run(self):
        response = self.mm.run(fake_input())
        self.assertIn('23036330', response)

class TestCitations(unittest.TestCase):

    def test_remove_email(self):
        self.assertEqual(citations._remove_email(
            'hello address eMail Electronic hel++lo@123-bbk.ac.uk more text'), 
            'hello     more text')

    # tests work but I keep changing how many lines to keep

    # def test_format_address(self):
    #   self.assertEqual(citations._format_address(
    #       'Department of Pharmacology, Faculty of Medical Sciences, Lagos State University College of Medicine, 1-5 Oba Akinjobi Way, G.R.A., Ikeja, Lagos State, Nigeria.'),
    #       ' G.R.A., Ikeja, Lagos State, Nigeria.')

    # def test_format_address_dummy(self):
    #   self.assertEqual(citations._format_address('this, is, a, test, string'), ' is, a, test, string')

    def test_unique_addresses(self): 
        self.assertEqual(citations.unique_addresses([{'AffiliationInfo' : [{'Affiliation' : 'some location;another location'}]},
                                                     {'AffiliationInfo' : [{'Affiliation' : 'some LOCATION'}]},
                                                     {'AffiliationInfo' : [{'Affiliation' : 'some !!location'}]},
                                                     {'AffiliationInfo' : [{'Affiliation' : '  some lo cat  io n'}]}]), 
            ['some location', 'another location'])

    # @mock.patch('citations.search')
    # def test_start_search(self, mock_search):
    #     mock_search.start_search('any query').assert

if __name__ == '__main__':
    unittest.main()