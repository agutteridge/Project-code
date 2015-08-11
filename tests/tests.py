import unittest
import os
import json

from app import metamap
from app import citations

def fake_json():
    with open(os.path.join('./tests/resources', 'eFetch_sample.json'), 'r') as datafile:
        obj = json.load(datafile)
        datafile.close()
        return obj

def load_read_close(path, filename):
    with open(os.path.join(path, filename), 'r') as datafile:
        txt = datafile.read()
        datafile.close()
        return txt

class TestMetaMap(unittest.TestCase):

    def test_write_file(self):
        flag = metamap.write_file('metamap_input.txt', fake_json())
        self.assertTrue(flag)
        
        test_txt = load_read_close('./app/static', 'metamap_input.txt')
        assert_txt = load_read_close('./tests/resources', 'metamap_input.txt')

        self.maxDiff = None # able to see all differences during testing
        self.assertEqual(test_txt, assert_txt)
        os.remove(os.path.join('./app/static', 'metamap_input.txt'))

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

if __name__ == '__main__':
    unittest.main()