import unittest
from unittest import mock
import os
import json
from unittest.mock import patch

from app import metamap, geocode, citations

def fake_json(filename):
    with open(os.path.join('./tests/resources', filename), 'r') as datafile:
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
        flag = metamap.write_file('metamap_input.txt', fake_json('eFetch_sample.json'))
        self.assertTrue(flag)
        
        test_txt = load_read_close('./app/static', 'metamap_input.txt')
        assert_txt = load_read_close('./tests/resources', 'metamap_input.txt')

        self.maxDiff = None # able to see all differences during testing
        self.assertEqual(test_txt, assert_txt)
        os.remove(os.path.join('./app/static', 'metamap_input.txt'))

    def test_format_results(self):
        self.assertEqual(metamap.format_results(
                ["23036330|Humans|C0086418|1000|CT|Breast Cancer;CT Treecode Lookup: C17.800.090 (Breast Cancer);CT Text Lookup: human||MM;RC"]),
            [{'PMID': '23036330',
              'concepts': [['Humans',
                           'C0086418',
                           '1000',
                           'CT',
                           'Breast Cancer;CT Treecode Lookup: C17.800.090 (Breast Cancer);CT Text Lookup: human',
                           '',
                           'MM;RC']]}])

class TestGeocode(unittest.TestCase):

    def test_remove_email(self):
        self.assertEqual(geocode.remove_email(
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
        self.assertEqual(geocode.unique_addresses([{'AffiliationInfo' : [{'Affiliation' : 'some location;another location'}]},
                                                   {'AffiliationInfo' : [{'Affiliation' : 'some LOCATION'}]},
                                                   {'AffiliationInfo' : [{'Affiliation' : 'some !!location'}]},
                                                   {'AffiliationInfo' : [{'Affiliation' : '  some lo cat  io n'}]}]), 
            ['some location', 'another location'])

    @patch('app.geocode.request')
    def test_run(self, mock_request):
        # patch geocode.request with example available from API docs
        mock_request.return_value = fake_json('places_search_output.json')
        # get_location_output.json contains the first element of the results list
        place_result = fake_json('get_location_output.json')

        expected = [
            {'PMID': '00000000',
             'places': [place_result]},
            {'PMID': '00000001',
             'places': [place_result]}
        ]

        observed = geocode.run(fake_json('eFetch_sample.json'))
        self.maxDiff = None
        self.assertEqual(observed, expected)

    @patch('app.geocode.request')
    def test_retrieve(self, mock_request):
        # patch geocode.request with example available from API docs
        mock_request.return_value = fake_json('places_details_output.json')

        expected = [{
            'PMID': "00000000",
            'results': [{
                'name': 'Google Sydney',
                'geometry': {
                    'location': {
                        "lat": -33.8669710,
                        "lng": 151.1958750
                    }
                }
            }]
        }]

        observed = geocode.retrieve(fake_json('cache_example.json'))
        self.maxDiff = None
        self.assertEqual(observed, expected)

# class TestCitations(unittest.TestCase):

if __name__ == '__main__':
    unittest.main()