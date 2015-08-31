import unittest
from unittest import mock
from unittest.mock import patch
import os
import json

from app import metamap, geocode, citations, umls, config

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
                'C0086418']]}])

class TestGeocode(unittest.TestCase):

    def test_remove_department(self):
        self.assertEqual(['Seoul National University College of Medicine', 
                          'Seoul National University Bundang Hospital'],
            geocode.remove_dept('Department of Orthopaedic Surgery, Seoul National University College of Medicine, Seoul National University Bundang Hospital'))

    def test_remove_email(self):
        self.assertEqual(geocode.remove_email(
            'hello address eMail Electronic hel++lo@123-bbk.ac.uk more text'), 
            'hello     more text')

    def test_format_address(self):
        self.assertEqual(geocode.format_address(
            'Department of Pharmacology, Faculty of Medical Sciences, Lagos State University College of Medicine, 1-5 Oba Akinjobi Way, G.R.A., Ikeja, Lagos State, Nigeria.'),
            ' Faculty of Medical Sciences, Lagos State University College of Medicine, 1-5 Oba Akinjobi Way, G.R.A., Ikeja, Lagos State, Nigeria.')

    def test_format_address_dummy(self):
        self.assertEqual(geocode.format_address('this, is, a, test, string'), ' is, a, test, string')

    @patch('app.geocode.format_address')
    def test_unique_addresses(self, mock_format_address): 
        mock_format_address.side_effect = ['some location',
                                            'another location',
                                            'some LOCATION',
                                            'some !!location',
                                            '  some lo cat  io n']

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

        expected = {
            'results': [
                {
                    'PMID': '00000000',
                    'place': place_result
                },
                {
                    'PMID': '00000001',
                    'place': place_result
                }
            ],
           'for_cache': [
                {
                    'PMID': '00000000',
                    'placeids': [
                        'ChIJyWEHuEmuEmsRm9hTkapTCrk'
                    ]
                },
                {
                    'PMID': '00000001',
                    'placeids': [
                        'ChIJyWEHuEmuEmsRm9hTkapTCrk'
                    ]
                }
            ]
        }

        observed = geocode.run(fake_json('eFetch_sample.json'))
        
        self.assertEqual(mock_request.call_count, 2)
        self.maxDiff = None
        self.assertEqual(observed, expected)

    @patch('app.geocode.request')
    def test_retrieve(self, mock_request):
        # patch geocode.request with example available from API docs
        mock_request.return_value = fake_json('places_details_output.json')

        expected = [{'PMID': "00000000",
            'place': {
                'name': 'Google Sydney',
                'geometry': {
                    'location': {
                        "lat": -33.8669710,
                        "lng": 151.1958750 }}}}]

        observed = geocode.retrieve(fake_json('cache_example.json'))

        self.assertEqual(mock_request.call_count, 1)
        self.maxDiff = None
        self.assertEqual(observed, expected)

class TestUmls(unittest.TestCase):
    def test_format_results(self):
        pmids_names = umls.organise(fake_json('cache_example.json'))['PMIDs']

        observed = umls.format_json(
            pmids_names,
            [{
                'CHILD_CUI': 'C0086418',
                'S_TYPE': 'Species', 
                'PARENT_CUI': 'C0086417', 
                'PARENT_STR': 'Homo Sapiens'
            }])

        expected = {
            'name': 'flare',
            'children': [{
                'name': 'Species',
                'children': [{
                    'name': 'Homo Sapiens',
                    'children': [{
                        'name': 'imaginaryconcept', 
                        'PMIDs': [
                            '00000000' ]}]}]}]}

        self.maxDiff = None
        self.assertEqual(observed, expected)

    def test_organise(self):
        observed = umls.organise(fake_json('cache_example.json'))
        expected = {'PMIDs': {'C0086418': 
                {'child_name': 'imaginaryconcept',
                'PMIDs': ['00000000']}},
            'CUIs': ['C0086418']}
        self.assertEqual(observed, expected)

    def test_group_other(self):
        expected = [{"PARENT_CUI": "C0010298", "CHILD_CUI": "C0000941", "S_TYPE": "Individual concepts", "PARENT_STR": "Credentialing"},
                    {"PARENT_CUI": "C0001779", "CHILD_CUI": "C0001675", "S_TYPE": "Age Group", "PARENT_STR": "Age"}, 
                    {"PARENT_CUI": "C0001792", "CHILD_CUI": "C0001795", "S_TYPE": "Age Group", "PARENT_STR": "Elderly (population group)"}]
        observed = umls.group_other(fake_json('umls_sql_output.json'))
        self.maxDiff = None
        self.assertEqual(observed, expected)

class TestCitations(unittest.TestCase):
    def test_format_papers_et_al(self):
        observed = citations.format_papers(fake_json('cache_example.json'))
        expected = [{'PMID': '00000000',
            'title': 'Example of a title',
            'authors': 'Heredia et al.',
            'date': '2013',
            'journal': 'Methods'
        }]

        self.maxDiff = None
        self.assertEqual(observed, expected)

if __name__ == '__main__':
    unittest.main()