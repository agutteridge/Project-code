import unittest, mock, os
import citations, metamap 

def fake_metamap_run():
    with open(os.path.join('/resources', 'loremipsum_for_metamap.txt', 
            'rb')) as datafile:
        return json.load(datafile)

class TestMetaMap(unittest.TestCase):
    def setUp(self):
        self.patcher = patch('metamap.run', fake_urlopen)
        self.patcher.start()
        self.

    def tearDown(self):
        self.patcher.stop()

    def test_run(self):
        response = metamap.run(user)
        self.assertIn('name', response)
        self.assertEqual(response['name'], 'Test User')


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
            {'some location', 'another location'})

    # @mock.patch('citations.search')
    # def test_start_search(self, mock_search):
    #     mock_search.start_search('any query').assert

if __name__ == '__main__':
    unittest.main()