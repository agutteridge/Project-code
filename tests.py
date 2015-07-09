import unittest
import printing_affiliations

class TestPrinting_Affiliations(unittest.TestCase):

    def test_remove_email(self):
        self.assertEqual(printing_affiliations._remove_email(
        	'hello address email electronic hel++lo@123-bbk.ac.uk more text'), 
        	'hello     more text')

    def test_format_address(self):
    	self.assertEqual(printing_affiliations._format_address(
    		'Department of Pharmacology, Faculty of Medical Sciences, Lagos State University College of Medicine, 1-5 Oba Akinjobi Way, G.R.A., Ikeja, Lagos State, Nigeria.'),
    		' G.R.A., Ikeja, Lagos State, Nigeria.')

    def test_format_address_dummy(self):
    	self.assertEqual(printing_affiliations._format_address('this, is, a, test, string'), ' is, a, test, string')


if __name__ == '__main__':
    unittest.main()