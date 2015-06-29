import unittest
import printing_affiliations

class TestPrinting_Affiliations(unittest.TestCase):

    def test_remove_email(self):
        self.assertEqual(printing_affiliations.remove_email('address hel++lo@123-bbk.ac.uk more text'), 'address more text')

if __name__ == '__main__':
    unittest.main()