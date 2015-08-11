import os
import json
import app.metamap
from app.metamap import MetaMap

# # testing..
if __name__ == "__main__":
	txt = []
	with open(os.path.join('./tests/resources', 'metamap_output.txt'), 'r') as datafile:
		txt = datafile.read()
	txt = txt.split('\n')
	print(app.metamap.format_results(txt))
