import os
import json
import app.metamap
from app.metamap import MetaMap

# # testing..
if __name__ == "__main__":
    with open(os.path.join('./tests/resources', 'eFetch_sample.json'), 
            'r') as datafile:
        mm = MetaMap()
        mm.run(json.load(datafile))
