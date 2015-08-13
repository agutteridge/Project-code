import os
import json
import app.metamap
import app.umls
from app.metamap import MetaMap

def load_read_close(path, filename):
    with open(os.path.join(path, filename), 'r') as datafile:
        txt = datafile.read()
        datafile.close()
        return txt

# # testing..
if __name__ == "__main__":
    # with open(os.path.join('./tests/resources', 'eFetch_sample.json'), 'r') as datafile:
    #     mm = MetaMap()
    #     mm.run(json.load(datafile))
    #     datafile.close()

    test_txt = load_read_close('./tests/resources', 'metamap_output.txt').split('\n')
    paper_terms = app.metamap.format_results(test_txt)
    print(app.umls.run_with_data(paper_terms))