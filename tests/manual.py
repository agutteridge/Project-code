import os

from app import citations

def init():
    with open(os.path.join('./tests/resources', 'eFetch_random.json'), 'a') as datafile:
        count = 0

        while count < 100:
            id_list = citations._search('biology')['IdList']
            fetch_results = fetch_details(id_list)
            count += len(fetch_results)
            datafile.write(results)