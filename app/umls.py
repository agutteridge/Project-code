import pymysql.cursors
from collections import defaultdict

from app import config

# Visualisation of concepts is crowded if many have their own hierarchy.
# This function groups all concepts with unique semantic types, by changing
# the S_TYPE to 'Individual concepts'. 
def group_other(results):
    dd = defaultdict(int)

    for r in results:
        dd[r['S_TYPE']] += 1

    for r in results:
        s_type = r['S_TYPE']

        if dd[s_type] == 1:
            r['S_TYPE'] = 'Individual concepts'

    return results

# Formatting Python dict for D3
def format_json(pmids_names, results):
    dict0 = dict()
    # level 0: all data
    dict0['name'] = 'flare'
    dict0['children'] = []

    for r in results:
        # level 1: Semantic Types
        st_key = r['S_TYPE']

        found = False
        for c in dict0['children']:
            if c['name'] == st_key:
                # No need to add Semantic type
                found = True
                dict1 = c
                
        if not found:
            dict1 = dict()
            dict1['children'] = []
            dict0['children'].append(dict1)
            dict1['name'] = st_key

        # level 2: Parent names
        par_key = r['PARENT_STR']

        found = False
        for c in dict1['children']:
            if c['name'] == par_key:
                # No need to add Parent name                
                found = True
                dict2 = c

        if not found:
            dict2 = dict()
            dict2['children'] = []
            dict1['children'].append(dict2)
            dict2['name'] = par_key

        # level 3: Child names
        c_cui = r['CHILD_CUI']
        c_key = pmids_names[c_cui]['child_name']
        pmids = pmids_names[c_cui]['PMIDs']

        dict3 = dict()
        dict3['PMIDs'] = pmids
        dict2['children'].append(dict3)
        dict3['name'] = c_key
        
    return dict0

# Arg is a list of CUIs to retrieve information for
def execute_sql(cui_list, connection):
    try:
        with connection.cursor() as cursor:

            sql = """SELECT * FROM (
                        SELECT S.`CUI` AS CHILD_CUI, M0.`CUI2` AS PARENT_CUI, C0.`STR` AS PARENT_STR, S.`STY` AS S_TYPE
                        FROM `MRSTY` AS S
                        LEFT OUTER JOIN `MRREL` AS M0 on M0.`CUI1` = S.`CUI`AND M0.`RELA`=%s
                        LEFT OUTER JOIN `MRCONSO` AS C0 on C0.`CUI` = M0.`CUI2` AND C0.`STT`=%s AND C0.`ISPREF`=%s
                        LEFT OUTER JOIN `MRRANK` AS R0 ON (C0.`SAB` = R0.`SAB` AND C0.`TTY` = R0.`TTY`)
                        WHERE S.`CUI` IN %s
                        ORDER BY R0.`RANK` DESC, PARENT_CUI ASC, S.STN ASC) AS TT
                        GROUP BY TT.`CHILD_CUI`"""

            cursor.execute(sql, ('inverse_isa', 'PF', 'Y', cui_list))
            return cursor.fetchall() 

    finally:
        connection.close()

# organise data for input to SQL query, as well as for matching results to PMIDs
def organise(input_data):
    pmids_names = dict()
    cui_list = set() # set ensures no duplicates

    for i in input_data:
        if 'MedlineCitation' in i:
            pmid = i['MedlineCitation']['PMID']
        else:
            pmid = i['PMID']

        concepts = i['concepts']

        for c in concepts:
            cui_list.add(c[1]) # add CUI to set

            # only append PMID as child name should be the same
            if c[1] in pmids_names:
                pmids_names[c[1]]['PMIDs'].append(pmid)
            else:

                # Some key terms have asterisks preceding the name (STR)
                c[0] = c[0].replace('*', '') 

                data = {'child_name': c[0],
                        'PMIDs': [pmid]}
                pmids_names[c[1]] = data

    return {
        'PMIDs': pmids_names, 
        'CUIs': list(cui_list)
    }

def run(input_data):
    print('in umls.run')
    organised = organise(input_data)

    # Connect to the database
    connection = pymysql.connect(host='localhost',
                                 user=config.db_un,
                                 password=config.db_pwd,
                                 db='umls',
                                 charset='utf8mb4',
                                 cursorclass=pymysql.cursors.DictCursor)

    output = execute_sql(organised['CUIs'], connection)

    results = format_json(organised['PMIDs'], group_other(output))
    print('returning from umls.run')
    return results
    