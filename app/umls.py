import pymysql.cursors
import json
import os

import config

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
        dict3['size'] = 500
        
    return dict0

# Arg is a list of CUIs to retrieve information for
def execute_sql(cui_list, connection):
    try:
        with connection.cursor() as cursor:

            sql = """SELECT * FROM (
                        SELECT M0.`CUI1` AS CHILD_CUI, M0.`CUI2` AS PARENT_CUI, C0.`STR` AS PARENT_STR, S.`STY` AS S_TYPE
                        FROM `MRCONSO` AS C0
                        INNER JOIN `MRRANK` AS R0 ON (C0.`SAB` = R0.`SAB` AND C0.`TTY` = R0.`TTY`)
                        INNER JOIN `MRREL` AS M0 on M0.`CUI2` = C0.`CUI`
                        INNER JOIN `MRSTY` AS S on S.`CUI` = M0.`CUI1`
                            WHERE M0.`CUI1` IN %s
                            AND (M0.`RELA`=%s
                            OR M0.`REL`=%s
                            OR M0.`REL`=%s)
                            AND C0.`STT`=%s
                            AND C0.`ISPREF`=%s
                        ORDER BY R0.`RANK` DESC, PARENT_CUI ASC) AS TT
                        GROUP BY TT.`CHILD_CUI`"""

            cursor.execute(sql, (cui_list, 'inverse_isa', 'PAR', 'RB', 'PF', 'Y'))
            return cursor.fetchall() 

    finally:
        connection.close()

# organise data for input to SQL query, as well as for matching results to PMIDs
def organise(input_data):
    pmids_names = dict()
    cui_list = set()

    for i in input_data:
        pmid = i['MedlineCitation']['PMID']
        concepts = i['concepts']

        for c in concepts:
            cui_list.add(c[1]) # add CUI to set

            # only append PMID as child name should be the same
            if c[1] in pmids_names:
                pmids_names[c[1]]['PMIDs'].append(pmid)
            else:

                # Some key terms have asterisks preceding the name (STR)
                c[0] = c[0].replace('*', '') 

                data = {
                    'child_name': c[0],
                    'PMIDs': [
                        pmid
                    ]
                }
                pmids_names[c[1]] = data

    return(pmids_names, list(cui_list))

def q_run(input_data, q):
    q.put(run(input_data))

def run(input_data):
    (pmids_names, cui_list) = organise(input_data)

    # Connect to the database
    connection = pymysql.connect(host='localhost',
                                 user=config.db_un,
                                 password=config.db_pwd,
                                 db='umls',
                                 charset='utf8mb4',
                                 cursorclass=pymysql.cursors.DictCursor)

    output = execute_sql(cui_list, connection)

    results = format_json(pmids_names, output)

    return results