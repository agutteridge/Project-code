import pymysql.cursors
import json
import os

import config

def format_results(pmids_names, results):
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
        
    return dict0

# Arg is a list of CUIs to retrieve information for
def execute_sql(cui_list, connection):
    try:
        with connection.cursor() as cursor:

            sql = """SELECT DISTINCT T2.`CHILD_CUI`, T2.`S_TYPE`, T2.`PARENT_CUI`, C2.`STR` AS PARENT_STR
                      FROM `MRCONSO` AS C2
                      INNER JOIN (SELECT DISTINCT T1.`CHILD_CUI`, T1.`PARENT_CUI`, MAX(R1.`RANK`) AS MAXRANK, T1.`S_TYPE`
                                  FROM `MRCONSO` AS C1
                                  INNER JOIN `MRRANK` AS R1 ON (C1.`SAB` = R1.`SAB` AND C1.`TTY` = R1.`TTY`)
                                  INNER JOIN (SELECT DISTINCT M1.`CUI1` AS CHILD_CUI, M1.`CUI2` AS PARENT_CUI, S.`STY` AS S_TYPE
                                              FROM `MRREL` AS M1
                                              INNER JOIN `MRSTY` AS S ON S.`CUI` = M1.`CUI1`
                                              WHERE M1.`CUI1` IN %s
                                              AND (M1.`RELA`=%s
                                                   OR M1.`REL`=%s
                                                   OR M1.`REL`=%s)) AS T1 ON T1.`PARENT_CUI` = C1.`CUI` 
                                  WHERE C1.`STT`=%s
                                  AND C1.`ISPREF`=%s
                                  GROUP BY C1.`CUI`) AS T2 ON T2.`PARENT_CUI` = C2.`CUI`
                      INNER JOIN `MRRANK` AS R3 ON (R3.`RANK` = T2.`MAXRANK` AND C2.`SAB` = R3.`SAB` AND C2.`TTY` = R3.`TTY`)
                      WHERE C2.`STT`=%s
                      AND C2.`ISPREF`=%s
                      ORDER BY CHILD_CUI ASC"""

            cursor.execute(sql, (cui_list, 'inverse_isa', 'PAR', 'RB', 'PF', 'Y', 'PF', 'Y'))
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
                data = {
                    'child_name': c[0],
                    'PMIDs': [
                        pmid
                    ]
                }
                pmids_names[c[1]] = data

    return(pmids_names, list(cui_list))    

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

    with open(os.path.join('./tests/resources', 'umls_db_output.json'), 'w') as datafile:
            datafile.write(json.dumps(output))

    with open(os.path.join('./tests/resources', 'umls_output.json'), 'w') as datafile:
        datafile.write(json.dumps(format_results(pmids_names, output)))
