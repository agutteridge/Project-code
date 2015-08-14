import pymysql.cursors
import json

import config

def run_with_data(paper_terms):
    cui_list = []
    for p in paper_terms:
        concepts = p['concepts']
        for c in concepts:
            cui_list.append(c[1])

    # Connect to the database
    connection = pymysql.connect(host='localhost',
                                 user=config.db_un,
                                 password=config.db_pwd,
                                 db='umls',
                                 charset='utf8mb4',
                                 cursorclass=pymysql.cursors.DictCursor)

    # result = '' # JSON..?
    try:
        with connection.cursor() as cursor:

            sql = """SELECT DISTINCT T2.`CHILD_CUI`, T2.`S_TYPE`, T2.`PARENT1_CUI`, C2.`STR` AS PARENT1_STR
                      FROM `MRCONSO` AS C2
                      INNER JOIN (SELECT DISTINCT T1.`CHILD_CUI`, T1.`PARENT1_CUI`, MAX(R1.`RANK`) AS MAXRANK, T1.`S_TYPE`
                                  FROM `MRCONSO` AS C1
                                  INNER JOIN `MRRANK` AS R1 ON (C1.`SAB` = R1.`SAB` AND C1.`TTY` = R1.`TTY`)
                                  INNER JOIN (SELECT DISTINCT M1.`CUI1` AS CHILD_CUI, M1.`CUI2` AS PARENT1_CUI, S.`STY` AS S_TYPE
                                              FROM `MRREL` AS M1
                                              INNER JOIN `MRSTY` AS S ON S.`CUI` = M1.`CUI1`
                                              WHERE M1.`CUI1` IN %s
                                              AND (M1.`RELA`=%s
                                                   OR M1.`REL`=%s
                                                   OR M1.`REL`=%s)) AS T1 ON T1.`PARENT1_CUI` = C1.`CUI` 
                                  WHERE C1.`STT`=%s
                                  AND C1.`ISPREF`=%s
                                  GROUP BY C1.`CUI`) AS T2 ON T2.`PARENT1_CUI` = C2.`CUI`
                      INNER JOIN `MRRANK` AS R3 ON (R3.`RANK` = T2.`MAXRANK` AND C2.`SAB` = R3.`SAB` AND C2.`TTY` = R3.`TTY`)
                      WHERE C2.`STT`=%s
                      AND C2.`ISPREF`=%s"""

            cursor.execute(sql, (cui_list, 'inverse_isa', 'PAR', 'RB', 'PF', 'Y', 'PF', 'Y'))
            results = cursor.fetchall() #parent IDs

    finally:
        connection.close()

    # for r in result:
    #     for p in paper_terms:
    #         concepts = p['concepts']
    #         for c in concepts:
    #             if r['CUI'] == c[1]:
    #                 r.setdefault('PMID', []).append(p['PMID'])
    return results