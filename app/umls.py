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

            var = 'C0008659'

            sql1 = """SELECT C1.`CUI`, MRREL.`CUI2`
                        FROM `MRCONSO` AS C1
                        INNER JOIN `MRREL` ON MRREL.`CUI1` = C1.`CUI`
                        INNER JOIN `MRCONSO` AS C2 ON C2.`CUI` = MRREL.`CUI2`
                        WHERE (MRREL.`RELA`=%s
                        OR MRREL.`REL`=%s
                        OR MRREL.`REL`=%s)
                        AND C1.`CUI` IN %s"""

            cursor.execute(sql1, ('inverse_isa', 'PAR', 'RB', cui_list))
            result = cursor.fetchall()

            use stored procedure instead recursively call
            either a set number of times, or until a certain stopping point?
            sql = """SELECT C1.`CUI`, C1.`STR`, 
                     FROM (SELECT C2.`CUI`, MAX(R2.`RANK`) AS MAXRANK
                           FROM `MRCONSO` AS C2 
                           INNER JOIN `MRRANK` AS R2 ON (C2.`SAB`=R2.`SAB` AND C2.`TTY`=R2.`TTY`)
                           WHERE C2.`ISPREF`=%s
                           AND C2.`STT`=%s
                           AND C2.`CUI` IN %s
                           GROUP BY C2.`CUI`) AS T
                     INNER JOIN `MRRANK` AS R1 ON R1.`RANK`=T.`MAXRANK`
                     INNER JOIN `MRCONSO` C1 ON (C1.`SAB`=R1.`SAB` AND C1.`TTY`=R1.`TTY`)
                     INNER JOIN `MRREL` M1 ON (M1.CUI1 = )
                     WHERE C1.`CUI` = T.`CUI`"""

            # cursor.execute(sql, ('Y', 'PF', cui_list))
            result = cursor.fetchall()
    finally:
        connection.close()

    # for r in result:
    #     for p in paper_terms:
    #         concepts = p['concepts']
    #         for c in concepts:
    #             if r['CUI'] == c[1]:
    #                 r.setdefault('PMID', []).append(p['PMID'])
    return result