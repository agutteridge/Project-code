import pymysql.cursors

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

    try:
        with connection.cursor() as cursor:
            sql = """SELECT DISTINCT COALESCE(
                        (SELECT C.`STR` 
                         FROM `MRCONSO` AS C
                         WHERE C.`CUI`=M.`CUI`
                         AND C.`SAB`=%s
                         AND C.`TTY`=%s 
                         LIMIT 1),
                        
                        (SELECT C.`STR` 
                         FROM `MRCONSO` AS C
                         WHERE C.`CUI`=M.`CUI` 
                         AND C.`SAB`=%s 
                         AND C.`TTY`=%s 
                         LIMIT 1),

                        (SELECT C.`STR` 
                         FROM `MRCONSO` AS C
                         WHERE C.`CUI`=M.`CUI` 
                         LIMIT 1)) `Concept name`, M.`CUI`
                    FROM `MRCONSO` AS M
                    WHERE M.`CUI` IN %s
                    AND M.`ISPREF`=%s
                    AND M.`STT`=%s"""
            print(sql)
            cursor.execute(sql, ('MSH', 'MH',
                                 'MTH', 'PN',
                                 cui_list, 'Y', 'PF'))
            result = cursor.fetchall()
            print(result) # turn into list extract best STR based on SAB and TTY.. then get PAR? or get PAR based on CUI?
    finally:
        connection.close()  

def run():
    # Connect to the database
    connection = pymysql.connect(host='localhost',
                                 user=config.db_un,
                                 password=config.db_pwd,
                                 db='umls',
                                 charset='utf8mb4',
                                 cursorclass=pymysql.cursors.DictCursor)

    example_id = 'C0027051'

    try:
        with connection.cursor() as cursor:
            # Read a single record
            sql = """SELECT C.`STR` FROM `MRCONSO` AS C WHERE `ISPREF`=%s AND `STT`=%s AND `CUI`=%s LIMIT 1"""
            cursor.execute(sql, ('Y', 'PF', example_id))
            result = cursor.fetchone()
            print(result)
    finally:
        connection.close()