import pymysql.cursors

import config

def run_with_data(paper_terms):
    # Connect to the database
    connection = pymysql.connect(host='localhost',
                                 user=config.db_un,
                                 password=config.db_pwd,
                                 db='umls',
                                 charset='utf8mb4',
                                 cursorclass=pymysql.cursors.DictCursor)

    try:
        with connection.cursor() as cursor:
            for p in paper_terms:
                concepts = p['concepts']
                for c in concepts:
                    cui = c[1]
                    print(cui)
                    sql = "SELECT `STR`, `SAB`, `TTY` FROM `MRCONSO` WHERE `ISPREF`=%s AND `STT`=%s AND `CUI`=%s"
                    cursor.execute(sql, ('Y', 'PF', cui))
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
            sql = "SELECT `STR` FROM `MRCONSO` WHERE `ISPREF`=%s AND `STT`=%s AND `CUI`=%s LIMIT 1"
            cursor.execute(sql, ('Y', 'PF', example_id))
            result = cursor.fetchone()
            print(result)
    finally:
        connection.close()