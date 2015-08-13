DELIMITER $$
CREATE PROCEDURE `GETRELS`()
SQL SECURITY INVOKER
BEGIN
SELECT C1.STR AS 'CHILD STR', C1.SAB AS 'CHILD SOURCE', C1.TTY AS 'CHILD TERM TYPE', C1.CUI AS 'CHILD CUI'
FROM MRCONSO C1
WHERE C1.CUI IN ('C0069515', 'C0086418', 'C0086418', 'C1702024', '3027822')
AND C1.ISPREF = 'Y'
AND C1.STT = 'PF'
AND C1.SAB = 'MTH';
END$$
DELIMITER;

export to results table, delete all rows that aren''t needed?


!!! too unspecific !!!
finding Semantic Type:
SELECT TUI, STN, STY, ATUI FROM MRSTY
WHERE CUI = 'C0027051';

+------+------------+---------------------+------------+
| TUI  | STN        | STY                 | ATUI       |
+------+------------+---------------------+------------+
| T047 | B2.2.1.2.1 | Disease or Syndrome | AT32679180 |
+------+------------+---------------------+------------+

finding parents of Semantic Type:
SELECT SRSTRE1.UI3 AS 'Concept TUI', 
	SRDEF.STY_RL AS 'Semantic Type'
FROM SRSTRE1 
	JOIN SRDEF ON SRDEF.UI = SRSTRE1.UI3
WHERE SRSTRE1.UI2 = 'T186'
AND SRSTRE1.UI1 = 'T047';
!!! too unspecific !!!

SELECT C.AUI, C.CUI, C.STR FROM MRCONSO C JOIN MRRANK R on C.SAB=R.SAB
WHERE CUI IN ('C0069515', 'C0086418', 'C0086418', 'C1702024', '3027822')
AND C.AUI IN 
(SELECT C2.AUI FROM MRCONSO C2 JOIN MRRANK R ON C2.SAB=R.SAB
WHERE R.TTY = C.TTY
ORDER BY R.RANK DESC
LIMIT 1);