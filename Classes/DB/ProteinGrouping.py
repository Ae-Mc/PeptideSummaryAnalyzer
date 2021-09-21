from Classes.Comparable import Comparable
from sqlite3 import Cursor


class ProteinGrouping:
    """Отвечает за работу Protein Grouping

    Attributes:
        cursor: экземляр класса Cursor для связи с базой данных
    """

    cursor: Cursor

    def __init__(self, cursor: Cursor) -> None:
        self.cursor = cursor

    def createFiltredAccessionTable(self, confidence: Comparable):
        if confidence.op is not None:
            self.cursor.execute(
                f"""CREATE TABLE filtred_peptide_accession AS SELECT *
                    FROM peptide_accession t3
                         LEFT JOIN peptide_row t4
                         ON t3.row_id = t4.id
                    WHERE EXISTS (
                        SELECT table_number, accession
                        FROM peptide_accession t1 
                             LEFT JOIN peptide_row t2
                             ON t1.row_id = t2.id
                        WHERE (
                            confidence {confidence.op} {confidence.val}
                            AND t2.table_number = t4.table_number
                            AND t3.accession = t1.accession
                        )
                        GROUP BY table_number, accession
                    )
                );"""
            )
        else:
            self.cursor.execute(
                """CREATE TABLE filtred_peptide_accession AS SELECT *
                    FROM peptide_accession t3
                        LEFT JOIN peptide_row t4
                        ON t3.row_id = t4.id;"""
            )

    def createCountTable(self):
        """CREATE TABLE accession_count_per_table AS
        SELECT table_number, accession, COUNT(*) AS count_per_table
        FROM filtred_peptide_accession
        GROUP BY table_number, accession
        ORDER BY table_number, accession;"""

    def createGeneralCountTable(self):
        """CREATE TABLE general_accession_count AS
        SELECT accession, COUNT(*) AS general_count
        FROM filtred_peptide_accession
        GROUP BY accession
        ORDER BY accession;"""
