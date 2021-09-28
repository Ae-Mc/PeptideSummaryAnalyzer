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
        """Clarification required"""
        if confidence.op is not None:
            self.cursor.execute(
                # Clarification required
                #
                # f"""CREATE TABLE filtred_peptide_accession AS SELECT *
                #     FROM peptide_accession t3
                #          LEFT JOIN peptide_row t4
                #          ON t3.row_id = t4.id
                #     WHERE EXISTS (
                #         SELECT table_number, accession
                #         FROM peptide_accession t1
                #              LEFT JOIN peptide_row t2
                #              ON t1.row_id = t2.id
                #         WHERE (
                #             confidence {confidence.op} {confidence.val}
                #             AND t2.table_number = t4.table_number
                #             AND t3.accession = t1.accession
                #         )
                #     );"""
                f"""CREATE TABLE filtred_peptide_accession AS SELECT *
                FROM peptide_accession t3
                    LEFT JOIN peptide_row t4
                    ON t3.row_id = t4.id
                WHERE confidence {confidence.op} {confidence.val};"""
            )
        else:
            self.cursor.execute(
                """CREATE TABLE filtred_peptide_accession AS SELECT *
                    FROM peptide_accession t3
                        LEFT JOIN peptide_row t4
                        ON t3.row_id = t4.id;"""
            )

    def createCountTable(self):
        self.cursor.execute(
            """CREATE TABLE accession_count_per_table AS
            SELECT table_number, accession, COUNT(*) AS count_per_table
            FROM filtred_peptide_accession
            GROUP BY table_number, accession
            ORDER BY table_number, accession;"""
        )

    def createGeneralCountTable(self):
        self.cursor.execute(
            """CREATE TABLE general_accession_count AS
            SELECT accession, COUNT(*) AS general_count
            FROM filtred_peptide_accession
            GROUP BY accession
            ORDER BY accession;"""
        )

    def applyProteinGrouping(self) -> None:
        self.cursor.execute(
            """INSERT INTO representative (row_id, table_number, representative)
            SELECT id as row_id, table_number, accession FROM (
                SELECT
                    joint.id,
                    joint.accession,
                    joint.table_number,
                    row_number() OVER (
                        PARTITION BY joint.id
                        ORDER BY
                            count_per_table DESC,
                            general_count DESC,
                            LENGTH(sequence.sequence) DESC,
                            joint.id DESC
                    ) AS seqnum
                FROM
                    peptide_joint joint
                    LEFT JOIN accession_count_per_table t_count
                        ON t_count.table_number = joint.table_number
                        AND t_count.accession = joint.accession
                    LEFT JOIN general_accession_count g_count
                        ON g_count.accession = joint.accession
                    LEFT JOIN sequence ON joint.accession = sequence.accession
            ) WHERE seqnum = 1;"""
        )
        self.cursor.execute(
            """INSERT INTO accession_group (
                representative_id,
                accession,
                count_in_table
            ) SELECT repr.id, joint.accession, count_per_table
            FROM peptide_joint joint
                INNER JOIN representative repr ON joint.id = repr.row_id
                INNER JOIN accession_count_per_table t_count
                        ON t_count.table_number = joint.table_number
                        AND t_count.accession = joint.accession;"""
        )

    def fillPeptideTable(self) -> None:
        self.cursor.execute(
            """INSERT INTO peptide (
                table_number,
                accession,
                confidence,
                score,
                peptide_intensity,
                sequence
            ) SELECT
                row.table_number,
                representative,
                confidence,
                score,
                peptide_intensity,
                sequence
            FROM peptide_row row JOIN representative repr ON row.id = row_id"""
        )
