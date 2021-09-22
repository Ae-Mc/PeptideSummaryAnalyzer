from Classes.Comparable import Comparable
from sqlite3.dbapi2 import Cursor


class Functions:
    """Отвечает за разные функции, которые нельзя объединить в одну общую группу

    Attributes:
        cursor: экземляр класса Cursor для связи с базой данных
    """

    cursor: Cursor

    def __init__(self, cursor: Cursor) -> None:
        self.cursor = cursor

    def testFastaDatabase(self) -> None:
        """Проверка, все ли Accession из peptide_accession присутствуют в .fasta БД.

        В случае, если какие-то Accession не найдены, вызывает IndexError, выводя все
        ненайденные Accession"""

        absenceAccessions = self.cursor.execute(
            r"""SELECT DISTINCT accession FROM peptide_accession
            WHERE (
                NOT IS_REVERSED(accession)
                AND accession NOT IN (SELECT accession FROM sequence)
            );"""
        ).fetchall()
        if len(absenceAccessions) > 0:
            raise IndexError(
                "Accessions not found in .fasta file:\n\t"
                + "\n\t".join(map(lambda x: x[0], absenceAccessions))
            )

    def applyExclusion(self) -> None:
        """Применение ID exclusion фильтра"""

        if (
            self.cursor.execute("SELECT COUNT(*) FROM exclusion LIMIT(1)").fetchall()[
                0
            ][0]
            > 0
        ):
            self.cursor.execute(
                """DELETE FROM peptide_accession
                WHERE accession IN (SELECT accession FROM exclusion);"""
            )
            self.cursor.execute(
                """DELETE FROM peptide_row
                WHERE id IN (SELECT peptide_row.id AS id
                    FROM peptide_row LEFT JOIN peptide_accession
                        ON row_id = peptide_row.id
                    WHERE row_id IS NULL);"""
            )

    def applyPeptideConfidenceValue(self, confidence: Comparable) -> None:
        if confidence.op is not None:
            self.cursor.execute(
                f"""DELETE FROM peptide_accession WHERE id IN (
                    SELECT t3.id
                    FROM peptide_accession t3
                         LEFT JOIN peptide_row t4
                         ON t3.row_id = t4.id
                    WHERE NOT EXISTS (
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
                );
                """
            )
            self.removeLeftoversFromPeptideRow()

    def applyPeptideConfidenceDefault(self) -> None:
        self.cursor.execute(
            """DELETE FROM peptide_accession AS paout WHERE id IN (
                SELECT DISTINCT pain.id
                FROM peptide_accession pain JOIN peptide_row pr ON pain.row_id = pr.id
                WHERE confidence < 99
                GROUP BY table_number, accession
                HAVING COUNT(*) < 2 OR MAX(confidence) < 95
            );"""
        )
        self.removeLeftoversFromPeptideRow()

    def removeLeftoversFromPeptideRow(self) -> None:
        self.cursor.execute(
            """DELETE FROM peptide_row
            WHERE id NOT IN (SELECT DISTINCT row_id FROM peptide_accession);"""
        )
