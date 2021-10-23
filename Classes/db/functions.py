"""Отвечает за различные функции, которые нельзя объединить в одну общую
группу."""
from sqlite3.dbapi2 import Cursor

from Classes.comparable import Comparable


class Functions:
    """Отвечает за разные функции, которые нельзя объединить в одну общую группу.

    Attributes:
        cursor: экземляр класса Cursor для связи с базой данных.
    """

    cursor: Cursor

    def __init__(self, cursor: Cursor) -> None:
        self.cursor = cursor

    def test_fasta_database(self) -> None:
        """Проверка, все ли accession из peptide_accession присутствуют в .fasta БД.

        В случае, если какие-то accession не найдены, вызывает IndexError,
        выводя все отсутствующие accession.

        Raises:
            IndexError: Ошибка при поиске некоторых accession."""

        absence_accessions = self.cursor.execute(
            """--sql
            SELECT DISTINCT accession FROM peptide_accession
            WHERE (
                NOT IS_REVERSED(accession)
                AND accession NOT IN (SELECT accession FROM sequence)
            );"""
        ).fetchall()
        if len(absence_accessions) > 0:
            raise IndexError(
                "Accessions not found in .fasta file:\n\t"
                + "\n\t".join(map(lambda x: x[0], absence_accessions))
            )

    def apply_exclusion(self) -> None:
        """Применение ID exclusion фильтра."""

        if (
            self.cursor.execute(
                "SELECT COUNT(*) FROM exclusion LIMIT(1)"
            ).fetchall()[0][0]
            > 0
        ):
            self.cursor.execute(
                """--sql
                DELETE FROM peptide_accession
                WHERE accession IN (SELECT accession FROM exclusion);"""
            )
            self.cursor.execute(
                """--sql
                DELETE FROM peptide_row
                WHERE id IN (SELECT peptide_row.id AS id
                    FROM peptide_row LEFT JOIN peptide_accession
                        ON row_id = peptide_row.id
                    WHERE row_id IS NULL);"""
            )

    def apply_peptide_confidence_value(self, confidence: Comparable) -> None:
        """Применение #Protein filter -> Peptide Confidence (value).

        Args:
            confidence: параметр confidence для сравнения."""

        if confidence.operation is not None:
            self.cursor.execute(
                f"""--sql
                DELETE FROM peptide_row WHERE id NOT IN (
                    SELECT DISTINCT pr.id
                    FROM peptide_accession pa JOIN peptide_row pr
                        ON pa.row_id = pr.id
                    WHERE (table_number, accession) IN (
                        SELECT DISTINCT table_number, accession
                        FROM peptide_joint
                        WHERE confidence {confidence.operation}{confidence.val}
                    )
                );
                """
            )

    def apply_peptide_confidence_default(self) -> None:
        """Применение #Protein filter -> Peptide Confidence (default)."""

        self.cursor.execute(
            """
            --sql
            DELETE FROM peptide_row WHERE id NOT IN (
                SELECT DISTINCT id
                FROM peptide_joint WHERE confidence >= 95
                GROUP BY table_number, accession
                HAVING COUNT(*) >= 2 OR MAX(confidence) >= 99
            );
            """
        )

    def apply_peptide_confidence_value2(self, confidence: Comparable) -> None:
        """Применение #Peptide filter -> Peptide Confidence (value)

        Args:
            confidence: параметр confidence для сравнения."""

        if confidence.operation is not None:
            self.cursor.execute(
                f"""--sql
                DELETE FROM peptide
                WHERE NOT (confidence {confidence.operation} {confidence.val});
                """
            )
