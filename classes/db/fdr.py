"""См. класс FDR."""

from sqlite3.dbapi2 import Cursor


class FDR:
    """Отвечает за работу FDR фильтра.

    Attributes:
        cursor: экземляр класса Cursor для связи с БД.
    """

    cursor: Cursor

    def __init__(self, cursor: Cursor) -> None:
        self.cursor = cursor

    def none(self) -> None:
        """Удаляет все перевёрнутые accession."""
        self.cursor.execute(
            """--sql
            DELETE FROM peptide_accession WHERE IS_REVERSED(accession);
            """
        )
        self.cursor.execute(
            """--sql
            DELETE FROM peptide_row WHERE row_id NOT IN (
                SELECT DISTINCT row_id FROM peptide_joint
            );
            """
        )

    def default(self) -> None:
        """Применяет default FDR фильтр."""

        self.cursor.execute(
            """--sql
            DELETE FROM peptide_accession WHERE IS_REVERSED(accession);"""
        )
        for row_id, table_num in self.cursor.execute(
            """--sql
            SELECT peptide_row.row_id, table_number
            FROM peptide_row LEFT JOIN peptide_accession
                ON peptide_row.row_id = peptide_accession.row_id
            WHERE peptide_accession.row_id IS NULL
            ORDER BY peptide_row.row_id;"""
        ).fetchall():
            self.cursor.execute(
                """--sql
                DELETE FROM peptide_row
                WHERE row_id >= (?) AND table_number = (?);
                """,
                [row_id, table_num],
            )

    def fill_main_accession(self) -> None:
        """Заполняет столбец main_accession в таблице protein_group"""

        self.cursor.execute(
            """--sql
            UPDATE protein_group
            SET main_accession = (
                SELECT accession FROM (
                    SELECT
                        group_id,
                        accession,
                        row_number() OVER (
                            PARTITION BY group_id
                            ORDER BY unused DESC, IS_REVERSED(accession) DESC
                        ) AS row_priority
                    FROM protein_row
                ) pr
                WHERE (
                    row_priority = 1
                    AND protein_group.group_id = pr.group_id
                )
            );
            """
        )
