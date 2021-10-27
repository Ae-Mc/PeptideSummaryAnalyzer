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
            DELETE FROM peptide_row WHERE id NOT IN (
                SELECT DISTINCT id FROM peptide_joint
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
      SELECT peptide_row.id AS id, table_number
      FROM peptide_row LEFT JOIN peptide_accession
        ON row_id = peptide_row.id
      WHERE row_id IS NULL
      ORDER BY peptide_row.id;"""
        ).fetchall():
            self.cursor.execute(
                """--sql
                DELETE FROM peptide_row WHERE id >= (?) AND table_number = (?);
                """,
                [row_id, table_num],
            )
