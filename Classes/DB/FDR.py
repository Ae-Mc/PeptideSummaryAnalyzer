from sqlite3.dbapi2 import Cursor


class FDR:
    cursor: Cursor

    def __init__(self, cursor: Cursor) -> None:
        self.cursor = cursor

    def default(self) -> None:
        self.cursor.execute(
            """DELETE FROM peptide_accession WHERE IS_REVERSED(accession);"""
        )
        for rowID, tableNum in self.cursor.execute(
            """SELECT peptide_row.id AS id, table_number
            FROM peptide_row LEFT JOIN peptide_accession
                ON row_id = peptide_row.id
            WHERE row_id IS NULL;"""
        ).fetchall():
            self.cursor.execute(
                f"""DELETE FROM peptide_row
                WHERE id >= {rowID} AND table_number = {tableNum};"""
            )