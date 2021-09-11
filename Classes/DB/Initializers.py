from sqlite3.dbapi2 import Cursor


class Initializers:
    """Отвечает за создание таблиц

    Attributes:
        cursor: cursor, через который создаются таблицы"""

    cursor: Cursor

    def __init__(self, cursor: Cursor) -> None:
        self.cursor = cursor

    def createAllTables(self):
        self.createSequenceTable()
        self.createPeptideRowTable()
        self.createPeptideAccessionTable()
        self.createAccessionTable()
        self.createExclusionTable()

    def createSequenceTable(self) -> None:
        self.cursor.execute(
            """CREATE TABLE sequence (
                id INTEGER NOT NULL PRIMARY KEY AUTOINCREMENT,
                accession TEXT NOT NULL UNIQUE,
                description TEXT,
                sequence TEXT NOT NULL
            );"""
        )

    def createAccessionTable(self) -> None:
        self.cursor.execute(
            """CREATE TABLE accession_table (
                id INTEGER NOT NULL PRIMARY KEY AUTOINCREMENT,
                table_number TEXT NOT NULL,
                accession TEXT NOT NULL,
                confidence FLOAT NOT NULL,
                score FLOAT NOT NULL,
                precursor_signal FLOAT NOT NULL,
                sequence TEXT NOT NULL
            );"""
        )

    def createPeptideRowTable(self) -> None:
        self.cursor.execute(
            """CREATE TABLE peptide_row (
                id INTEGER NOT NULL PRIMARY KEY AUTOINCREMENT,
                table_number TEXT NOT NULL,
                confidence FLOAT NOT NULL,
                score FLOAT NOT NULL,
                precursor_signal FLOAT NOT NULL,
                sequence TEXT NOT NULL
            );"""
        )

    def createPeptideAccessionTable(self) -> None:
        self.cursor.execute(
            """CREATE TABLE peptide_accession (
                id INTEGER NOT NULL PRIMARY KEY AUTOINCREMENT,
                row_id INTEGER NOT NULL,
                accession TEXT NOT NULL,
                FOREIGN KEY (row_id)
                    REFERENCES peptide_row (id)
                    ON DELETE CASCADE
            );"""
        )

    def createExclusionTable(self) -> None:
        self.cursor.execute(
            """CREATE TABLE exclusion (
                id INTEGER NOT NULL PRIMARY KEY AUTOINCREMENT,
                accession TEXT NOT NULL UNIQUE
            );"""
        )
