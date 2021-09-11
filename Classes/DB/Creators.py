from sqlite3.dbapi2 import Cursor


class Creators:
    """Отвечает за создание таблиц

    Attributes:
        cursor: cursor, через который создаются таблицы"""

    cursor: Cursor

    def __init__(self, cursor: Cursor) -> None:
        self.cursor = cursor

    def createAllTables(self):
        self.createSequenceTable()
        self.createExclusionTable()
        self.createPeptideRowTable()
        self.createPeptideAccessionTable()
        self.createPeptideTable()

    def createSequenceTable(self) -> None:
        self.cursor.execute(
            """CREATE TABLE sequence (
                id INTEGER NOT NULL PRIMARY KEY AUTOINCREMENT,
                accession TEXT NOT NULL UNIQUE,
                description TEXT,
                sequence TEXT NOT NULL
            );"""
        )

    def createExclusionTable(self) -> None:
        self.cursor.execute(
            """CREATE TABLE exclusion (
                id INTEGER NOT NULL PRIMARY KEY AUTOINCREMENT,
                accession TEXT NOT NULL UNIQUE
            );"""
        )

    def createPeptideRowTable(self) -> None:
        self.cursor.execute(
            """CREATE TABLE peptide_row (
                id INTEGER NOT NULL PRIMARY KEY AUTOINCREMENT,
                table_number TEXT NOT NULL,
                confidence FLOAT NOT NULL,
                score FLOAT NOT NULL,
                peptide_intensity FLOAT NOT NULL,
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

    def createPeptideTable(self) -> None:
        self.cursor.execute(
            """CREATE TABLE peptide (
                id INTEGER NOT NULL PRIMARY KEY AUTOINCREMENT,
                table_number TEXT NOT NULL,
                accession TEXT NOT NULL,
                confidence FLOAT NOT NULL,
                score FLOAT NOT NULL,
                peptide_intensity FLOAT NOT NULL,
                sequence TEXT NOT NULL
            );"""
        )
