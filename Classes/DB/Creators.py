from sqlite3.dbapi2 import Cursor


def IsReversed(accession: str) -> bool:
    return accession.startswith("RRRRR")


class Creators:
    """Отвечает за создание таблиц и дополнительных функций

    Attributes:
        cursor: cursor, через который создаются таблицы"""

    cursor: Cursor

    def __init__(self, cursor: Cursor) -> None:
        self.cursor = cursor

    def createAllFunctions(self) -> None:
        self.cursor.connection.create_function("IS_REVERSED", 1, IsReversed)

    def createAllTables(self):
        self.createSequenceTable()
        self.createExclusionTable()
        self.createPeptideRowTable()
        self.createPeptideAccessionTable()
        self.createPeptideTable()
        self.createRepresentativeTable()
        self.createGroupTable()
        self.createPeptideWithSum()

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

    def createRepresentativeTable(self) -> None:
        self.cursor.execute(
            """CREATE TABLE representative (
                id INTEGER NOT NULL PRIMARY KEY AUTOINCREMENT,
                table_number TEXT NOT NULL,
                accession TEXT NOT NULL
            );"""
        )

    def createGroupTable(self) -> None:
        self.cursor.execute(
            """CREATE TABLE accession_group (
                id INTEGER NOT NULL PRIMARY KEY AUTOINCREMENT,
                representative_id INTEGER NOT NULL,
                accession TEXT NOT NULL,
                count_in_table INTEGER NOT NULL
            );"""
        )

    def createPeptideWithSum(self) -> None:
        self.cursor.execute(
            """CREATE TABLE peptide_with_sum (
                id INTEGER NOT NULL PRIMARY KEY AUTOINCREMENT,
                table_number TEXT NOT NULL,
                accession TEXT NOT NULL,
                sc_sum FLOAT NOT NULL,
                peptide_intensity_sum FLOAT NOT NULL,
                sc_norm FLOAT NOT NULL,
                peptide_intensity_norm FLOAT NOT NULL,
                sc_norm_to_file_norm_ratio FLOAT NOT NULL,
                peptide_intensity_norm_to_file_norm_ratio FLOAT NOT NULL
            );"""
        )
