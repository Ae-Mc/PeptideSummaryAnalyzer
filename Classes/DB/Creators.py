from sqlite3.dbapi2 import Cursor
from typing import Optional


def IsReversed(accession: str) -> bool:
    return accession.startswith("RRRRR")


def GroupNumber(table_number: str) -> int:
    splitted = table_number.split(".")
    return int(splitted[0])


class Creators:
    """Отвечает за создание таблиц и дополнительных функций

    Attributes:
        cursor: cursor, через который создаются таблицы"""

    cursor: Cursor

    def __init__(self, cursor: Cursor) -> None:
        self.cursor = cursor

    def createAllFunctions(self) -> None:
        self.cursor.connection.create_function("IS_REVERSED", 1, IsReversed)
        self.cursor.connection.create_function("GET_GROUP_NUMBER", 1, GroupNumber)

    def createAllTables(self):
        self.createSequenceTable()
        self.createExclusionTable()
        self.createPeptideRowTable()
        self.createPeptideAccessionTable()
        self.createPeptideTable()
        self.createRepresentativeTable()
        self.createGroupTable()
        self.createPeptideWithSum()
        self.createJointPeptideTableView()

    def createSequenceTable(self) -> None:
        self.cursor.execute(
            """--sql
            CREATE TABLE sequence (
                id INTEGER NOT NULL PRIMARY KEY AUTOINCREMENT,
                accession TEXT NOT NULL UNIQUE,
                description TEXT,
                sequence TEXT NOT NULL,
                raw_sequence TEXT NOT NULL
            );"""
        )

    def createExclusionTable(self) -> None:
        self.cursor.execute(
            """--sql
            CREATE TABLE exclusion (
                id INTEGER NOT NULL PRIMARY KEY AUTOINCREMENT,
                accession TEXT NOT NULL UNIQUE
            );"""
        )

    def createPeptideRowTable(self) -> None:
        self.cursor.execute(
            """--sql
            CREATE TABLE peptide_row (
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
            """--sql
            CREATE TABLE peptide_accession (
                id INTEGER NOT NULL PRIMARY KEY AUTOINCREMENT,
                row_id INTEGER NOT NULL,
                accession TEXT NOT NULL,
                FOREIGN KEY (row_id) REFERENCES peptide_row (id) ON DELETE CASCADE
            );"""
        )

    def createJointPeptideTableView(self) -> None:
        self.cursor.execute(
            """--sql
            CREATE VIEW peptide_joint AS
            SELECT row.*, p_acc.accession
            FROM peptide_row row INNER JOIN peptide_accession p_acc
                ON row.id = p_acc.row_id;"""
        )

    def createPeptideTable(self) -> None:
        self.cursor.execute(
            """--sql
            CREATE TABLE peptide (
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
            """--sql
            CREATE TABLE representative (
                id INTEGER NOT NULL PRIMARY KEY AUTOINCREMENT,
                row_id INTEGER NOT NULL,
                table_number TEXT NOT NULL,
                representative TEXT NOT NULL,
                FOREIGN KEY (row_id) REFERENCES peptide_row (id) ON DELETE CASCADE
            );"""
        )

    def createGroupTable(self) -> None:
        self.cursor.execute(
            """--sql
            CREATE TABLE accession_group (
                id INTEGER NOT NULL PRIMARY KEY AUTOINCREMENT,
                representative_id INTEGER NOT NULL,
                accession TEXT NOT NULL,
                count_in_table INTEGER NOT NULL
            );"""
        )

    def createPeptideWithSum(self) -> None:
        self.cursor.execute(
            """--sql
            CREATE TABLE peptide_with_sum (
                id INTEGER NOT NULL PRIMARY KEY AUTOINCREMENT,
                table_number TEXT NOT NULL,
                accession TEXT NOT NULL,
                count INT NOT NULL,
                seq_length_sum INT NOT NULL,
                sc_sum FLOAT NOT NULL,
                peptide_intensity_sum FLOAT NOT NULL,
                sc_norm FLOAT NOT NULL,
                peptide_intensity_norm FLOAT NOT NULL,
                sc_norm_to_file_norm_ratio FLOAT NOT NULL,
                peptide_intensity_norm_to_file_norm_ratio FLOAT NOT NULL
            );"""
        )
