"""Отвечает за создание таблиц, которые используются в БД, за исключением
внутренних таблиц различных этапов обработки данных и за создание функций,
доступных в запросах SQL."""

from sqlite3.dbapi2 import Cursor


def is_reversed(accession: str) -> bool:
    """Проверяет, является ли accession перевёрнутым

    Args:
        accession (str): проверяемый accession

    Returns:
        bool: true, если перевёрнутый, иначе - false
    """
    return accession.startswith("RRRRR")


def group_number(table_number: str) -> int:
    """Ковертирует номер таблицы в номер группы, в которой эта таблица состоит.

    Например, 0.1 превретится в 0, 11.3 в 11 и т. д.

    Args:
        table_number (str): номер таблицы

    Returns:
        int: номер группы
    """
    splitted = table_number.split(".")
    return int(splitted[0])


class Creators:
    """Отвечает за создание таблиц и дополнительных функций

    Attributes:
        cursor: cursor, через который создаются таблицы"""

    cursor: Cursor

    def __init__(self, cursor: Cursor) -> None:
        self.cursor = cursor

    def create_all_functions(self) -> None:
        """Создаёт дополнительные функции, которые можно использовать в SQL
        запросах.
        """
        self.cursor.connection.create_function("IS_REVERSED", 1, is_reversed)
        self.cursor.connection.create_function(
            "GET_GROUP_NUMBER", 1, group_number
        )

    def create_all_tables(self):
        """Создаёт все таблицы, используемые программными модулями при передаче
        данных между собой и для считывания входных данных.
        """
        self._create_sequence_table()
        self._create_exclusion_table()
        self._create_peptide_row_table()
        self._create_peptide_accession_table()
        self._create_accession_count_per_table_table()
        self._create_peptide_table()
        self._create_representative_table()
        self._create_group_table()
        self._create_peptide_with_sum()
        self._create_joint_peptide_table_view()

    def _create_sequence_table(self) -> None:
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

    def _create_exclusion_table(self) -> None:
        self.cursor.execute(
            """--sql
            CREATE TABLE exclusion (
                id INTEGER NOT NULL PRIMARY KEY AUTOINCREMENT,
                accession TEXT NOT NULL UNIQUE
            );"""
        )

    def _create_peptide_row_table(self) -> None:
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

    def _create_peptide_accession_table(self) -> None:
        self.cursor.execute(
            """--sql
            CREATE TABLE peptide_accession (
                id INTEGER NOT NULL PRIMARY KEY AUTOINCREMENT,
                row_id INTEGER NOT NULL,
                accession TEXT NOT NULL,
                FOREIGN KEY (row_id) REFERENCES peptide_row (id)
                    ON DELETE CASCADE
            );"""
        )

    def _create_accession_count_per_table_table(self) -> None:
        self.cursor.execute(
            """
            --sql
            CREATE TABLE accession_count_per_table (
                id INTEGER NOT NULL PRIMARY KEY AUTOINCREMENT,
                table_number TEXT NOT NULL,
                accession TEXT NOT NULL,
                count INTEGER NOT NULL,
                UNIQUE (table_number, accession) ON CONFLICT FAIL
            );
            """
        )

    def _create_joint_peptide_table_view(self) -> None:
        self.cursor.execute(
            """--sql
            CREATE VIEW peptide_joint AS
            SELECT row.*, p_acc.accession
            FROM peptide_row row INNER JOIN peptide_accession p_acc
                ON row.id = p_acc.row_id;"""
        )

    def _create_peptide_table(self) -> None:
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

    def _create_representative_table(self) -> None:
        self.cursor.execute(
            """--sql
            CREATE TABLE representative (
                representative_id INTEGER NOT NULL PRIMARY KEY AUTOINCREMENT,
                representative TEXT NULL
            );"""
        )

    def _create_group_table(self) -> None:
        self.cursor.execute(
            """--sql
            CREATE TABLE accession_group (
                id INTEGER NOT NULL PRIMARY KEY AUTOINCREMENT,
                representative_id INTEGER NOT NULL,
                accession TEXT NOT NULL UNIQUE
            );"""
        )

    def _create_peptide_with_sum(self) -> None:
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
