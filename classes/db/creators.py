"""Отвечает за создание таблиц, которые используются в БД, за исключением
внутренних таблиц различных этапов обработки данных и за создание функций,
доступных в запросах SQL."""

from sqlite3.dbapi2 import Cursor
from math import log, exp
from typing import Optional


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


# pylint: disable=invalid-name
def ln(number: Optional[float]) -> Optional[float]:
    """Natural logarithm that accepts NULL as argument"""
    if number is None:
        return None
    return log(number)


def my_round(number: Optional[float], precision: int) -> Optional[str]:
    """Custom round. Converts float to string with custom precision."""
    if number is None:
        return None
    return (
        ("{:." + str(precision) + "f}").format(number).rstrip("0").rstrip(".")
    )


def exponent(number: Optional[float]) -> Optional[float]:
    """Returns e^number."""
    if number is None:
        return None
    return exp(number)


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
        self.cursor.connection.create_function("LN", 1, ln, deterministic=True)
        self.cursor.connection.create_function(
            "EXP", 1, exponent, deterministic=True
        )
        self.cursor.connection.create_function(
            "MY_ROUND", 2, my_round, deterministic=True
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
        self._create_protein_group_table()
        self._create_protein_row_table()
        self._create_fdr_data_table()
        self._create_fdr_summary_table()

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
                row_id INTEGER NOT NULL PRIMARY KEY AUTOINCREMENT,
                table_number TEXT NOT NULL,
                N INTEGER NOT NULL,
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
                FOREIGN KEY (row_id) REFERENCES peptide_row (row_id)
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
                SELECT peptide_row.*, peptide_accession.accession
                FROM peptide_row INNER JOIN peptide_accession USING(row_id);"""
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

    def _create_protein_group_table(self) -> None:
        self.cursor.execute(
            """--sql
            CREATE TABLE protein_group (
                group_id INTEGER NOT NULL PRIMARY KEY AUTOINCREMENT,
                table_number TEXT NOT NULL,
                N INTEGER NOT NULL,
                main_accession TEXT NULL,
                UNIQUE (table_number, N)
            );
            """
        )

    def _create_protein_row_table(self) -> None:
        self.cursor.execute(
            """--sql
            CREATE TABLE protein_row (
                row_id INTEGER NOT NULL PRIMARY KEY AUTOINCREMENT,
                group_id INTEGER NOT NULL,
                accession TEXT NOT NULL,
                unused FLOAT NOT NULL,
                FOREIGN KEY (group_id) REFERENCES protein_group (group_id)
            );
            """
        )

    def _create_fdr_data_table(self) -> None:
        self.cursor.execute(
            """
            --sql
            CREATE TABLE fdr_data (
                fdr_data_id INTEGER NOT NULL PRIMARY KEY AUTOINCREMENT,
                table_number TEXT NOT NULL,
                table_n INTEGER NOT NULL,
                unused FLOAT NOT NULL,
                accession TEXT NOT NULL,
                index_n INTEGER NOT NULL,
                observed_fdr FLOAT NOT NULL,
                UNIQUE (table_number, table_n),
                UNIQUE (table_number, index_n)
            );
            """
        )

    def _create_fdr_summary_table(self) -> None:
        self.cursor.execute(
            """
            --sql
            CREATE TABLE fdr_summary (
                table_number TEXT NOT NULL UNIQUE,
                target_count INTEGER NULL,
                decoy_count INTEGER NULL,
                k FLOAT NULL
                    AS (target_count * 1.0 / decoy_count) STORED,
                a FLOAT NULL,
                b FLOAT NULL,
                fdr_01 INTEGER NULL
                    AS (ROUND(ln(0.1 / a) / b, 0)) STORED,
                fdr_05 INTEGER NULL
                    AS (ROUND(ln(0.5 / a) / b, 0)) STORED,
                fdr_10 INTEGER NULL
                    AS (ROUND(ln(1.0 / a) / b, 0)) STORED,
                fdr_20 INTEGER NULL
                    AS (ROUND(ln(2.0 / a) / b, 0)) STORED,
                squared_R FLOAT NULL,
                MAE FLOAT NULL,
                MAPE FLOAT NULL
            );
            """
        )
