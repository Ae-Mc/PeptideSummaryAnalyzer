"""См. класс Fillers."""

from sqlite3.dbapi2 import Cursor
from typing import List

from classes import (
    PeptideRow,
    RawPeptideTable,
    RawTables,
    Sequence,
    SequenceDatabase,
)


class Fillers:
    """Отвечает за заполнение таблиц БД исходными данными, считанными из
    входных файлов.

    Attributes:
        cursor: Cursor, через который заполняются таблицы"""

    cursor: Cursor

    def __init__(self, cursor: Cursor) -> None:
        self.cursor = cursor

    def fill_sequence(self, seq_database: SequenceDatabase):
        """Заполняет таблицу sequence данными, считанными из .fasta файла в
        SequenceDatabase.

        Args:
            seq_database (SequenceDatabase): база данных Sequence, считанная из
                .fasta файла.
        """
        sequence: Sequence
        for sequence in seq_database.values():
            self.cursor.execute(
                """--sql
                INSERT INTO sequence (
                    accession, description, sequence, raw_sequence
                ) VALUES ((?), (?), (?), (?));
                """,
                [
                    sequence.accession,
                    sequence.desc,
                    sequence.seq,
                    sequence.raw_seq,
                ],
            )

    def fill_raw_peptide(self, peptide_tables: RawTables[RawPeptideTable]):
        """Заполняет таблицы peptide_row и peptide_accession данными,
        считанными из DistinctPeptideSummary файлов.

        Args:
            peptide_tables (RawTables[RawPeptideTable]): считанные данные
        """
        for table_num in peptide_tables.get_sorted_table_nums():
            row: PeptideRow
            for row in peptide_tables[table_num]:
                self.cursor.execute(
                    """--sql
                    INSERT INTO peptide_row (
                        table_number,
                        N,
                        confidence,
                        score,
                        peptide_intensity,
                        sequence
                    ) VALUES ((?), (?), (?), (?), (?), (?));
                    """,
                    [
                        table_num,
                        row.N,
                        float(row.confidence),
                        float(row.score),
                        float(row.peptide_intensity),
                        row.sequence,
                    ],
                )
                row_id = self.cursor.lastrowid
                for accession in row.accessions:
                    self.cursor.execute(
                        """--sql
                        INSERT INTO peptide_accession (row_id, accession)
                        VALUES ((?), (?));""",
                        [row_id, accession],
                    )

    def fill_exclusion(self, exclusiuon_list: List[str]):
        """Заполняет таблицу exclusion списком accession, которые подлежат
        исключению.

        Args:
            exclusion_list (List[str]): список accession, подлежащих
                исключению.
        """
        for accession in exclusiuon_list:
            self.cursor.execute(
                """--sql
                INSERT INTO exclusion (accession) VALUES ((?));
                """,
                [accession],
            )
