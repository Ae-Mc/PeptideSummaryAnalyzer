"""Отвечает за заполнение таблиц БД исходными данными, считанными из входных
файлов."""

from sqlite3.dbapi2 import Cursor
from typing import List

from Classes.PeptideRow import PeptideRow
from Classes.RawPeptideTables import RawPeptideTables
from Classes.Sequence import Sequence
from Classes.SequenceDatabase import SequenceDatabase


class Fillers:
    """Отвечает за заполнение таблиц

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
                r"""INSERT INTO sequence (
              accession, description, sequence, raw_sequence
          ) VALUES (
                  (?), (?), (?), (?)
          );""",
                [
                    sequence.accession,
                    sequence.desc,
                    sequence.seq,
                    sequence.rawSeq,
                ],
            )

    def fill_raw_peptide(self, peptide_tables: RawPeptideTables):
        """Заполняет таблицы peptide_row и peptide_accession данными,
        считанными из DistinctPeptideSummary файлов.

        Args:
            peptide_tables (RawPeptideTables): считанные данные
        """
        for table_num in peptide_tables.GetSortedTableNums():
            row: PeptideRow
            for row in peptide_tables[table_num]:
                self.cursor.execute(
                    """--sql
            INSERT INTO peptide_row (
                table_number, confidence, score, peptide_intensity, sequence
            ) VALUES (
                (?), (?), (?), (?), (?)
            );""",
                    [
                        table_num,
                        float(row.confidence),
                        float(row.sc),
                        float(row.precursorSignal),
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
                """INSERT INTO exclusion (accession) VALUES (?)""", [accession]
            )