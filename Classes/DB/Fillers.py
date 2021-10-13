from sqlite3.dbapi2 import Cursor
from typing import List

from Classes.PeptideRow import PeptideRow
from Classes.RawPeptideTables import RawPeptideTables
from Classes.Sequence import Sequence
from Classes.SequenceDatabase import SequenceDatabase


class Fillers:
    """Отвечает за заполнение таблиц

    Attributes:
        cursor: cursor, через который заполняются таблицы"""

    cursor: Cursor

    def __init__(self, cursor: Cursor) -> None:
        self.cursor = cursor

    def fillSequence(self, seqDB: SequenceDatabase):
        sequence: Sequence
        for sequence in seqDB.values():
            self.cursor.execute(
                r"""INSERT INTO sequence (
              accession, description, sequence, raw_sequence
          ) VALUES (
                  (?), (?), (?), (?)
          );""",
                [sequence.accession, sequence.desc, sequence.seq, sequence.rawSeq],
            )

    def fillRawPeptide(self, peptideTables: RawPeptideTables):
        for tableNum in peptideTables.GetSortedTableNums():
            row: PeptideRow
            for row in peptideTables[tableNum]:
                self.cursor.execute(
                    r"""INSERT INTO peptide_row (
                        table_number, confidence, score, peptide_intensity, sequence
                    ) VALUES (
                        (?), (?), (?), (?), (?)
                    );""",
                    [
                        tableNum,
                        float(row.confidence),
                        float(row.sc),
                        float(row.precursorSignal),
                        row.sequence,
                    ],
                )
                rowID = self.cursor.lastrowid
                for accession in row.accessions:
                    self.cursor.execute(
                        r"""INSERT INTO peptide_accession (
                            row_id, accession
                        ) VALUES (
                            (?), (?)
                        );""",
                        [rowID, accession],
                    )

    def fillExclusion(self, exclusiuonList: List[str]):
        for accession in exclusiuonList:
            self.cursor.execute(
                rf"""INSERT INTO exclusion (accession) VALUES ("{accession}")"""
            )
