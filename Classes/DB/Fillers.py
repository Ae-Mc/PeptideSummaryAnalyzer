from typing import List
from Classes.PeptideRow import PeptideRow
from Classes.RawPeptideTables import RawPeptideTables
from Classes.DB.db import DB
from Classes.Sequence import Sequence
from Classes.SequenceDatabase import SequenceDatabase


class Fillers:
    db: DB

    def __init__(self, db: DB) -> None:
        self.db = db

    def fillSequence(self, seqDB: SequenceDatabase):
        sequence: Sequence
        for sequence in seqDB.values():
            self.db.execute(
                r"""INSERT INTO sequence (
              accession, description, sequence
          ) VALUES (
                  (?), (?), (?)
          );""",
                [sequence.accession, sequence.desc, sequence.seq],
            )

    def fillPeptide(self, peptideTables: RawPeptideTables):
        for tableNum in peptideTables.GetSortedTableNums():
            row: PeptideRow
            for row in peptideTables[tableNum]:
                self.db.execute(
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
                rowID = self.db.cursor.lastrowid
                for accession in row.accessions:
                    self.db.execute(
                        r"""INSERT INTO peptide_accession (
                            row_id, accession
                        ) VALUES (
                            (?), (?)
                        );""",
                        [rowID, accession],
                    )

    def fillExclusion(self, exclusiuonList: List[str]):
        for accession in exclusiuonList:
            self.db.execute(
                rf"""INSERT INTO exclusion (accession) VALUES ("{accession}")"""
            )
