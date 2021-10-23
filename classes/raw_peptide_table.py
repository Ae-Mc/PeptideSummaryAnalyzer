"""См. класс RawPeptideTable."""

from decimal import Decimal
from typing import List

from .peptide_row import PeptideRow
from .base_classes import Table
from .peptide_columns import PeptideColumns


class RawPeptideTable(Table):
    """Считывает PeptideTable в список вида List[PeptideRow]."""

    columns: PeptideColumns

    def __init__(
        self,
        tableFilename: str = None,
        unsafeFlag: bool = False,
        columns: PeptideColumns = PeptideColumns(),
    ):
        """Initializes table settings and reads table from file"""
        self.columns = columns
        super().__init__(tableFilename, unsafeFlag)

    def load(self, table_filename) -> List[PeptideRow]:
        super().load(table_filename)
        self.columns.test_column_names(self.pop(0))
        for i, line in enumerate(self):
            self[i] = PeptideRow(
                accessions=list(
                    map(
                        lambda x: x.strip(),
                        line[self.columns.accession[0]].split(";"),
                    )
                ),
                confidence=Decimal(line[self.columns.confidence[0]]),
                score=Decimal(line[self.columns.score[0]]),
                peptide_intensity=(
                    Decimal(line[self.columns.peptide_intensity[0]])
                    if len(line[self.columns.peptide_intensity[0]].strip()) > 0
                    else Decimal(0)
                ),
                sequence=line[self.columns.sequence[0]],
            )
        return self

    def __str__(self):
        return self.__repr__()

    def __repr__(self):
        return "[\n " + "\n ".join([str(e) for e in self]) + "\n]"
