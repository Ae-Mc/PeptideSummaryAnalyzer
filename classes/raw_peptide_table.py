"""См. класс RawPeptideTable."""

from decimal import Decimal
from typing import List

from .peptide_row import PeptideRow
from .base_classes import TableWithHeaders
from .peptide_columns import PeptideColumns


class RawPeptideTable(TableWithHeaders):
    """Считывает PeptideTable в список вида List[PeptideRow]."""

    columns: PeptideColumns

    def __init__(
        self,
        table_filename: str,
        unsafe_flag: bool = False,
        columns: PeptideColumns = PeptideColumns(),
    ):
        super().__init__(table_filename, unsafe_flag, columns)

    def load(self, table_filename) -> List[PeptideRow]:
        super().load(table_filename)
        new_elements: List[PeptideRow] = []
        for line in self:
            if line[self.columns.n[0]] == "":
                break
            new_elements.append(
                PeptideRow(
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
                        if len(line[self.columns.peptide_intensity[0]].strip())
                        > 0
                        else Decimal(0)
                    ),
                    sequence=line[self.columns.sequence[0]],
                )
            )
        self.clear()
        self.extend(new_elements)
        return self
