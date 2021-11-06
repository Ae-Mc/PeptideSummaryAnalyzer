"""См. класс RawPeptideTable."""

from decimal import Decimal
from typing import List

from .peptide_row import PeptideRow
from .base_classes import TableWithHeaders
from .peptide_columns import PeptideColumns


class RawPeptideTable(TableWithHeaders):
    """Считывает PeptideTable в список вида List[PeptideRow]."""

    columns: PeptideColumns

    def load(self, table_filename) -> List[PeptideRow]:
        super().load(table_filename)
        filtered_elements: List[PeptideRow] = []
        for line in self:
            if line[self.columns.N[0]].strip() == "":
                break
            filtered_elements.append(
                PeptideRow(
                    N=int(line[self.columns.N[0]]),
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
        self.extend(filtered_elements)
        return self
