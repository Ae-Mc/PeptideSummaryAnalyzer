"""См. класс RawProteinTable."""

from decimal import Decimal
from typing import List

from .protein_row import ProteinRow
from .base_classes import TableWithHeaders
from .protein_columns import ProteinColumns


class RawProteinTable(TableWithHeaders):
    """Считывает ProteinTable в список вида List[ProteinRow]."""

    columns: ProteinColumns

    def load(self, table_filename) -> List[ProteinRow]:
        super().load(table_filename)
        latest_unused = Decimal(0)
        for i, line in enumerate(self):
            unused = Decimal(line[self.columns.unused[0]])
            if unused != Decimal(0):
                latest_unused = unused

            self[i] = ProteinRow(
                n=int(line[self.columns.n[0]]),
                accession=line[self.columns.accession[0]],
                unused=latest_unused,
            )
        return self
