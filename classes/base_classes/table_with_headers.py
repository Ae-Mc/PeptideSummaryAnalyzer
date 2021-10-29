"""См. класс TableWithHeaders."""

from typing import List

from .table import Table
from .column_names import ColumnNames


class TableWithHeaders(Table):
    """Считывает таблицу в список строк, проверяет заголовки и образает их."""

    columns: ColumnNames

    def __init__(
        self,
        tableFilename: str,
        unsafe_flag: bool = False,
        columns: ColumnNames = ColumnNames(),
    ):
        """Initializes table settings and reads table from file"""
        self.columns = columns
        super().__init__(tableFilename, unsafe_flag)

    def load(self, table_filename) -> List:
        super().load(table_filename)
        self.columns.test_column_names(self.pop(0))
        return self
