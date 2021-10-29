"""Модуль, содержащий абстрактные классы."""

from .table import Table
from .column_names import ColumnNames
from .table_with_headers import TableWithHeaders


__all__ = ["Table", "ColumnNames", "TableWithHeaders"]
