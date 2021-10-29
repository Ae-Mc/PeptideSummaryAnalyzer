"""См. класс RawTables."""

from os import listdir
from re import match
from typing import Generic, List, Type, TypeVar
from classes.base_classes.column_names import ColumnNames

from .base_classes import TableWithHeaders


Table = TypeVar("Table", bound=TableWithHeaders)


class RawTables(dict, Generic[Table]):
    """Словарь вида {
        "номер таблицы": TableType,
    }

    Attributes:
        column_names: имена стобцов (заголовки стобцов)
    """

    column_names: ColumnNames
    filename_pattern: str
    table_to_instantiate: Type[Table]

    def __init__(
        self,
        input_dir: str,
        filename_pattern: str,
        column_names: ColumnNames,
        table_to_instantiate: Type[Table],
    ) -> None:
        """
        Args:
            columnNames: имена заголовков
            inputDir: путь, из которого считываются таблицы
        """

        self.column_names = column_names
        self.filename_pattern = filename_pattern
        self.table_to_instantiate = table_to_instantiate

        super().__init__()
        self.read_tables(input_dir)

    def read_tables(self, input_dir: str) -> None:
        """Считывает все файлы в словарь

        Args:
            inputDir: путь, из которого считываются таблицы
        """
        for filename in listdir(input_dir):
            if match(self.filename_pattern, filename):
                table_num = filename.split("_")[0]
                self[table_num] = self.table_to_instantiate(
                    input_dir + "/" + filename,
                    unsafe_flag=True,
                    columns=self.column_names,
                )

    def get_sorted_table_nums(self) -> List[str]:
        """Получает отсортированный список номеров таблиц
        Returns:
            отсортированный список номеров таблиц
        """
        return sorted(self.keys(), key=float)

    # pylint: disable=useless-super-delegation
    def __getitem__(self, k: str) -> Table:
        return super().__getitem__(k)

    def __setitem__(self, k: str, v: Table) -> None:
        return super().__setitem__(k, v)
