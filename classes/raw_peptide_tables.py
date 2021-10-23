"""См. класс RawPeptideTables."""

from os import listdir
from typing import List

from .raw_peptide_table import RawPeptideTable
from .peptide_columns import PeptideColumns


class RawPeptideTables(dict):
    """Словарь вида {
        "номер таблицы": RawPeptideTable,
    }

    Attributes:
        column_names: имена стобцов (заголовки стобцов)
    """

    column_names: PeptideColumns

    def __init__(
        self, column_names: PeptideColumns, input_dir: str = None
    ) -> None:
        """
        Args:
            columnNames: имена заголовков
            inputDir: путь, из которого считываются таблицы
        """

        self.column_names = column_names

        super().__init__()
        if input_dir is not None:
            self.read_peptide_summaries(input_dir)

    def read_peptide_summaries(self, input_dir: str) -> None:
        """Считывает все PeptideSummary файлы в словарь

        Args:
            inputDir: путь, из которого считываются таблицы
        """
        for filename in listdir(input_dir):
            if "Peptide" in filename:
                table_num = filename.split("_")[0]
                self[table_num] = RawPeptideTable(
                    input_dir + "/" + filename,
                    unsafeFlag=True,
                    columns=self.column_names,
                )

    def get_sorted_table_nums(self) -> List[str]:
        """Получает отсортированный список номеров таблиц
        Returns:
            отсортированный список номеров таблиц
        """
        return sorted(self.keys(), key=float)

    # pylint: disable=useless-super-delegation
    def __getitem__(self, k: str) -> RawPeptideTable:
        return super().__getitem__(k)

    def __setitem__(self, k: str, v: RawPeptideTable) -> None:
        return super().__setitem__(k, v)
