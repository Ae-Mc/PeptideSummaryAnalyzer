"""Отвечает за работу с базой данных SQLite."""

from sqlite3 import Connection, Cursor, connect
from typing import Any, Iterable

from Classes.db.creators import Creators
from Classes.db.fdr import FDR
from Classes.db.fillers import Fillers
from Classes.db.functions import Functions
from Classes.db.output import Output
from Classes.db.output_grouping import OutputGrouping
from Classes.db.protein_grouping import ProteinGrouping
from Classes.db.summaries import Summaries
from Classes.input import Input


class DB:
    """Отвечает за создание базы данных и функций.

    Содержит другие классы дл взаимодействия с БД.

    Attributes:
        connection: соединение с БД
        cursor: connection.cursor()
        inputParams: входные параметры скрипта
        initializers: класс отвечающий за создание таблиц
        proteinGrouping: класс, отвечающий за работу protein grouping фильтра
        fdr: класс с FDR фильтрами
        summaries: класс, отвечающий за заполнение последней таблицы -
            peptide_with_sum
    """

    connection: Connection
    cursor: Cursor
    input_params: Input
    initializers: Creators
    functions: Functions
    fillers: Fillers
    param_filters: ProteinGrouping
    fdr: FDR
    summaries: Summaries
    output: Output
    output_grouping: OutputGrouping

    def __init__(self, input_params: Input) -> None:
        self.input_params = input_params

        self.connection = connect(":memory:")
        self.connection.execute("PRAGMA foreign_keys = 1;")
        self.cursor = self.connection.cursor()

        self.initializers = Creators(self.cursor)
        self.initializers.create_all_functions()
        self.initializers.create_all_tables()

        self.functions = Functions(self.cursor)
        self.fillers = Fillers(self.cursor)
        self.protein_grouping = ProteinGrouping(self.cursor)
        self.fdr = FDR(self.cursor)
        self.summaries = Summaries(self.cursor)
        self.output_grouping = OutputGrouping(self.cursor)
        self.output = Output(self.cursor, input_params)

    def execute(self, sql: str, parameters: Iterable[Any] = None) -> Cursor:
        """Сокращение для self.cursor.execute."""
        return self.cursor.execute(sql, parameters or [])

    def pretty_fetch_all(
        self, sql: str, parameters: Iterable[Any] = None
    ) -> list:
        """Работает как сокращение от self.cursor.execute(...).fetchall(), но
        добавляет имена столбцов первым элементом в возвращаемом списке строк.
        """
        cursor = self.cursor.execute(sql, parameters or [])
        result = cursor.fetchall()
        result.insert(0, tuple(map(lambda x: x[0], cursor.description)))
        return result

    def __enter__(self):
        return self

    def __exit__(self, *exc_info):
        self.connection.close()
