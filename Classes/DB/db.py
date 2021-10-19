from Classes.DB.OutputGrouping import OutputGrouping
from Classes.DB.Summaries import Summaries
from sqlite3.dbapi2 import Connection, Cursor, connect
from sys import stdout
from typing import Any, Dict, Iterable, List, TextIO, Tuple

from Classes.DB.Creators import Creators
from Classes.DB.FDR import FDR
from Classes.DB.Fillers import Fillers
from Classes.DB.Functions import Functions
from Classes.DB.ProteinGrouping import ProteinGrouping
from Classes.DB.Output import Output
from Classes.Input import Input


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
        summaries: класс, отвечающий за заполнение последней таблицы - peptide_with_sum
    """

    connection: Connection
    cursor: Cursor
    inputParams: Input
    initializers: Creators
    functions: Functions
    fillers: Fillers
    paramFilters: ProteinGrouping
    fdr: FDR
    summaries: Summaries
    output: Output
    outputGrouping: OutputGrouping

    def __init__(self, inputParams: Input) -> None:
        self.inputParams = inputParams

        self.connection = connect(":memory:")
        self.cursor = self.connection.cursor()

        self.initializers = Creators(self.cursor)
        self.initializers.createAllFunctions()
        self.initializers.createAllTables()

        self.functions = Functions(self.cursor)
        self.fillers = Fillers(self.cursor)
        self.proteinGrouping = ProteinGrouping(self.cursor)
        self.fdr = FDR(self.cursor)
        self.summaries = Summaries(self.cursor)
        self.outputGrouping = OutputGrouping(self.cursor)
        self.output = Output(self.cursor, inputParams)

    def __enter__(self):
        return self

    def getDB(self) -> List[tuple]:
        return self.prettyFetchAll(
            """--sql
            SELECT row_id, table_number, accession, confidence
            FROM
                (SELECT id, table_number, confidence
                 FROM peptide_row) AS peptide_row
                JOIN peptide_accession ON row_id = peptide_row.id;"""
        )

    def execute(self, sql: str, parameters: Iterable[Any] = []) -> Cursor:
        return self.cursor.execute(sql, parameters)

    def prettyFetchAll(self, sql: str, parameters: Iterable[Any] = []) -> list:
        cursor = self.cursor.execute(sql, parameters)
        result = cursor.fetchall()
        result.insert(0, tuple(map(lambda x: x[0], cursor.description)))
        return result

    def close(self) -> None:
        self.connection.close()

    def __exit__(self, *exc_info):
        self.close()
