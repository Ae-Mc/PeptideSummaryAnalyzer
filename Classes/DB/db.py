from sqlite3.dbapi2 import Connection, Cursor, connect
from sys import stdout
from typing import Any, Dict, Iterable, List, TextIO, Tuple

from Classes.DB.Creators import Creators
from Classes.DB.FDR import FDR
from Classes.DB.Fillers import Fillers
from Classes.DB.Functions import Functions
from Classes.DB.ProteinGrouping import ProteinGrouping


class DB:
    """Отвечает за создание базы данных и функций.

    Содержит другие классы дл взаимодействия с БД.

    Attributes:
        connection: соединение с БД
        cursor: connection.cuesor()
        initializers: класс отвечающий за создание таблиц
        proteinGrouping: класс, отвечающий за работу protein grouping фильтра
        fdr: класс с FDR фильтрами"""

    connection: Connection
    cursor: Cursor
    initializers: Creators
    functions: Functions
    fillers: Fillers
    paramFilters: ProteinGrouping
    fdr: FDR

    def __init__(self) -> None:
        self.connection = connect(":memory:")
        self.cursor = self.connection.cursor()
        self.initializers = Creators(self.cursor)
        self.initializers.createAllFunctions()
        self.initializers.createAllTables()
        self.functions = Functions(self.cursor)
        self.fillers = Fillers(self.cursor)
        self.proteinGrouping = ProteinGrouping(self.cursor)
        self.fdr = FDR(self.cursor)

    def getDB(self) -> List[Tuple[Any]]:
        return self.execute(
            """SELECT row_id, table_number, accession, confidence
            FROM
                (SELECT id, table_number, confidence
                 FROM peptide_row) AS peptide_row
                JOIN peptide_accession ON row_id = peptide_row.id"""
        ).fetchall()

    def execute(self, sql: str, parameters: Iterable[Any] = []) -> Cursor:
        return self.cursor.execute(sql, parameters)

    def prettyFetchAll(self, sql: str, parameters: Iterable[Any] = []) -> list:
        cursor = self.cursor.execute(sql, parameters)
        result = cursor.fetchall()
        result.insert(0, tuple(map(lambda x: x[0], cursor.description)))
        return result

    def prettyPrintRepresentatives(self, file: TextIO = stdout) -> None:
        rows = self.execute(
            """SELECT repr.*, acc_g.accession, acc_g.count_in_table
            FROM representative repr
                INNER JOIN accession_group acc_g ON repr.id = acc_g.representative_id
            WHERE (
                SELECT COUNT(*)
                FROM accession_group acc_g
                WHERE acc_g.representative_id = repr.id
            ) > 1;
            """
        ).fetchall()
        tableNumbers = tuple(
            map(
                lambda x: x[0],
                self.execute(
                    """SELECT DISTINCT table_number
            FROM representative repr
            WHERE (
                SELECT COUNT(*)
                FROM accession_group acc_g
                WHERE acc_g.representative_id = repr.id
            ) > 1
            ORDER BY CAST(table_number AS DECIMAL);"""
                ).fetchall(),
            )
        )
        print("Repr\tAcc\t", "\t".join(tableNumbers), file=file)
        groups: Dict[str, Dict[str, List[str]]] = {}
        for row in rows:
            if row[3] not in groups:
                groups[row[3]] = {}
            if row[4] not in groups[row[3]]:
                groups[row[3]][row[4]] = ["0" for _ in tableNumbers]
            groups[row[3]][row[4]][tableNumbers.index(row[2])] = str(row[5])

        for representative, accessions in groups.items():
            print(
                representative,
                end="\t\t" + "\t".join(accessions[representative]) + "\n",
                file=file,
            )

            for accession, counts in filter(
                lambda p: p[0] != representative, accessions.items()
            ):
                print("\t".join(["", accession, *counts]), file=file)

    def close(self) -> None:
        self.connection.close()
