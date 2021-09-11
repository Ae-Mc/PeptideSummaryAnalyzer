from Classes.DB.db import DB
from Classes.PeptideRow import PeptideRow
from sqlite3 import Connection, Cursor
from typing import Any, Iterable, List, Tuple
from Classes.PeptideColumns import PeptideColumns
from Classes.RawPeptideTables import RawPeptideTables
from Classes.Input import Input
from pprint import pprint
from Classes.Functions import FindFastaFile, GetFileLines
from Classes.SequenceDatabase import SequenceDatabase
from Classes.Sequence import Sequence


class tempDB:

    connection: Connection
    cursor: Cursor

    def getDefaultConfidenceFiltredAccessions(self) -> List[Tuple[str, str]]:
        return self.execute(
            """SELECT DISTINCT table_number, accession FROM (
            SELECT table_number, accession
            FROM peptide_row JOIN peptide_accession
                ON peptide_row.id = row_id
            WHERE confidence >= 95 AND confidence < 99
            GROUP BY table_number, accession
            HAVING COUNT(*) >= 2
            UNION SELECT table_number, accession
            FROM peptide_row JOIN peptide_accession
                ON peptide_row.id = row_id
            WHERE confidence >= 99
            GROUP BY table_number, accession);"""
        ).fetchall()

    # TODO
    def addPrecursorAndScoreAccessionSums(self) -> None:
        pass

    # TODO
    def getTableSums(self) -> None:
        self.execute(
            """SELECT
                SUM(precursor_signal) AS precursor_signal_sum,
                SUM(score) AS score_sum"""
        )

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

    def close(self) -> None:
        self.connection.close()


if __name__ == "__main__":
    inp = Input()
    inp.rootPath = "./Presets/preset_IDexcl"
    inp.inputPath = inp.rootPath + "/Input"
    peptideTables = RawPeptideTables(PeptideColumns(), inp.inputPath)
    db = DB()
    inp.seqDB = SequenceDatabase.fromFile(FindFastaFile(inp.rootPath))
    # Заполняем таблицу для хранения информации о последовательностях
    sequence: Sequence
    for sequence in inp.seqDB.values():
        db.execute(
            r"""INSERT INTO sequence (
                accession, description, sequence
            ) VALUES (
                    (?), (?), (?)
            );""",
            [sequence.accession, sequence.desc, sequence.seq],
        )
    for tableNum in peptideTables.GetSortedTableNums():
        row: PeptideRow
        for row in peptideTables[tableNum]:
            db.execute(
                """INSERT INTO peptide_row (
                    table_number, confidence, score, precursor_signal, sequence
                ) VALUES (
                    (?), (?), (?), (?), (?)
                );""",
                [
                    tableNum,
                    float(row.confidence),
                    float(row.sc),
                    float(row.precursorSignal),
                    row.sequence,
                ],
            )
            rowID = db.cursor.lastrowid
            for accession in row.accessions:
                db.execute(
                    """INSERT INTO peptide_accession (
                        row_id, accession
                    ) VALUES (
                        (?), (?)
                    );""",
                    [rowID, accession],
                )
    # Получение всех accession, не присутствующих в .fasta бд
    absenceAccessions = db.execute(
        """SELECT DISTINCT accession FROM peptide_accession
        WHERE (
            NOT IS_REVERSED(accession)
            AND accession NOT IN (SELECT accession FROM sequence)
        );"""
    ).fetchall()
    if len(absenceAccessions) > 0:
        raise IndexError(
            "Accessions not found in .fasta file:\n\t"
            + "\n\t".join(absenceAccessions)
        )

    try:
        for accession in GetFileLines(inp.rootPath + "/IDexcl.txt") or []:
            db.execute(
                f"""INSERT INTO exclusion (accession) VALUES ("{accession}")"""
            )
    except FileNotFoundError:
        print("ID exclusion file not found")

    db.fdr.default()

    # Применение ID exclusion фильтра
    if (
        db.execute("SELECT COUNT(*) FROM exclusion LIMIT(1)").fetchall()[0][0]
        > 0
    ):
        db.execute(
            """DELETE FROM peptide_accession
            WHERE accession IN (SELECT accession FROM exclusion);"""
        )
        db.execute(
            """DELETE FROM peptide_row
            WHERE id IN (SELECT peptide_row.id AS id
                FROM peptide_row LEFT JOIN peptide_accession
                    ON row_id = peptide_row.id
                WHERE row_id IS NULL);"""
        )

    print("Количество появлений каждого accession по всем таблицам")
    pprint(
        db.execute(
            """CREATE TABLE overall_accession_count AS
            SELECT accession, COUNT(*) AS count
            FROM peptide_row JOIN peptide_accession
                ON peptide_row.id = peptide_accession.row_id
            GROUP BY accession
            ORDER BY accession"""
        ).fetchall()
    )
    print("Количество появлений каждого accession в каждой таблице")
    pprint(
        db.execute(
            """CREATE TABLE accession_count_per_table AS
            SELECT table_number, accession, COUNT(*) AS count
            FROM peptide_row JOIN peptide_accession
                ON peptide_row.id = peptide_accession.row_id
            GROUP BY table_number, accession
            ORDER BY table_number, accession"""
        ).fetchall()
    )

    print(
        "Выбор репрезентативных accession для строк и заполнение таблицы"
        " accession_table"
    )
    # TODO
    # pprint(
    #     db.execute(
    #         """SELECT row_id, accession FROM peptide_accession
    #         WHERE """))

    print("Получение accession, не соответствующих confidence default фильтру")
    pprint(db.paramFilters.confidenceDefault())

    pprint(db.getDB())
    db.close()
