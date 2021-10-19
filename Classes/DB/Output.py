from pprint import PrettyPrinter
from os import makedirs, path
from sqlite3 import Cursor
from typing import Any, Dict, List, Optional, Tuple

from Classes.Input import Input


class Output:
    """Отвечает за вывод полученной информации из БД в файлы.

    Attributes:
        cursor (Cursor): экземляр класса Cursor, через который происходит связь с БД
        inputParams (Input): входные параметры скрипта
        tableNumbers (List[str]): номера таблиц, отсортированные по возрастанию"""

    cursor: Cursor
    inputParams: Input
    tableNumbers: List[str]

    def __init__(self, cursor: Cursor, inputParams: Input) -> None:
        self.cursor = cursor
        self.inputParams = inputParams

    def GenerateOutputFiles(self) -> None:
        """Отвечает за создание и заполнение всех выходных файлов."""

        if not path.exists(self.inputParams.outputPath):
            makedirs(self.inputParams.outputPath)

        self._fillTableNumbers()

        # поле типа bool обозначает нужно ли добавлять столбцы Description и
        # Sequence Length
        floatColumnsToFiles: Tuple[Tuple[str, str, bool], ...] = (
            ("sc_norm_to_file_norm_ratio", "Sc_norm.txt", True),
            ("sc_sum", "Sc_summ.txt", True),
            (
                "peptide_intensity_norm_to_file_norm_ratio",
                "Pep_intensity_norm.txt",
                True,
            ),
            ("peptide_intensity_sum", "Pep_intensity_summ.txt", True),
        )
        for column, filename, additionalColumns in floatColumnsToFiles:
            self.GenerateTableFileByFloatColumn(filename, column, additionalColumns)

        columnsToFiles = (
            ("count", "Pep_counts.txt"),
            ("seq_length_sum", "Pep_seq_length_summ.txt"),
        )

        for column, filename in columnsToFiles:
            self.GenerateTableFileByColumn(filename, column, False)

        self.GenerateSequencesFiles(("Sequences.fasta", "Sequences.txt"))
        self.GenerateSettingsFile("Settings.txt")

    def _fillTableNumbers(self) -> None:
        """Заполняет аттрибут tableNumbers"""

        self.tableNumbers = list(
            map(
                lambda x: x[0],
                self.cursor.execute(
                    """--sql
                    SELECT DISTINCT table_number FROM peptide_with_sum
                    ORDER BY CAST(table_number AS FLOAT);"""
                ).fetchall(),
            )
        )

    def GenerateTableFileByFloatColumn(
        self, filename: str, column: str, includeAdditionalColumns: bool
    ) -> None:
        self._GenerateTableFileByColumn(
            filename,
            includeAdditionalColumns,
            f"RTRIM(RTRIM(ROUND(peptide_with_sum.{column}, 9), '0'), '.')",
        )

    def GenerateTableFileByColumn(
        self, filename: str, column: str, includeAdditionalColumns: bool
    ) -> None:
        self._GenerateTableFileByColumn(
            filename, includeAdditionalColumns, f"peptide_with_sum.{column}"
        )

    def _GenerateTableFileByColumn(
        self, filename: str, includeAdditionalColumns: bool, columnSql: str
    ) -> None:
        with open(self.GetJointOutputFilename(filename), "w") as outFile:
            outFile.write(
                "\t".join(
                    [
                        "Accession",
                        *(
                            ["Description", "Sequence length"]
                            if includeAdditionalColumns
                            else []
                        ),
                        *self.tableNumbers,
                    ]
                )
            )

            previousAccession: Optional[str] = None
            accession: str
            tableNumber: str
            results: Dict[str, Tuple[str, int, List[Any]]] = {}
            for (
                accession,
                tableNumber,
                description,
                sequenceLength,
                value,
            ) in self.cursor.execute(
                f"""--sql
                SELECT
                    accession,
                    table_number,
                    description,
                    LENGTH(sequence.sequence),
                    {columnSql}
                FROM peptide_with_sum JOIN sequence USING(accession)
                ORDER BY accession, CAST(table_number AS FLOAT);"""
            ).fetchall():
                if accession not in results:
                    results[accession] = (
                        description,
                        sequenceLength,
                        ["0" for _ in self.tableNumbers],
                    )

                results[accession][2][self.tableNumbers.index(tableNumber)] = str(value)

            for accession, columns in results.items():
                outFile.write(f"\n{accession}\t")
                    if includeAdditionalColumns:
                    outFile.write(f"{columns[0]}\t{columns[1]}\t")
                outFile.write("\t".join(columns[2]))

    def GenerateSequencesFiles(self, filenames: Tuple[str, str]) -> None:
        """Генерирует списки последовательностей по найденным Accession

        Args:
            filenames: имена выходных файлов - первый элемент для fasta файла,
                второй - для txt файла, который можно будет импортировать в
                Excel
        """
        fastaFilename = filenames[0]
        txtFilename = filenames[1]
        with open(self.GetJointOutputFilename(fastaFilename), "w") as outFile:
            accession: str
            rawSequence: str
            description: Optional[str]
            for accession, rawSequence, description in self.cursor.execute(
                """--sql
                SELECT DISTINCT sequence.accession, sequence.raw_sequence, description
                FROM sequence JOIN peptide_with_sum USING(accession)
                ORDER BY accession;"""
            ).fetchall():
                outFile.write(">" + accession)
                if description is not None:
                    outFile.write(" " + description)
                outFile.write(f"\n{rawSequence}\n")

        with open(self.GetJointOutputFilename(txtFilename), "w") as outFile:
            outFile.write("Accession\tSequence\n")
            sql = """--sql
                SELECT DISTINCT
                    sequence.accession,
                    sequence.raw_sequence
                FROM sequence JOIN peptide_with_sum USING(accession)
                ORDER BY accession;"""
            rows = ["\t".join(row) for row in self.cursor.execute(sql).fetchall()]
            outFile.write("\n".join(rows))

    def GenerateSettingsFile(self, filename: str) -> None:
        """Создаёт файл с информацией о параметрах запуска скрипта.

        Args:
            filename: имя выходного файла."""

        with open(self.GetJointOutputFilename(filename), "w") as outFile:
            text = f""""ProteinPilot summary analyzer"
                #Protein filter
                Global FDR critical value (<% k or default): {
                    self.inputParams.getFDRStr()
                }
                ID exclusion list: {(self.inputParams.blackList or [""])[0]}
                Peptide confidence (value or default): {
                    self.inputParams.proteinConfidence or ""
                }
                Protein grouping (conf): {
                    self.inputParams.proteinGroupingConfidence
                }
                #Peptide filter
                Peptide confidence (value): {self.inputParams.confPeptide}
                #Output filter
                Min groups with ID: {self.inputParams.minGroupsWithAccession}
                Max missing values per group: {self.inputParams.maxGroupLack}"""
            outFile.write("\n".join([m.strip() for m in text.split("\n")]))

    def GetJointOutputFilename(self, filename: str):
        return path.join(self.inputParams.outputPath, filename)
