from os import makedirs, path
from sqlite3 import Cursor
from typing import Optional, Tuple

from Classes.Input import Input


class Output:
    """Отвечает за вывод полученной информации из БД в файлы."""

    cursor: Cursor
    inputParams: Input

    def __init__(self, cursor: Cursor, inputParams: Input) -> None:
        self.cursor = cursor
        self.inputParams = inputParams

    def GenerateOutputFiles(self) -> None:
        """Отвечает за создание и заполнение всех выходных файлов."""

        if not path.exists(self.inputParams.outputPath):
            makedirs(self.inputParams.outputPath)

        columnsToFiles: Tuple[Tuple[str, str], ...] = (
            ("count", "Pep_counts.txt"),
            ("seq_length_sum", "Pep_seq_length_summ.txt"),
            ("sc_norm_to_file_norm_ratio", "Sc_norm.txt"),
            ("sc_sum", "Sc_summ.txt"),
            ("peptide_intensity_norm_to_file_norm_ratio", "Pep_intensity_norm.txt"),
            ("peptide_intensity_sum", "Pep_intensity_summ.txt"),
        )

        for column, filename in columnsToFiles:
            self.GenerateTableFileByColumn(filename, column)

        if self.inputParams.shouldExtractSequences:
            self.GenerateSequencesFiles(("Sequences.fasta", "Sequences.txt"))
        self.GenerateDescriptionFile("Description.txt")
        self.GenerateSettingsFile("Settings.txt")

    def GenerateTableFileByColumn(self, filename: str, column: str) -> None:
        with open(self.GetJointOutputFilename(filename), "w") as outFile:
            tableNumbers = list(
                map(
                    lambda x: x[0],
                    self.cursor.execute(
                        """--sql
                        SELECT DISTINCT table_number FROM peptide_with_sum
                        ORDER BY CAST(table_number AS FLOAT);"""
                    ).fetchall(),
                )
            )
            outFile.write("\t".join(["Accession", *tableNumbers]))

            previousAccession: Optional[str] = None
            accession: str
            for (accession, value) in self.cursor.execute(
                f"""--sql
                SELECT accession, {column} FROM peptide_with_sum
                ORDER BY accession, CAST(table_number AS FLOAT);"""
            ).fetchall():
                if accession != previousAccession:
                    outFile.write(f"\n{accession}")
                    previousAccession = accession
                outFile.write(f"\t{value}")

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
            row: Tuple[str, str, Optional[str]]
            for row in self.cursor.execute(
                """--sql
                SELECT DISTINCT sequence.accession, sequence.sequence, description
                FROM sequence JOIN peptide_with_sum USING(accession)
                ORDER BY accession;"""
            ).fetchall():
                outFile.write(">" + row[0])
                if row[2] is not None:
                    outFile.write(" " + row[2])
                outFile.write("\n")
                # Разбиваем строку на строки по 60 символов, как в оригинальном .fasta
                # файле.
                for i in range(0, len(row[1]), 60):
                    outFile.write(row[1][i : min(i + 60, len(row[1]))])
                    outFile.write("\n")

        with open(self.GetJointOutputFilename(txtFilename), "w") as outFile:
            outFile.write("ID\tSequence\n")
            outFile.write(
                "\n".join(
                    [
                        "\t".join(row)
                        for row in self.cursor.execute(
                            """--sql
                            SELECT DISTINCT
                                sequence.accession,
                                sequence.sequence,
                                IFNULL(description, "")
                            FROM sequence JOIN peptide_with_sum USING(accession)
                            ORDER BY accession;"""
                        ).fetchall()
                    ]
                )
            )

    def GenerateDescriptionFile(self, filename: str) -> None:
        """Создаёт файл с Accession и Description из .fasta БД

        Args:
            filename (str): имя выходного файла.
        """

        with open(self.GetJointOutputFilename(filename), "w") as outFile:
            outFile.write("Accession\tDescription\n")
            outFile.write(
                "\n".join(
                    [
                        "\t".join(row)
                        for row in self.cursor.execute(
                            """--sql
                            SELECT DISTINCT sequence.accession, IFNULL(description, "")
                            FROM sequence JOIN peptide_with_sum USING(accession)
                            ORDER BY accession;"""
                        ).fetchall()
                    ]
                )
            )

    def GenerateSettingsFile(self, filename: str) -> None:
        """Создаёт файл с информацией о параметрах запуска скрипта.

        Args:
            filename: имя выходного файла."""

        with open(self.GetJointOutputFilename(filename), "w") as outFile:
            text = f""""ProteinPilot summary analyzer"
                #Protein filter
                Global FDR critical value (<% k or default): {self.inputParams.fdr}
                ID exclusion list: {(self.inputParams.blackList or [""])[0]}
                Peptide confidence (value or default): {
                    self.inputParams.isProteinConfidence or ""
                }
                Protein grouping (conf): {
                    self.inputParams.proteinGroupingConfidence
                }
                #Peptide filter
                Peptide confidence (value): {self.inputParams.confPeptide}
                #Output filter
                Min groups with ID: {self.inputParams.minGroupsWithAccession}
                Max missing values per group: {self.inputParams.maxGroupLack}
                ID Extract sequences? (Y/N): {
                    "Y" if self.inputParams.shouldExtractSequences else "N"}"""
            outFile.write("\n".join([m.strip() for m in text.split("\n")]))

    def GetJointOutputFilename(self, filename: str):
        return path.join(self.inputParams.outputPath, filename)
