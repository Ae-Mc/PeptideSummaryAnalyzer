"""Отвечает за создание выходных файлов"""
from os import makedirs, path
from sqlite3 import Cursor
from typing import Any, Dict, List, Optional, TextIO, Tuple

from Classes.Input import Input


class Output:
    """Отвечает за вывод полученной информации из БД в файлы.

    Attributes:
        cursor (Cursor): экземляр класса Cursor, через который происходит связь
            с БД.
        input_params (Input): входные параметры скрипта.
        table_numbers (List[str]): номера таблиц, отсортированные по
            возрастанию.
    """

    cursor: Cursor
    input_params: Input
    table_numbers: List[str]

    def __init__(self, cursor: Cursor, inputParams: Input) -> None:
        self.cursor = cursor
        self.input_params = inputParams

    def generate_output_files(self) -> None:
        """Отвечает за создание и заполнение всех выходных файлов."""

        if not path.exists(self.input_params.outputPath):
            makedirs(self.input_params.outputPath)

        self._fill_table_numbers()

        # поле типа bool обозначает нужно ли добавлять столбцы Description и
        # Sequence Length
        float_columns_to_files: Tuple[Tuple[str, str, bool], ...] = (
            ("sc_norm_to_file_norm_ratio", "Sc_norm.txt", True),
            ("sc_sum", "Sc_summ.txt", True),
            (
                "peptide_intensity_norm_to_file_norm_ratio",
                "Pep_intensity_norm.txt",
                True,
            ),
            ("peptide_intensity_sum", "Pep_intensity_summ.txt", True),
        )
        for column, filename, additional_columns in float_columns_to_files:
            self.generate_table_file_by_float_column(
                filename, column, additional_columns
            )

        columns_to_files = (
            ("count", "Pep_counts.txt"),
            ("seq_length_sum", "Pep_seq_length_summ.txt"),
        )

        for column, filename in columns_to_files:
            self.generate_table_file_by_column(filename, column, False)

        self.generate_sequences_files(("Sequences.fasta", "Sequences.txt"))
        self.generate_settings_file("Settings.txt")
        self.generate_protein_groups_file("Protein groups.txt")
        self.generate_proteins_in_groups_file("Proteins in groups.txt")

    def _fill_table_numbers(self) -> None:
        """Заполняет аттрибут tableNumbers"""

        self.table_numbers = list(
            map(
                lambda x: x[0],
                self.cursor.execute(
                    """--sql
                    SELECT DISTINCT table_number FROM peptide_with_sum
                    ORDER BY CAST(table_number AS FLOAT);"""
                ).fetchall(),
            )
        )

    def generate_table_file_by_float_column(
        self, filename: str, column: str, include_additional_columns: bool
    ) -> None:
        """Создаёт выходной файл на основе определённого столбца таблицы
        peptide_with_sum, округляя его до 9 знаков после запятой.

        Args:
            filename (str): имя выходного файла.
            column (str): имя столбца.
            include_additional_columns (bool): включать ли стобцы
                Sequence Length и Description в выходной файл.
        """
        self._generate_table_file_by_column(
            filename,
            include_additional_columns,
            f"RTRIM(RTRIM(ROUND(peptide_with_sum.{column}, 9), '0'), '.')",
        )

    def generate_table_file_by_column(
        self, filename: str, column: str, include_additional_columns: bool
    ) -> None:
        """Создаёт выходной файл на основе определённого столбца таблицы
        peptide_with_sum.

        Args:
            filename (str): имя выходного файла.
            column (str): имя столбца.
            include_additional_columns (bool): включать ли стобцы
                Sequence Length и Description в выходной файл.
        """
        self._generate_table_file_by_column(
            filename, include_additional_columns, f"peptide_with_sum.{column}"
        )

    def _generate_table_file_by_column(
        self, filename: str, include_additional_columns: bool, column_sql: str
    ) -> None:
        with open(
            self.get_joint_output_filename(filename), "w", encoding="utf-8"
        ) as out_file:
            out_file.write(
                "\t".join(
                    [
                        "Accession",
                        *(
                            ["Description", "Sequence length"]
                            if include_additional_columns
                            else []
                        ),
                        *self.table_numbers,
                    ]
                )
            )

            accession: str
            table_number: str
            description: str
            sequence_length: int
            results: Dict[str, Tuple[str, int, List[Any]]] = {}
            for (
                accession,
                table_number,
                description,
                sequence_length,
                value,
            ) in self.cursor.execute(
                f"""--sql
                SELECT
                    accession,
                    table_number,
                    description,
                    LENGTH(sequence.sequence),
                    {column_sql}
                FROM peptide_with_sum JOIN sequence USING(accession)
                ORDER BY accession, CAST(table_number AS FLOAT);"""
            ).fetchall():
                if accession not in results:
                    results[accession] = (
                        description,
                        sequence_length,
                        ["0" for _ in self.table_numbers],
                    )

                results[accession][2][
                    self.table_numbers.index(table_number)
                ] = str(value)

            for accession, columns in results.items():
                out_file.write(f"\n{accession}\t")
                if include_additional_columns:
                    out_file.write(f"{columns[0]}\t{columns[1]}\t")
                out_file.write("\t".join(columns[2]))

    def generate_protein_groups_file(self, filename: str) -> None:
        """Создаёт файл со списком Protein групп.

        Args:
            filename (str): имя выходного файла.
        """
        rows = self.cursor.execute(
            """--sql
            SELECT repr.representative, accession, table_number, count
            FROM representative repr
                JOIN accession_group acc_g ON repr.id = acc_g.representative_id
                JOIN accession_count_per_table USING (accession)
            WHERE (
                SELECT COUNT(*)
                FROM accession_group acc_g
                WHERE acc_g.representative_id = repr.id
            ) > 1;
            """
        ).fetchall()
        with self.open_output_file(filename) as out_file:
            out_file.write(
                "Representative\tAccession\t" + "\t".join(self.table_numbers)
            )
            groups: Dict[str, Dict[str, List[str]]] = {}
            for row in rows:
                if row[0] not in groups:
                    groups[row[0]] = {}
                if row[1] not in groups[row[0]]:
                    groups[row[0]][row[1]] = ["0" for _ in self.table_numbers]
                groups[row[0]][row[1]][self.table_numbers.index(row[2])] = str(
                    row[3]
                )

            for representative, accessions in groups.items():
                out_file.write(
                    # pylint: disable=consider-using-f-string
                    "\n{}\t\t{}\n".format(
                        representative, "\t".join(accessions[representative])
                    )
                )

                accessions_without_representative = {
                    accession: counts
                    for accession, counts in accessions.items()
                    if accession != representative
                }

                for (
                    accession,
                    counts,
                ) in accessions_without_representative.items():
                    out_file.write("\t".join(["", accession, *counts]))

    def generate_proteins_in_groups_file(self, filename: str) -> None:
        """Создаёт файл с количеством accession, существующих в каждой группе.

        Имеет формат:
        Accession           ID in group
        [representative]    [count]

        Args:
            filename (str): имя выходного файла.
        """
        rows = self.cursor.execute(
            """--sql
            SELECT representative, COUNT(*) FROM (
                SELECT DISTINCT repr.representative, acc_g.accession
                FROM representative repr JOIN accession_group acc_g
                    ON repr.id = acc_g.representative_id
            )
            GROUP BY representative
            HAVING representative IN (
                SELECT DISTINCT accession FROM peptide_with_sum
            );"""
        ).fetchall()
        with self.open_output_file(filename) as out_file:
            out_file.write("Accession\tID in group")
            representative: str
            accession_count: int
            for representative, accession_count in rows:
                out_file.write(f"\n{representative}\t{accession_count}")

    def generate_sequences_files(self, filenames: Tuple[str, str]) -> None:
        """Генерирует списки последовательностей по найденным Accession

        Args:
            filenames: имена выходных файлов - первый элемент для fasta файла,
                второй - для txt файла, который можно будет импортировать в
                Excel
        """
        fasta_filename = filenames[0]
        txt_filename = filenames[1]
        with self.open_output_file(fasta_filename) as out_file:
            accession: str
            raw_sequence: str
            description: Optional[str]
            for accession, raw_sequence, description in self.cursor.execute(
                """--sql
                SELECT DISTINCT
                    sequence.accession,
                    sequence.raw_sequence,
                    description
                FROM sequence JOIN peptide_with_sum USING(accession)
                ORDER BY accession;"""
            ).fetchall():
                out_file.write(">" + accession)
                if description is not None:
                    out_file.write(" " + description)
                out_file.write(f"\n{raw_sequence}\n")

        with self.open_output_file(txt_filename) as out_file:
            out_file.write("Accession\tSequence\n")
            sql = """--sql
                SELECT DISTINCT
                    sequence.accession,
                    sequence.raw_sequence
                FROM sequence JOIN peptide_with_sum USING(accession)
                ORDER BY accession;"""
            rows = [
                "\t".join(row) for row in self.cursor.execute(sql).fetchall()
            ]
            out_file.write("\n".join(rows))

    def generate_settings_file(self, filename: str) -> None:
        """Создаёт файл с информацией о параметрах запуска скрипта.

        Args:
            filename: имя выходного файла."""

        with self.open_output_file(filename) as out_file:
            text = f""""ProteinPilot summary analyzer"
                #Protein filter
                Global FDR critical value (<% k or default): {
                    self.input_params.getFDRStr()
                }
                ID exclusion list: {(self.input_params.blackList or [""])[0]}
                Peptide confidence (value or default): {
                    self.input_params.proteinConfidence or ""
                }
                Protein grouping (conf): {
                    self.input_params.proteinGroupingConfidence
                }
                #Peptide filter
                Peptide confidence (value): {self.input_params.confPeptide}
                #Output filter
                Min groups with ID: {self.input_params.minGroupsWithAccession}
                Max missing values per group: {self.input_params.maxGroupLack}
                """
            out_file.write("\n".join([m.strip() for m in text.split("\n")]))

    def get_joint_output_filename(self, filename: str) -> str:
        """Возвращает путь к файлу с учётом input_params.output_path.

        Args:
            filename (str): имя файла.

        Returns:
            str: путь к файлу с учётом input_params.output_path.
        """
        return path.join(self.input_params.outputPath, filename)

    def open_output_file(self, filename: str) -> TextIO:
        """Открывает выходной файл на запись.

        Args:
            filename (str): имя выходного файла.

        Returns:
            TextIO: файловый поток.
        """
        return open(
            self.get_joint_output_filename(filename), "w", encoding="utf-8"
        )
