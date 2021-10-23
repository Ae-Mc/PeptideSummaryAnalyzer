"""Модуль для тестирования программы на подготовленных пресетах."""

from traceback import print_exc
from os import listdir, path, remove
from sys import argv, exit as sys_exit
from typing import List
from Classes.input import Input
from Classes.SequenceDatabase import SequenceDatabase
from Classes.Functions import GetFileLines, FindFastaFile
from sql import main as protein_main

PRESETS_FOLDER = "Presets"
INPUTPATH = ""


class Preset:
    """Класс для проверки пресетов и сравнения их с результатами тестов

    Attributes:
        folder: папка с пресетом
        presetOutputDir: папка, в которой хранятся "правильные" выходные файлы
        settings: настройки пресета
        errorCode: код ошибки, если она была
    """

    folder: str
    preset_output_dir: str
    settings: Input
    error_code: int = 0

    def __init__(self, folder: str):
        """
        Args:
            folder: папка с пресетом
        """
        self.folder = folder
        self.preset_output_dir = path.join(self.folder, "TrueOutput")

    def run(self) -> None:
        """Запуск теста пресета"""
        if self.error_code == 0:
            try:
                print(self.folder)
                protein_main(self.settings)
            except Exception:  # pylint: disable=broad-except
                print(f"Error running preset {self.folder}")
                print_exc()
            self.test_result()

    def read_settings(self) -> None:
        """Считывание настроек пресета из папки с пресетом"""
        self.settings: Input = Input()
        self.settings.inputPath = path.join(self.folder, "Input")
        self.settings.output_path = path.join(self.folder, "Output")
        try:
            self.test_preset_file_existance()
        except FileNotFoundError as error:
            print(error.args[0])
            return

        with open(
            path.join(self.folder, "preset.txt"), encoding="utf-8"
        ) as preset_file:
            preset_file_lines = preset_file.read().split("\n")
            try:
                self.test_preset_file(preset_file_lines)
            except IOError as error:
                print(
                    f"Error reading preset {path.split(self.folder)[1]}'s"
                    f" {error.args[0]}"
                )
                self.error_code = 1
                return
            preset_file_values = [
                line.split(":")[1].strip()
                for line in preset_file_lines
                if len(line.split(":")) > 1
            ]

        self.settings.set_fdr(preset_file_values[0])
        self.settings.exclusion_list = None
        black_list = GetFileLines(
            path.join(self.folder, preset_file_values[1])
        )
        if len(preset_file_values[1].strip()):
            self.settings.exclusion_list = (
                None
                if black_list is None
                else (preset_file_values[1], black_list)
            )
        self.settings.seq_db = SequenceDatabase.fromFile(
            FindFastaFile(self.folder)
        )
        self.settings.set_protein_confidence(preset_file_values[2])
        self.settings.set_protein_grouping_confidence(preset_file_values[3])
        self.settings.set_conf_peptide(preset_file_values[4])
        if (
            preset_file_values[5] is not None
            and len(preset_file_values[5].strip()) > 0
        ):
            self.settings.min_groups_with_accession = int(
                preset_file_values[5]
            )
        if (
            preset_file_values[6] is not None
            and len(preset_file_values[6].strip()) > 0
        ):
            self.settings.max_group_lack = int(preset_file_values[6])

    @staticmethod
    def test_preset_file(preset_file_lines: List[str]) -> None:
        """Проверка того, правильно ли написан файл с настроками пресета"""
        neccessary_string_parts = [
            ("ProteinPilot summary analyzer", "Program name"),
            ("#Protein filter", "Protein filter header"),
            ("FDR", "Global FDR value"),
            ("ID exclusion list", "ID blacklist"),
            ("Peptide confidence", "Protein group filter peptide confidence"),
            (
                "Protein grouping (conf)",
                "Protein group filter peptide confidence",
            ),
            ("#Peptide filter", "Peptide filter header"),
            ("Peptide confidence", "Peptide filter peptide confidence"),
            ("#Output filter", "Output filter header"),
            ("Min groups", "Min groups"),
            ("Max missing values", "Max missing values per group"),
        ]

        for i, (string_part, description) in enumerate(
            neccessary_string_parts
        ):
            if string_part.lower() not in preset_file_lines[i].lower():
                raise IOError(description)

    def test_preset_file_existance(self) -> None:
        """Проверяет, существует ли файл preset.txt в папке с пресетом.

        Raises:
            FileNotFoundError: вызывается, если файл не существует.
        """
        if not path.exists(path.join(self.folder, "preset.txt")):
            self.error_code = 1
            raise FileNotFoundError(
                f"ERROR!!! preset.txt file not found in preset {self.folder}"
            )

    def test_result(self) -> None:
        """Проверка выходных файлов"""
        for filename in listdir(self.preset_output_dir):
            try:
                self.test_output_file(filename)
            except FileNotFoundError:
                print(f'Error! File "{filename}" not found in Output dir!')
                self.error_code = 1

    def test_output_file(self, filename: str) -> None:
        """Проверка одного выходного файла

        Args:
            filename: имя файла для проверки
        """
        paths = (
            path.join(self.preset_output_dir, filename),
            path.join(self.settings.output_path, filename),
        )
        with open(paths[0], encoding="utf-8") as file1:
            with open(paths[1], encoding="utf-8") as file2:
                files_content = (
                    file1.read().strip().split("\n"),
                    file2.read().strip().split("\n"),
                )
        if len(files_content[0]) != len(files_content[1]):
            print(f"{filename} files have different lengths")
            self.error_code = 2
        else:
            for i in range(0, len(files_content[0])):
                if files_content[0][i] != files_content[1][i]:
                    print(f"Line {i+1} is not equals in file {filename}")
                    self.error_code = 3
                    break

    def clear_output_folder(self):
        """Очищает выходную папку"""
        if path.exists(self.settings.output_path):
            for filename in listdir(self.settings.output_path):
                remove(path.join(self.settings.output_path, filename))


def get_presets_folders() -> List[str]:
    """Получает список папок с пресетами

    Args:
        presetFolder: папка, в которой расположены пресеты

    Returns:
        Список путей к папкам с пресетами
    """
    presets_folders = []
    for folder_name in listdir(PRESETS_FOLDER):
        if folder_name.startswith("preset"):
            presets_folders.append(path.join(PRESETS_FOLDER, folder_name))
    return sorted(presets_folders, key=lambda x: x.lower())


def main():
    """Main function (Program entrypoint)."""

    error_code = 0
    presets = ["FDRdefault", "IDexcl", "MinMax", "pepfilter"]
    if len(argv) == 2:
        folder = argv[1]
        preset = Preset(folder)
        preset.read_settings()
        preset.clear_output_folder()
        preset.run()
        if preset.error_code:
            error_code = preset.error_code
    else:
        presets_folders = get_presets_folders()
        for folder in presets_folders:
            if any((p.lower() in folder.lower() for p in presets)):
                try:
                    preset = Preset(folder)
                    preset.read_settings()
                    preset.clear_output_folder()
                    preset.run()
                    if preset.error_code:
                        error_code = preset.error_code
                    print()
                except Exception:
                    print(
                        f"Unexpected error catched in preset {preset.folder}"
                    )
                    raise
    sys_exit(error_code)


if __name__ == "__main__":
    main()
