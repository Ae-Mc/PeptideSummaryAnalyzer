"""Содержит функции, которые нельзя объеденить в общий класс."""

from os import listdir, path
from sys import argv
from typing import List, Union

from classes.sequence_database import SequenceDatabase
from classes.input import Input


def get_file_lines(filename: str) -> Union[List[str], None]:
    """Считывает файл и возвращает список строк без символа переноса строки

    Args:
        filename: Имя файла, из которого происходит считывание строк

    Returns:
        Список строк файла, без символов переноса строки
    """

    if not (
        filename.endswith("\\")
        or filename.endswith("/")
        or len(filename.strip()) == 0
    ):
        with open(filename, encoding="utf-8") as tfile:
            return tfile.read().split("\n")
    return None


def find_fasta_file(input_path: str) -> str:
    """Ищет fasta файл в папке inputPath

    Args:
        input_path: путь к папке, в которой будет производиться поиск fasta
            файла

    Returns:
        Путь к fasta файлу, составленный из inputPath и имени fasta файла
    """
    for file in listdir(input_path):
        if file.endswith(".fasta"):
            return path.join(input_path, file)
    raise FileNotFoundError(f'Fasta file not found in folder "{input_path}"!')


def get_input() -> Input:
    """Получение параметров для запуска обработки

    Returns:
        Класс Input, содержащий все нужные параметры для запуска обработки
    """
    input_params = Input()
    input_params.root_path = "."
    input_params.inputPath = "./Input"
    input_params.seq_db = SequenceDatabase.from_file(
        find_fasta_file(input_params.root_path)
    )
    if len(argv) == 9:
        input_params.set_fdr(argv[1])
        black_list_lines = get_file_lines(argv[2])
        input_params.exclusion_list = (
            (argv[2], black_list_lines)
            if black_list_lines is not None
            else None
        )
        input_params.set_protein_confidence(argv[3])
        input_params.set_protein_grouping_confidence(argv[4])
        input_params.set_conf_peptide(argv[5])
        input_params.min_groups_with_accession = int(argv[6])
        input_params.max_group_lack = int(argv[7])
    else:
        print('"ProteinPilot summary analyzer"')
        print("#Protein filter")
        input_params.set_fdr(
            input("Global FDR critical value (<% k or default): ")
        )
        black_list_file = input("ID exclusion list: ")
        black_list_lines = get_file_lines(black_list_file)
        input_params.exclusion_list = (
            (black_list_file, black_list_lines)
            if black_list_lines is not None
            else None
        )
        input_params.set_protein_confidence(
            input("Peptide confidence (value or default): ")
        )
        input_params.set_protein_grouping_confidence(
            input("Protein grouping (conf): ")
        )
        print("#Peptide filter")
        input_params.set_conf_peptide(input("Peptide confidence (value): "))
        print("#Output filter")
        input_params.min_groups_with_accession = int(
            input("Min groups with ID: ")
        )
        input_params.max_group_lack = int(
            input("Max missing values per group: ")
        )
    return input_params
