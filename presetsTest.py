#!/bin/env python
from traceback import print_exc
from os import listdir, path, remove
from typing import List
from Classes.Input import Input
from Classes.Comparable import Comparable
from Classes.Functions import ReadSeqDB, GetFileLines, FindFastaFile
from PeptideSummaryAnalyzer import main as proteinMain


presetsFolder = "Presets"
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
    presetOutputDir: str
    settings: Input
    errorCode: int = 0

    def __init__(self, folder: str):
        """Args:
            folder: папка с пресетом
        """
        self.folder = folder
        self.presetOutputDir = path.join(self.folder, "TrueOutput")

    def Run(self) -> None:
        """Запуск теста пресета"""
        if self.errorCode == 0:
            try:
                print(self.folder)
                proteinMain(self.settings)
            except Exception:
                print(f"Error running preset {self.folder}")
                print_exc()
            self.TestResult()

    def ReadSettings(self) -> None:
        """Считывание настроек пресета из папки с пресетом"""
        self.settings: Input = Input()
        self.settings.inputPath = path.join(self.folder, "Input")
        self.settings.outputPath = path.join(self.folder, "Output")
        try:
            self.TestPresetFileExistance()
        except FileNotFoundError as e:
            print(e.args[0])
            return

        with open(path.join(self.folder, "preset.txt")) as presetFile:
            presetFileLines = presetFile.read().split("\n")
            try:
                self.TestPresetFile(presetFileLines)
            except IOError as e:
                print(f"Error reading preset {path.split(self.folder)[1]} "
                      f"{e.args[0]}")
                self.errorCode = 1
                return
            presetFileValues = [
                line.split(':')[1].strip() for line in presetFileLines
                if len(line.strip()) > 0]

        self.settings.blackList = None
        if len(presetFileValues[0].strip()):
            self.settings.blackList = GetFileLines(
                path.join(self.folder, presetFileValues[0]))
        self.settings.isProteinGroupFilter = presetFileValues[1].lower()
        self.settings.seqDB = ReadSeqDB(FindFastaFile(self.folder))
        self.settings.unused = Comparable(presetFileValues[2])
        self.settings.confID = presetFileValues[3]
        self.settings.confPeptide = presetFileValues[4]
        self.settings.minGroupsWithAccession = int(presetFileValues[5])
        self.settings.maxGroupAbsence = int(presetFileValues[6])

    def TestPresetFile(self, presetFileLines: List[str]) -> None:
        """Проверка того, правильно ли написан файл с настроками пресета"""
        neccessaryStringParts = [
            ("ID exclusion list", "ID blacklist"),
            ("Protein group filter", "Protein group filter"),
            ("Unused", "Unused"),
            ("Confidence ID", "Confidence ID"),
            ("Confidence peptide", "Confidence peptide"),
            ("Min groups", "Min groups"),
            ("Max missing values", "Max missing values per group")
        ]

        for i, (stringPart, description) in enumerate(neccessaryStringParts):
            if stringPart.lower() not in presetFileLines[i].lower():
                raise IOError(description)

    def TestPresetFileExistance(self) -> None:
        if not path.exists(path.join(self.folder, "preset.txt")):
            self.errorCode = 1
            raise FileNotFoundError(
                f"ERROR!!! preset.txt file not found in preset {self.folder}")

    def TestResult(self) -> None:
        """Проверка выходных файлов"""
        for filename in listdir(self.presetOutputDir):
            try:
                self.TestOutputFile(filename)
            except FileNotFoundError:
                print(f"Error! File \"{filename}\" not found in Output dir!")
                self.errorCode = 1

    def TestOutputFile(self, filename: str) -> None:
        """Проверка одного выходного файла

        Args:
            filename: имя файла для проверки
        """
        paths = (path.join(self.presetOutputDir, filename),
                 path.join(self.settings.outputPath, filename))
        with open(paths[0]) as file1:
            with open(paths[1]) as file2:
                filesContent = (file1.read().strip().split("\n"),
                                file2.read().strip().split("\n"))
        if(len(filesContent[0]) != len(filesContent[1])):
            print(f"{paths[0]} != {paths[1]}")
            self.errorCode = 2
        else:
            for i in range(0, len(filesContent[0])):
                if(filesContent[0][i] != filesContent[1][i]):
                    print(f"Line {i+1} is not equals in file {filename}")
                    self.errorCode = 3
                    break

    def ClearOutputFolder(self):
        """Очищает выходную папку"""
        if path.exists(self.settings.outputPath):
            for filename in listdir(self.settings.outputPath):
                remove(path.join(self.settings.outputPath, filename))


def GetPresetsFolders(presetFolder: str) -> List[str]:
    """Получает список папок с пресетами

    Args:
        presetFolder: папка, в которой расположены пресеты

    Returns:
        Список путей к папкам с пресетами
    """
    presetsFolders = []
    for folderName in listdir(presetsFolder):
        if folderName.startswith("preset"):
            presetsFolders.append(path.join(presetsFolder, folderName))
    return sorted(presetsFolders, key=lambda x: x.lower())


def main():
    presetsFolders = GetPresetsFolders(presetsFolder)
    errorCode = 0
    for folder in presetsFolders:
        preset = Preset(folder)
        preset.ReadSettings()
        preset.ClearOutputFolder()
        preset.Run()
        if preset.errorCode:
            errorCode = preset.errorCode
        print()
    exit(errorCode)


if __name__ == "__main__":
    main()
