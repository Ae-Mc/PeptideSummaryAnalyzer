#!/bin/env python
from os import listdir, path, remove
from typing import List
from Classes.Input import Input
from Classes.Comparable import Comparable
from Classes.Functions import ReadSeqDB, GetFileLines
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
        try:
            print(self.folder)
            proteinMain(self.settings)
        except Exception as e:
            print(f"Error running preset {self.folder}")
            raise e
        self.TestResult()

    def ReadSettings(self) -> None:
        """Считывание настроек пресета из папки с пресетом"""
        self.settings: Input = Input()
        self.settings.inputPath = path.join(self.folder, "Input")
        self.settings.outputPath = path.join(self.folder, "Output")
        with open(path.join(self.folder, "preset.txt")) as presetFile:
            presetFileLines = presetFile.read().split("\n")
            try:
                self.TestPresetFile(presetFileLines)
            except IOError as e:
                print(f"Error reading preset {path.split(self.folder)[1]} "
                      f"{e.args[0]}")
                exit()
            presetFileValues = [
                line.split(':')[1].strip() for line in presetFileLines
                if len(line.strip()) > 0]

        self.settings.proteinPilotVersion = presetFileValues[0]
        self.settings.blackList = None
        self.settings.whiteList = None
        if len(presetFileValues[1].strip()):
            self.settings.whiteList = GetFileLines(
                path.join(self.folder, presetFileValues[1]))
        if len(presetFileValues[2].strip()):
            self.settings.blackList = GetFileLines(
                path.join(self.folder, presetFileValues[2]))
        self.settings.isProteinGroupFilter = presetFileValues[3]
        self.settings.skipReversedIfSecondary = presetFileValues[4]
        self.settings.seqDB = ReadSeqDB(
            path.join(self.folder, presetFileValues[5]))
        self.settings.unused = Comparable(presetFileValues[6])
        self.settings.contrib = Comparable(presetFileValues[7])
        self.settings.confID = presetFileValues[8]
        self.settings.confPeptide = presetFileValues[9]
        self.settings.minGroupsWithAccession = int(presetFileValues[10])
        self.settings.maxGroupAbsence = int(presetFileValues[11])

    def TestPresetFile(self, presetFileLines: List[str]) -> None:
        """Проверка того, правильно ли написан файл с настроками пресета"""
        neccessaryStringParts = [
            ("ProteinPilot", "ProteinPilot Version"),
            ("ID list file", "ID whitelist"),
            ("ID exclusion list", "ID blacklist"),
            ("Protein group filter", "Protein group filter"),
            ("Skip REVERSED", "Skip REVERSED if secondary is not REVERSED"),
            ("Database", "Database"),
            ("Unused", "Unused"),
            ("Contribution", "Contribution"),
            ("Confidence ID", "Confidence ID"),
            ("Confidence peptide", "Confidence peptide"),
            ("Min groups", "Min groups"),
            ("Max missing values", "Max missing values per group")
        ]

        i = 0
        for stringPart, description in neccessaryStringParts:
            if stringPart.lower() not in presetFileLines[i].lower():
                raise IOError(description)
            i += 1

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
