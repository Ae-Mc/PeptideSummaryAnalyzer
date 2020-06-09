#!/bin/env python
from os import listdir, path, remove
from typing import List
from Classes.Input import Input
from Classes.Comparable import Comparable
from PeptideSummaryAnalyzer import ReadSeqDB, GetFileLines
from PeptideSummaryAnalyzer import main as proteinMain


presetsFolder = "Presets"
INPUTPATH = ""


class Preset:
    folder: str
    presetOutputDir: str
    settings: Input

    def __init__(self, folder: str):
        self.folder = folder
        self.presetOutputDir = path.join(self.folder, "TrueOutput")

    def Run(self) -> None:
        try:
            print(self.folder)
            proteinMain(self.settings)
        except Exception as e:
            print(f"Error running preset {self.folder}")
            raise e
        self.TestResult()

    def ReadSettings(self) -> None:
        self.settings: Input = Input()
        self.settings.inputPath = path.join(self.folder, "Input")
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
        self.settings.blackList = GetFileLines(presetFileValues[1])
        self.settings.whiteList = GetFileLines(presetFileValues[2])
        self.settings.isProteinGroupFilter = presetFileValues[3].lower()
        self.settings.seqDB = ReadSeqDB(
            path.join(self.folder, presetFileValues[4]))
        self.settings.unused = Comparable(presetFileValues[5])
        self.settings.contrib = Comparable(presetFileValues[6])
        self.settings.conf = presetFileValues[7]
        self.settings.minGroupsWithAccession = int(presetFileValues[8])
        self.settings.maxGroupAbsence = int(presetFileValues[9])

    def TestPresetFile(self, presetFileLines: List[str]) -> None:
        neccessaryStringParts = [
            ("ProteinPilot", "ProteinPilot Version"),
            ("ID list file", "ID whitelist"),
            ("ID exclusion list", "ID blacklist"),
            ("Protein group filter", "Protein group filter"),
            ("Database", "Database"),
            ("Unused", "Unused"),
            ("Contribution", "Contribution"),
            ("Confidence", "Confidence"),
            ("Min groups", "Min groups"),
            ("Max missing values", "Max missing values per group")
        ]

        i = 0
        for stringPart, description in neccessaryStringParts:
            if stringPart.lower() not in presetFileLines[i].lower():
                raise IOError(description)
            i += 1

    def TestResult(self):
        for filename in listdir(self.presetOutputDir):
            try:
                self.TestOutputFile(filename)
            except FileNotFoundError:
                print(f"Error! File \"{filename}\" not found in Output dir!")

    def TestOutputFile(self, filename: str) -> None:
        paths = (path.join(self.presetOutputDir, filename),
                 path.join("Output", filename))
        with open(paths[0]) as file1:
            with open(paths[1]) as file2:
                filesContent = (file1.read().split("\n"),
                                file2.read().split("\n"))
        if(len(filesContent[0]) != len(filesContent[1])):
            print(f"{paths[0]} != {paths[1]}")

    @staticmethod
    def ClearOutputFolder():
        for filename in listdir("Output"):
            remove(path.join("Output", filename))


def GetPresetsFolders(presetFolder: str):
    presetsFolders = []
    for folderName in listdir(presetsFolder):
        if folderName.startswith("preset"):
            presetsFolders.append(path.join(presetsFolder, folderName))
    return presetsFolders


def main():
    presetsFolders = GetPresetsFolders(presetsFolder)
    # RunPreset(presetsFolders[4])
    for folder in presetsFolders:
        preset = Preset(folder)
        Preset.ClearOutputFolder()
        preset.ReadSettings()
        preset.Run()


if __name__ == "__main__":
    main()
