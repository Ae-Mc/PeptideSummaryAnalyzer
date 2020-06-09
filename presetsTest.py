#!/bin/env python
from os import listdir, path
from typing import List
from Classes.Input import Input
from Classes.Comparable import Comparable
from PeptideSummaryAnalyzer import ReadSeqDB, GetFileLines
from PeptideSummaryAnalyzer import main as proteinMain


presetsFolder = "Presets"
INPUTPATH = ""


def TestPresetFile(presetFileLines: List[str]):
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


def GetSettingsFromFolder(folder: str):
    settings: Input = Input()
    with open(path.join(folder, "preset.txt")) as presetFile:
        presetFileLines = presetFile.read().split("\n")
        try:
            TestPresetFile(presetFileLines)
        except IOError as e:
            print(f"Error reading preset {path.split(folder)[1]} "
                  f"{e.args[0]}")
            exit()
        presetFileValues = [
            line.split(':')[1].strip() for line in presetFileLines
            if len(line.strip()) > 0]

    settings.proteinPilotVersion = presetFileValues[0]
    settings.blackList = GetFileLines(presetFileValues[1])
    settings.whiteList = GetFileLines(presetFileValues[2])
    settings.isProteinGroupFilter = presetFileValues[3].lower()
    settings.seqDB = ReadSeqDB(path.join(folder, presetFileValues[4]))
    settings.unused = Comparable(presetFileValues[5])
    settings.contrib = Comparable(presetFileValues[6])
    settings.conf = presetFileValues[7]
    settings.minGroupsWithAccession = int(presetFileValues[8])
    settings.maxGroupAbsence = int(presetFileValues[9])
    return settings


def RunPreset(presetFolder: str):
    settings: Input = GetSettingsFromFolder(presetFolder)
    settings.inputPath = path.join(presetFolder, "Input")
    try:
        print(presetFolder)
        proteinMain(settings)
    except Exception as e:
        print(f"Error running preset {presetFolder}")
        raise e


def GetPresetsFolders(presetFolder: str):
    presetsFolders = []
    for folderName in listdir(presetsFolder):
        if folderName.startswith("preset"):
            presetsFolders.append(path.join(presetsFolder, folderName))
    return presetsFolders


def main():
    presetsFolders = GetPresetsFolders(presetsFolder)
    RunPreset(presetsFolders[4])
    # for folder in presetsFolders:
    #     RunPreset(folder)


if __name__ == "__main__":
    main()
