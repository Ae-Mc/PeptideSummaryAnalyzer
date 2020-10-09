#!/usr/bin/env python3
from decimal import FloatOperation, getcontext

from Classes.AccessionTables import AccessionTables
from Classes.ColumnNames import ColumnNames
from Classes.Functions import (ApplyBlackList, ApplyConfidenceIDFilter,
                               ApplyGroupFilter, ApplyParamsFilter,
                               ApplyWhiteList, CalculateAccessionsNormRatios,
                               GetInput, GetScPsigAndNormFilesSumm)
from Classes.Input import Input
from Classes.Output import Output
from Classes.PeptideTables import PeptideTables
from Classes.ProteinAccessionsDB import ProteinAccessionsDB
from Classes.ProteinGroupsDB import ProteinGroupsDB
from Classes.ProteinPerTableList import ProteinPerTableList

## @mainpage Peptide Summary Analyzer
# @section intro_sec Краткое описание
# Peptide summary analyzer — это программа для обработки массивов файлов,
# полученных c помощью программы ProteinPilot (Sciex) в ходе протеомных
# исследований. Используется в Институте Цитологии РАН группой протеомики и
# масс-спектрометрии.
# @section usage_sec Использование
# Программа поддерживает как диалоговый вид, так и работу через передачу
# аргументов командной строки. Для передачи аргументов нужно сохранять
# следующий порядок их следования:
# -# Версия ProteinPilot (4 или 5)
# -# Файл белого списка
# -# Файл чёрного списка
# -# Нужно ли применять Protein group фильтр (Y или N)
# -# Файл бд для подсчёта длин последовательностей
# -# Unused параметр
# -# Contribution параметр
# -# Confidence ID параметр
# -# Confidence peptide параметр
# -# Минимальное количество групп с ID
# -# Максимальное количество отсутствий ID в каждой группе
#
# Параметры Unused, Contribution и Confidence peptide поддерживают следующие
# форматы (вместо 10 может быть любое число):
# - < 10 — меньше 10
# - > 10 — больше 10
# - >= 10 — больше либо равно 10
# - <= 10 — меньше либо равно 10
# - == 10 — равно 10
# - != 10 — не равно 10
#
# Пример работы (`""` — пустой аргумент): @n
# <b> `python PeptideSummaryAnalyzer.py 5 "" "" N test.fasta "" ">= 90" ""
# ">= 95" 1 1`</b>


def main(inputParams: Input = None) -> None:
    getcontext().traps[FloatOperation] = True

    if inputParams is None:
        inputParams = GetInput()
    if inputParams.proteinPilotVersion == '5':
        columnNames = ColumnNames(precursorSignal="Intensity (Peptide)")
    else:
        columnNames = ColumnNames()

    peptideTables = PeptideTables(columnNames, inputDir=inputParams.inputPath)
    proteinGroupsDB = None
    if inputParams.isProteinGroupFilter:
        proteinAccessionsDB = ProteinAccessionsDB(inputParams.inputPath)
        proteinGroupsDB = ProteinGroupsDB(proteinAccessionsDB,
                                          inputParams.seqDB,
                                          inputParams.inputPath)
        proteinPerTableList = ProteinPerTableList(proteinGroupsDB)
        peptideTables.ApplyProteinPerTableList(proteinPerTableList)
        peptideTables.ApplyProteinReplacements(
            proteinGroupsDB.GetReplacementsPerTable())

    if inputParams.blackList:
        ApplyBlackList(peptideTables.peptideTables,
                       inputParams.blackList,
                       columnNames)

    if inputParams.isConfID:
        ApplyConfidenceIDFilter(inputParams.confID,
                                peptideTables.peptideTables,
                                columnNames)
    ApplyParamsFilter(inputParams.unused,
                      inputParams.contrib,
                      inputParams.confPeptide,
                      peptideTables.peptideTables,
                      columnNames)

    accessionTables = AccessionTables(inputParams.seqDB,
                                      peptideTables,
                                      columnNames=columnNames)
    accessionTables.sortedTableNums = peptideTables.GetSortedTableNums()
    filesSumms = GetScPsigAndNormFilesSumm(accessionTables.accessionsPerTable)

    accessionTables.GetAccessionsPerTable(inputParams.seqDB, peptideTables)
    CalculateAccessionsNormRatios(accessionTables, filesSumms)

    if inputParams.whiteList:
        ApplyWhiteList(peptideTables.peptideTables,
                       inputParams.whiteList,
                       columnNames)

    ApplyGroupFilter(accessionTables,
                     inputParams.maxGroupAbsence,
                     inputParams.minGroupsWithAccession)

    Output("Output/",
           seqDB=inputParams.seqDB,
           accessionTables=accessionTables,
           proteinGroupsDB=proteinGroupsDB)


if __name__ == "__main__":
    main()
