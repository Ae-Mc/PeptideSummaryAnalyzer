#!/usr/bin/env python3
from decimal import FloatOperation, getcontext

from Classes import AccessionTables
from Classes import (ApplyBlackList, ApplyConfidenceIDFilter,
                     ApplyGroupFilter, ApplyParamsFilter,
                     ApplyWhiteList, CalculateAccessionsNormRatios, GetInput,
                     GetScPsigAndNormFilesSumm)
from Classes import Input
from Classes import Output
from Classes import PeptideTables
from Classes import ProteinAccessionsDB
from Classes import ProteinGroupsDB
from Classes import ProteinPerTableList
from Classes import PeptideColumns, PeptideColumns5

"""Peptide Summary Analyzer

Peptide summary analyzer — это программа для обработки массивов файлов,
полученных c помощью программы ProteinPilot (Sciex) в ходе протеомных
исследований. Используется в Институте Цитологии РАН группой протеомики и
масс-спектрометрии.

Использование:
    Программа поддерживает как диалоговый вид, так и работу через передачу
    аргументов командной строки. Для передачи аргументов нужно сохранять
    следующий порядок их следования:
        - Версия ProteinPilot (4 или 5)
        - Файл белого списка
        - Файл чёрного списка
        - Нужно ли применять Protein group фильтр (Y или N)
        - Файл бд для подсчёта длин последовательностей
        - Unused параметр
        - Contribution параметр
        - Confidence ID параметр
        - Confidence peptide параметр
        - Минимальное количество групп с ID
        - Максимальное количество отсутствий ID в каждой группе

    Параметры Unused, Contribution и Confidence peptide поддерживают следующие
    форматы (вместо 10 может быть любое число):
        - < 10 — меньше 10
        - > 10 — больше 10
        - >= 10 — больше либо равно 10
        - <= 10 — меньше либо равно 10
        - == 10 — равно 10
        - != 10 — не равно 10

Пример работы (`""` — пустой аргумент):
    `python PeptideSummaryAnalyzer.py 5 "" "" N test.fasta "" ">= 90" ""
        ">= 95" 1 1`
"""


def main(inputParams: Input = None) -> None:
    getcontext().traps[FloatOperation] = True

    if inputParams is None:
        inputParams = GetInput()
    columnNames = PeptideColumns()
    if inputParams.proteinPilotVersion == '5':
        columnNames = PeptideColumns5()

    peptideTables = PeptideTables(columnNames, inputDir=inputParams.inputPath)
    proteinGroupsDB = None
    if inputParams.isProteinGroupFilter:
        proteinTables = ProteinAccessionsDB.GetProteinTables(
            inputParams.inputPath)
        proteinAccessionsDB = ProteinAccessionsDB(proteinTables)
        proteinGroupsDB = ProteinGroupsDB(proteinAccessionsDB,
                                          inputParams.seqDB,
                                          proteinTables)
        proteinPerTableList = ProteinPerTableList(proteinGroupsDB)
        peptideTables.ApplyProteinPerTableList(proteinPerTableList)
        peptideTables.ApplyProteinReplacements(
            proteinGroupsDB.GetReplacementsPerTable())

    if inputParams.blackList:
        ApplyBlackList(peptideTables,
                       inputParams.blackList)

    if inputParams.isConfID:
        ApplyConfidenceIDFilter(inputParams.confID, peptideTables)
    ApplyParamsFilter(inputParams.unused,
                      inputParams.contrib,
                      inputParams.confPeptide,
                      peptideTables)

    accessionTables = AccessionTables(inputParams.seqDB, peptideTables)
    accessionTables.sortedTableNums = peptideTables.GetSortedTableNums()
    filesSumms = GetScPsigAndNormFilesSumm(accessionTables)
    CalculateAccessionsNormRatios(accessionTables, filesSumms)

    if inputParams.whiteList:
        ApplyWhiteList(peptideTables, inputParams.whiteList)

    ApplyGroupFilter(accessionTables,
                     inputParams.maxGroupAbsence,
                     inputParams.minGroupsWithAccession)

    Output(inputParams.outputPath,
           seqDB=inputParams.seqDB,
           accessionTables=accessionTables,
           proteinGroupsDB=proteinGroupsDB)


if __name__ == "__main__":
    main()
