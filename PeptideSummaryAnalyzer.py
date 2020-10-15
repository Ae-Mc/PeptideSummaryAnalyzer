#!/usr/bin/env python3
from decimal import FloatOperation, getcontext

from Classes import AccessionTables
from Classes import ColumnNames
from Classes import (ApplyBlackList, ApplyConfidenceIDFilter,
                     ApplyGroupFilter, ApplyParamsFilter,
                     ApplyWhiteList, CalculateAccessionsNormRatios,
                     GetInput, GetScPsigAndNormFilesSumm)
from Classes import Input
from Classes import Output
from Classes import PeptideTables
from Classes import ProteinGroupsDB
from Classes import ProteinPerTableList

"""Peptide Summary Analyzer

См. README.md
"""


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
        proteinGroupsDB = ProteinGroupsDB(inputParams.seqDB,
                                          inputParams.skipReversedIfSecondary,
                                          inputParams.inputPath)
        proteinPerTableList = ProteinPerTableList(proteinGroupsDB)
        peptideTables.ApplyProteinPerTableList(proteinPerTableList)
        peptideTables.ApplyProteinReplacements(
            proteinGroupsDB.GetReplacementsPerTable())

    if inputParams.blackList:
        ApplyBlackList(peptideTables,
                       inputParams.blackList,
                       columnNames)

    if inputParams.isConfID:
        ApplyConfidenceIDFilter(inputParams.confID,
                                peptideTables,
                                columnNames)
    ApplyParamsFilter(inputParams.unused,
                      inputParams.contrib,
                      inputParams.confPeptide,
                      peptideTables,
                      columnNames)

    accessionTables = AccessionTables(inputParams.seqDB,
                                      peptideTables,
                                      columnNames=columnNames)
    accessionTables.sortedTableNums = peptideTables.GetSortedTableNums()
    filesSumms = GetScPsigAndNormFilesSumm(accessionTables)
    CalculateAccessionsNormRatios(accessionTables, filesSumms)

    if inputParams.whiteList:
        ApplyWhiteList(peptideTables,
                       inputParams.whiteList,
                       columnNames)

    ApplyGroupFilter(accessionTables,
                     inputParams.maxGroupAbsence,
                     inputParams.minGroupsWithAccession)

    Output(inputParams.outputPath,
           seqDB=inputParams.seqDB,
           accessionTables=accessionTables,
           proteinGroupsDB=proteinGroupsDB)


if __name__ == "__main__":
    main()
