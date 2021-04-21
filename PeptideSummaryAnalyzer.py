#!/usr/bin/env python3
from decimal import FloatOperation, getcontext

from Classes import AccessionTables
from Classes import (ApplyBlackList, ApplyConfidenceIDFilter,
                     ApplyGroupFilter, ApplyParamsFilter,
                     CalculateAccessionsNormRatios, GetInput,
                     GetScPsigAndNormFilesSumm)
from Classes import Input
from Classes import Output
from Classes import PeptideTables
from Classes import ProteinAccessionsDB
from Classes import ProteinGroupsDB
from Classes import ProteinPerTableList
from Classes import PeptideColumns

"""См. README"""


def main(inputParams: Input = None) -> None:
    getcontext().traps[FloatOperation] = True

    if inputParams is None:
        inputParams = GetInput()
    columnNames = PeptideColumns()

    peptideTables = PeptideTables(columnNames, inputDir=inputParams.inputPath)
    proteinGroupsDB = None
    if inputParams.isProteinGroupFilter is not None:
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

    if inputParams.blackList is not None:
        ApplyBlackList(peptideTables,
                       inputParams.blackList[1])

    if inputParams.isConfID:
        ApplyConfidenceIDFilter(inputParams.confID, peptideTables)
    ApplyParamsFilter(inputParams.unused,
                      inputParams.confPeptide,
                      peptideTables)

    accessionTables = AccessionTables(inputParams.seqDB, peptideTables)
    accessionTables.sortedTableNums = peptideTables.GetSortedTableNums()
    filesSumms = GetScPsigAndNormFilesSumm(accessionTables)
    CalculateAccessionsNormRatios(accessionTables, filesSumms)

    ApplyGroupFilter(accessionTables,
                     inputParams.maxGroupAbsence,
                     inputParams.minGroupsWithAccession)

    Output(inputParams,
           seqDB=inputParams.seqDB,
           accessionTables=accessionTables,
           proteinGroupsDB=proteinGroupsDB)


if __name__ == "__main__":
    main()
