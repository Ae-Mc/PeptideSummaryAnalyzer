#!/usr/bin/env python3
from Classes.FDRFilter import FDRFilter
from Classes.RawPeptideTables import RawPeptideTables
from decimal import FloatOperation, getcontext

from Classes import AccessionTables
from Classes import (
    ApplyBlackList,
    ApplyPeptideConfidenceFilter,
    ApplyGroupFilter,
    ApplyProteinConfidenceFilter,
    CalculateAccessionsNormRatios,
    GetInput,
    GetScPsigAndNormFilesSumm,
    TestFastaAccessions,
)
from Classes import Input
from Classes import Output
from Classes import PeptideTables
from Classes import PeptideColumns

"""См. README"""


def main(inputParams: Input = None) -> None:
    getcontext().traps[FloatOperation] = True

    if inputParams is None:
        inputParams = GetInput()
    columnNames = PeptideColumns()

    rawPeptideTables = RawPeptideTables(
        columnNames, inputDir=inputParams.inputPath
    )

    # TODO: FDR filter
    FDRFilter(rawPeptideTables=rawPeptideTables).ApplyDefaultFilter()

    if inputParams.blackList is not None:
        ApplyBlackList(rawPeptideTables, inputParams.blackList[1])

    TestFastaAccessions(inputParams.seqDB, rawPeptideTables)
    peptideTables = PeptideTables(rawPeptideTables, seqDB=inputParams.seqDB)

    if inputParams.isProteinConfidence is True:
        ApplyProteinConfidenceFilter(
            inputParams.proteinConfidence, peptideTables
        )
    ApplyPeptideConfidenceFilter(inputParams.confPeptide, peptideTables)

    accessionTables = AccessionTables(inputParams.seqDB, peptideTables)
    accessionTables.sortedTableNums = peptideTables.GetSortedTableNums()
    filesSumms = GetScPsigAndNormFilesSumm(accessionTables)
    CalculateAccessionsNormRatios(accessionTables, filesSumms)

    ApplyGroupFilter(
        accessionTables,
        inputParams.maxGroupAbsence,
        inputParams.minGroupsWithAccession,
    )

    Output(
        inputParams,
        seqDB=inputParams.seqDB,
        accessionTables=accessionTables,
        proteinGroupsDB=None,
    )


if __name__ == "__main__":
    main()
