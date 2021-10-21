from Classes.DB.db import DB
from Classes.Functions import GetInput
from Classes.Input import Input
from Classes.PeptideColumns import PeptideColumns
from Classes.RawPeptideTables import RawPeptideTables


def main(inputParams: Input = None):
    if inputParams is None:
        inputParams = GetInput()
    peptideTables = RawPeptideTables(PeptideColumns(), inputParams.inputPath)
    with DB(inputParams=inputParams) as db:
        # Заполняем таблицы
        db.fillers.fillSequence(inputParams.seqDB)
        db.fillers.fillRawPeptide(peptideTables)
        db.functions.testFastaDatabase()

        if inputParams.getFDRStr() == "default":
            db.fdr.default()
        if inputParams.blackList:
            db.fillers.fillExclusion(inputParams.blackList[1])
            db.functions.applyExclusion()

        if inputParams.isProteinConfidence:
            if inputParams.proteinConfidence is None:
                db.functions.applyPeptideConfidenceDefault()
            else:
                db.functions.applyPeptideConfidenceValue(inputParams.proteinConfidence)

        db.proteinGrouping.createFiltredAccessionTable(
            inputParams.proteinGroupingConfidence
        )
        db.proteinGrouping.createCountTable()
        db.proteinGrouping.createGeneralCountTable()
        db.proteinGrouping.applyProteinGrouping()
        db.proteinGrouping.fillPeptideTable()

        db.summaries.fillPeptideWithSum()

        db.functions.applyPeptideConfidenceValue2(inputParams.confPeptide)

        if (
            inputParams.maxGroupLack is not None
            and inputParams.minGroupsWithAccession is not None
        ):
            db.outputGrouping.applyGroupFilter(
                inputParams.maxGroupLack, inputParams.minGroupsWithAccession
            )

        db.output.GenerateOutputFiles()


if __name__ == "__main__":
    main()
