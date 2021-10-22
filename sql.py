"""Главный модуль. Отвечает за работу программы."""

from Classes.db import DB
from Classes.Functions import GetInput
from Classes.Input import Input
from Classes.PeptideColumns import PeptideColumns
from Classes.RawPeptideTables import RawPeptideTables


def main(input_params: Input = None):
    """Main function (program entrypoint). Could be imported from testing
    module.

    Args:
        inputParams (Input, optional): Входные параметры. Если не заданы, то
            будут запрошены через GetInput.
    """
    if input_params is None:
        input_params = GetInput()
    peptide_tables = RawPeptideTables(PeptideColumns(), input_params.inputPath)
    with DB(input_params=input_params) as database:
        # Заполняем таблицы
        database.fillers.fill_sequence(input_params.seqDB)
        database.fillers.fill_raw_peptide(peptide_tables)
        database.functions.test_fasta_database()

        if input_params.getFDRStr() == "default":
            database.fdr.default()
        if input_params.blackList:
            database.fillers.fill_exclusion(input_params.blackList[1])
            database.functions.apply_exclusion()

        if input_params.isProteinConfidence:
            if input_params.proteinConfidence is None:
                database.functions.apply_peptide_confidence_default()
            else:
                database.functions.apply_peptide_confidence_value(
                    input_params.proteinConfidence
                )

        database.protein_grouping.create_filtred_accession_table(
            input_params.proteinGroupingConfidence
        )
        database.protein_grouping.fill_count_table()
        database.protein_grouping.create_general_count_table()
        database.protein_grouping.apply_protein_grouping()
        database.protein_grouping.fill_peptide_table()

        database.summaries.fill_peptide_with_sum()

        database.functions.apply_peptide_confidence_value2(
            input_params.confPeptide
        )

        if (
            input_params.maxGroupLack is not None
            and input_params.minGroupsWithAccession is not None
        ):
            database.output_grouping.apply_group_filter(
                input_params.maxGroupLack, input_params.minGroupsWithAccession
            )

        database.output.generate_output_files()


if __name__ == "__main__":
    main()
