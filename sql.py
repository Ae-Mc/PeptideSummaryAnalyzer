"""Главный модуль. Отвечает за работу программы. Описание программы см. в
README.md."""

from classes import (
    FDRtype,
    Input,
    PeptideColumns,
    ProteinConfidenceType,
    RawPeptideTables,
    get_input,
)
from classes.db import DB


def main(input_params: Input = None):
    """Main function (program entrypoint). Could be imported from testing
    module.

    Args:
        inputParams (Input, optional): Входные параметры. Если не заданы, то
            будут запрошены через GetInput.
    """
    if input_params is None:
        input_params = get_input()
    peptide_tables = RawPeptideTables(PeptideColumns(), input_params.inputPath)
    with DB(input_params=input_params) as database:
        # Заполняем таблицы
        database.fillers.fill_sequence(input_params.seq_db)
        database.fillers.fill_raw_peptide(peptide_tables)
        database.functions.test_fasta_database()

        if input_params.fdr_type == FDRtype.NONE:
            database.fdr.none()
        elif input_params.fdr_type == FDRtype.DEFAULT:
            database.fdr.default()
        if input_params.exclusion_list:
            database.fillers.fill_exclusion(input_params.exclusion_list[1])
            database.functions.apply_exclusion()

        if (
            input_params.protein_confidence_type
            == ProteinConfidenceType.DEFAULT
        ):
            database.functions.apply_peptide_confidence_default()
        elif (
            input_params.protein_confidence_type == ProteinConfidenceType.VALUE
        ):
            database.functions.apply_peptide_confidence_value(
                input_params.protein_confidence
            )

        database.protein_grouping.create_filtred_accession_table(
            input_params.protein_grouping_confidence
        )
        database.protein_grouping.fill_count_table()
        database.protein_grouping.create_general_count_table()
        database.protein_grouping.apply_protein_grouping()
        database.protein_grouping.fill_peptide_table()

        database.functions.apply_peptide_confidence_value2(
            input_params.conf_peptide
        )

        database.summaries.fill_peptide_with_sum()

        if (
            input_params.max_group_lack is not None
            and input_params.min_groups_with_accession is not None
        ):
            database.output_grouping.apply_group_filter(
                input_params.max_group_lack,
                input_params.min_groups_with_accession,
            )

        database.output.generate_output_files()


if __name__ == "__main__":
    main()
