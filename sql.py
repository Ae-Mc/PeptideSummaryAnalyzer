"""Главный модуль. Отвечает за работу программы. Описание программы см. в
README.md."""

from classes import (
    FDRtype,
    Input,
    PeptideColumns,
    ProteinConfidenceType,
    ProteinColumns,
    RawTables,
    RawPeptideTable,
    RawProteinTable,
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
    peptide_tables = RawTables(
        input_params.inputPath,
        r"\d+\.\d+_.*DistinctPeptideSummary\.txt",
        PeptideColumns(),
        RawPeptideTable,
    )
    with DB(input_params=input_params) as database:
        # Заполняем таблицы
        database.fillers.fill_sequence(input_params.seq_db)
        database.fillers.fill_raw_peptide(peptide_tables)
        database.functions.test_fasta_database()

        if input_params.fdr_type == FDRtype.DEFAULT:
            database.fdr.default()
        elif input_params.fdr_type in [
            FDRtype.FDR_K_RANGE,
            FDRtype.FDR_KEST_RANGE,
        ]:
            protein_tables = RawTables(
                input_params.inputPath,
                r"\d+\.\d+_.*ProteinSummary\.txt",
                ProteinColumns(),
                RawProteinTable,
            )
            database.fillers.fill_protein(protein_tables)
            database.fdr.fill_main_accession()
            database.fdr.initialize_fdr_summary_table()
            database.fdr.fill_target_and_decoy_columns()

            database.fdr.fill_fdr_data(input_params.fdr_k or None)

            database.fdr.fill_a_and_b_columns(input_params.fdr_range)
            database.fdr.fill_additional_statistic_columns(
                input_params.fdr_range
            )
            database.fdr.apply_fdr_k(input_params.fdr_fdr)
        database.fdr.none()

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
