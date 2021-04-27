class ColumnNotFoundError(Exception):
    """Если столбец не был найден"""
    def __init__(self, columnName: str, filename: str = None):
        if columnName is not None:
            if filename is not None:
                self.message = (
                    f"Column {columnName} not found in file {filename}!")
            else:
                self.message = f"Column {columnName} not found!"
        else:
            self.message = None


class RepresentativeAccessionNotFoundError(Exception):
    """Если для Protein группы не был определён репрезентативный ID"""
    def __init__(self, group):
        if group is not None:
            self.message = (
                "Representative accession for group with accessions "
                f"{group.accessions} not found!"
            )
        else:
            self.message = None
