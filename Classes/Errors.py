class AccessionNotFoundError(Exception):

    def __init__(self, message):
        self.message = message


class ColumnNotFoundError(Exception):

    def __init__(self, columnName: str, filename: str = None):
        if columnName is not None:
            self.message = f"Column {columnName} not found!"
            if filename is not None:
                self.message = (
                    f"Column {columnName} not found in file {filename}!")
        else:
            self.message = None


class RepresentativeAccessionNotFoundError(Exception):

    def __init__(self, group):
        if group is not None:
            self.message = (
                "Representative accession for group with accessions "
                f"{group.accessions} not found!"
            )
        else:
            self.message = None


class EmptyGroupError(Exception):

    def __init__(self, message):
        self.message = message
