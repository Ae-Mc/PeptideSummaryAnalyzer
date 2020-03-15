class AccessionNotFoundError(Exception):

    def __init__(self, message):
        self.message = message


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
