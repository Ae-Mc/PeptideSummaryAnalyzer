from typing import List, Optional


class ProteinGroup:
    representativeAccession: Optional[str]
    accessions: List[str]
    unused: float

    def __init__(
            self,
            unused: float = None,
            accessions: List[str] = None,
            representativeAccession: str = None) -> None:

        if unused is not None:
            self.unused = unused
        if accessions is not None:
            self.accessions = accessions
        self.representativeAccession = representativeAccession

    def __repr__(self):
        return (f"{self.representativeAccession} "
                f"({self.unused}): {self.accessions}")

    def __str__(self):
        return (f"{self.representativeAccession} "
                f"({self.unused}): {self.accessions}")

    def __eq__(self, other):
        if not isinstance(other, ProteinGroup):
            return NotImplemented
        return self.__dict__ == other.__dict__
