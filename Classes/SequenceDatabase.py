from .Sequence import Sequence


class SequenceDatabase(dict):
    """База последовательностей из файла БД с последовательностями в виде
    словаря Sequence вида {"Accession": Sequence}
    """

    @staticmethod
    def fromFile(seqDBFilename: str) -> "SequenceDatabase":
        """Считывание последовательностей из файла

        Считывание длин последовательностей из файла БД с последовательностями
        в словарь классов Sequence вида {"Accession": Sequence}

        Args:
            seqDBFilename: Имя файла с последовательностями
        Returns:
            Словарь вида {"Accession": Sequence}
        """

        with open(seqDBFilename) as seqDBFile:
            strings = seqDBFile.read().strip().split("\n")
            seqDBFile.close()
            seqDB = SequenceDatabase()
            i = 0
            while i < len(strings):
                if len(strings[i].strip()) > 0:
                    if strings[i].strip()[0] == ">":
                        seqID = strings[i].strip().split(" ")[0][1:]
                        seqDB[seqID] = Sequence(accession=seqID)
                        if len(strings[i].strip().split(" ")) > 1:
                            seqDB[seqID].desc = strings[i].strip().split(" ")[1]
                            for word in strings[i].strip().split(" ")[2:]:
                                if word.startswith("OS="):
                                    break
                                seqDB[seqID].desc += " " + word
                        i += 1
                        while (i < len(strings)) and strings[i].strip()[0] != ">":
                            seqDB[seqID].rawSeq += strings[i].strip()
                            i += 1
                        i -= 1
                        if not seqDB[seqID].len:
                            input("Error! Length of sequence with id " f"{seqID} = 0")
                            raise (IndexError)
                i += 1

            return seqDB

    def __getitem__(self, key) -> Sequence:
        if key not in self:
            raise KeyError(f"ERROR! Missing sequence: {key}")
        return super().__getitem__(key)

    def __setitem__(self, k: str, v: Sequence) -> None:
        return super().__setitem__(k, v)
