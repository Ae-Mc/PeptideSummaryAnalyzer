"""См. класс SequenceDatabase."""

from .sequence import Sequence


class SequenceDatabase(dict):
    """База последовательностей из файла БД с последовательностями в виде
    словаря Sequence вида {"Accession": Sequence}
    """

    @staticmethod
    def from_file(sequence_db_filename: str) -> "SequenceDatabase":
        """Считывание последовательностей из файла

        Считывание длин последовательностей из файла БД с последовательностями
        в словарь классов Sequence вида {"Accession": Sequence}

        Args:
            seqDBFilename: Имя файла с последовательностями
        Returns:
            Словарь вида {"Accession": Sequence}
        """

        with open(sequence_db_filename, encoding="utf-8") as sequence_db_file:
            strings = list(
                map(
                    lambda x: x.strip(),
                    sequence_db_file.read().strip().split("\n"),
                )
            )
            sequence_db_file.close()
            sequence_db = SequenceDatabase()
            i = 0
            while i < len(strings):
                if len(strings[i]) > 0:
                    if strings[i][0] == ">":
                        accession = strings[i].split(" ")[0][1:]
                        sequence_db[accession] = Sequence(accession=accession)
                        if len(strings[i].split(" ")) > 1:
                            sequence_db[accession].desc = strings[i].split(
                                " "
                            )[1]
                            for word in strings[i].split(" ")[2:]:
                                if word.startswith("OS="):
                                    break
                                sequence_db[accession].desc += " " + word
                        i += 1
                        while (i < len(strings)) and strings[i][0] != ">":
                            sequence_db[accession].raw_seq += strings[i]
                            i += 1
                        i -= 1
                        if not sequence_db[accession].len:
                            input(
                                "Error! Length of sequence with id "
                                f"{accession} = 0"
                            )
                            raise IndexError
                i += 1

            return sequence_db

    def __getitem__(self, key) -> Sequence:
        if key not in self:
            raise KeyError(f"ERROR! Missing sequence: {key}")
        return super().__getitem__(key)

    # pylint: disable=useless-super-delegation
    def __setitem__(self, k: str, v: Sequence) -> None:
        return super().__setitem__(k, v)
