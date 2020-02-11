from typing import Dict, List


def ReadTable(tableFilename: str, sep='\t') -> Dict[str, List[str]]:
    """ Считываем файл в словарь, где ключом является заголовок, а
    значением — список значений в столбце """

    with open(tableFilename) as tempFile:
        strings = tempFile.read().split('\n')
        tempFile.close()
        table: Dict[str, List[str]] = {}
        columns = strings[0].split(sep)

        for column in columns:
            table[column] = []
        for string in strings[1:]:
            if len(string.strip()):
                i = 0
                for value in string.split(sep):
                    table[columns[i]].append(value)
                    i += 1

        return table
    return None
