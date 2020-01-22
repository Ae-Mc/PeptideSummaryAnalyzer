import sys
from os import listdir, mkdir
from os.path import exists

INPUTPATH = "Input"
OUTPUTPATH = "Output/"

# Номера столбцов, в которых содержаться нужные данные
positions = {
    "Unused": 1,
    "Id": 6,
    "Contribution": 10,
    "Confidence": 11,
    "Sequence": 12,
    "Sc": 21,
    "Precursor Signal": 24
}
# Список параметров, для которых создаются выходные файлы
params = ["Unused", "seq_length_summ", "counts", "Sc_summ", "Psignal_summ",
          "Sc_norm", "Psignal_norm", "SP_2"]
# Список параметров, присутствующих в output.txt
outParams = ["Unused", "seq_length_summ", "counts", "Sc_summ", "Psignal_summ",
             "Sc_norm", "Psignal_norm", "SP_2", "seq_length"]


# Класс, реализующий параметры базового id
class ID(dict):
    def __init__(self):
        for param in outParams:
            self[param] = 0


# Функция для проверки является ли строка числом
def isFloat(string):
    try:
        float(string)
        return True
    except ValueError:
        return False


# Ввод начальных данных
def inputParams(idFileName, idExclusionFileName, dbFileName, input_unused,
                input_contribution, input_confidence, input_minGroupsWithId,
                input_maxGroupAbsence):
    # Функция, получающая из строки метод сравнения и числовое значение
    def getOpData(string):
        string = string.strip()
        if not len(string):
            return "", None
        op = ""
        for ch in string:
            if ch in "!<>=":
                op += ch
            else:
                break
        if op == "":
            return "", None
        return op, float(string[len(op):])

    IDs = []
    if idFileName != "":
        tempFile = open(idFileName)
        IDs = tempFile.read().split("\n")
        tempFile.close()

    IDExclusionList = []
    if idExclusionFileName != "":
        tempFile = open(idExclusionFileName)
        IDExclusionList = tempFile.read().split("\n")
        tempFile.close()

    if dbFileName != "":
        tempFile = open(dbFileName)
        strings = tempFile.read().split("\n")
        tempFile.close()
        seqDB = {}  # Хранит длины последовательностей для каждого id
        seqDescDB = {}  # Хранит описания для каждого id
        i = 0
        while (i < len(strings)):
            if len(strings[i]):
                if strings[i][0] == '>':
                    seqID = strings[i].split(' ')[0][1:]
                    seqDB[seqID] = 0
                    if len(strings[i].split(" ")) > 1:
                        seqDescDB[seqID] = strings[i].split(" ")[1]
                        for word in strings[i].split(" ")[2:]:
                            if word.startswith("OS="):
                                break
                            seqDescDB[seqID] += " " + word
                    else:
                        seqDescDB[seqID] = ""
                    i += 1
                    while ((i < len(strings)) and
                            ((not len(strings[i])) or strings[i][0] != '>')):
                        seqlen = 0
                        for ch in strings[i]:
                            if ch in "ABCDEFGHIKLMNOPQRSTUVWXYZ":
                                seqlen += 1
                        seqDB[seqID] += seqlen
                        i += 1
                    i -= 1
            i += 1

        for i in seqDB:
            if seqDB[i] == 0:
                input("!!!ERROR!!!\nLength of sequence with id " + i + " = 0")
                raise()

        unusedOp = "".join([ch for ch in input_unused if ch in "!<>="])
        unused = input_unused[len(unusedOp):].strip()
        if isFloat(unused):
            unused = float(unused)
        else:
            if len(unusedOp):
                tempFile = open(unused.strip())
                strings = tempFile.read().replace(' ', '\t').split('\n')
                tempFile.close()
                unused = {}
                for string in strings:
                    unused[string.split('\t')[0]] = string.split('\t')[-1]
        exit(0)
        contribOp, contrib = getOpData(input_contribution)
        confOp, conf = getOpData(input_confidence)
        minGroupsWithId = int(input_minGroupsWithId)
        maxGroupAbsence = input_maxGroupAbsence
        return (idFileName, IDs, IDExclusionList, seqDB, seqDescDB, unusedOp,
                unused, contribOp, contrib, confOp, conf, minGroupsWithId,
                maxGroupAbsence)
    else:
        input("!!!ERROR!!!\nNo database file!")
        raise()


def readInputFiles():
    """
    Считывание и занесение в оперативную память информации
    из файлов из папки INPUTPATH
    """
    fileNames = []
    for fileName in listdir(path=INPUTPATH):
        fileNames.append(fileName)
    db = {}
    # Сортировка массива с именами файлов
    newArr = [float(val.split("_")[0]) for val in fileNames]
    for i in range(0, len(newArr)):
        for j in range(i, len(newArr)):
            if newArr[i] > newArr[j]:
                newArr[i] += newArr[j]
                newArr[j] = newArr[i] - newArr[j]
                newArr[i] -= newArr[j]
                newArr[i] = round(newArr[i], 5)
                newArr[j] = round(newArr[j], 5)

    shortFileNames = [str(val) for val in newArr]

    for name in fileNames:
        tempFile = open(INPUTPATH + "/" + name)
        tempArray = []
        for string in tempFile.read().split("\n")[1:]:
            tempArray.append(string.split("\t"))
        if tempArray[-1] == ['']:
            tempArray.pop()
        db[name.split("_")[0]] = tempArray
        del(tempArray)
        tempFile.close()

    rawOutput = {}
    for fileName in shortFileNames:
        if fileName != "":
            rawOutput[fileName] = {}

    return fileNames, db, shortFileNames, rawOutput


# Структурирование данных из файлов, загруженных
# в оперативную память. Первичный параметр - имя файла
def getInfoFromFiles(idFileName, IDs, IDExclusionList, seqDB, db, rawOutput,
                     confOp, conf, contribOp, contrib, unusedOp, unused):
    # Сравнивает arg1 и arg2 с помощью операции, указанной в op
    def compare(op, arg1, arg2):
        if op == "" or arg2 == "":
            return True
        if op == "=":
            return arg1 == arg2
        return(eval(str(arg1) + op + str(arg2)))

    if isFloat(unused):
        curUnused = unused
    for fileName in db:
        if not isFloat(unused):
            try:
                curUnused = unused[fileName]
            except KeyError:
                curUnused = ""
        rawOutput[fileName]["Sc_file_summ"] = 0
        rawOutput[fileName]["Psignal_file_summ"] = 0

        for line in db[fileName]:
            if line == "":
                continue
            lineId = line[positions["Id"]].split("; ")[0]
            if(compare(confOp, line[positions["Confidence"]], conf) and
               compare(contribOp, line[positions["Contribution"]], contrib) and
               compare(unusedOp, line[positions["Unused"]], curUnused)):
                # Если ID находится в белом списке, или если файл с
                # белым списком не был указан
                if (lineId in IDs) or idFileName == "":
                    # Если ID не находится в чёрном списке
                    if not (lineId in IDExclusionList):
                        # Если такого ID ещё не было, добавляем его в базу
                        if not (lineId in IDs):
                            IDs.append(lineId)
                        # Если такого ID ещё не было, добавляем его в базу
                        if not (lineId in rawOutput[fileName]):
                            rawOutput[fileName][lineId] = ID()
                        # Заносим значения из файла в объект
                        rawOutput[fileName][lineId]["counts"] += 1
                        rawOutput[fileName][lineId]["seq_length_summ"] += (
                            len(line[positions["Sequence"]]))
                        rawOutput[fileName][lineId]["Sc_summ"] += (
                            float(line[positions["Sc"]]))
                        rawOutput[fileName][lineId]["Psignal_summ"] += (
                            float(line[positions["Precursor Signal"]]))
                        rawOutput[fileName][lineId]["Unused"] = (
                            float(line[positions["Unused"]]))
                if not (lineId in IDExclusionList):
                    # Тест на присутствие id в Database file
                    if not (lineId in seqDB):
                        print("ERROR:Incomplete database!")
                        print("File: " + fileName + "\nId: " + lineId)
                        raise()

                    rawOutput[fileName]["Sc_file_summ"] += (
                        float(line[positions["Sc"]]))
                    rawOutput[fileName]["Psignal_file_summ"] += (
                        float(line[positions["Precursor Signal"]]))

    for name in rawOutput:
        rawOutput[name]["sc_norm_summ"] = 0
        rawOutput[name]["psignal_norm_summ"] = 0
    #  Подсчёт отношений Psignal_summ и Sc_summ к длине последовательности
    # из Database file для каждого ID
    for id in IDs:
        for name in rawOutput:
            if id in rawOutput[name]:
                try:
                    rawOutput[name][id]["Sc_norm"] = (
                        rawOutput[name][id]["Sc_summ"] / seqDB[id])
                    rawOutput[name][id]["Psignal_norm"] = (
                        rawOutput[name][id]["Psignal_summ"] / seqDB[id])
                except Exception as e:
                    print(name, id, sep="\n")
                    print(seqDB[id], end="\n\n")
                    raise e
                rawOutput[name]["sc_norm_summ"] += (
                    rawOutput[name][id]["Sc_norm"])
                rawOutput[name]["psignal_norm_summ"] += (
                    rawOutput[name][id]["Psignal_norm"])
    return IDs, rawOutput


# Конвертация в trueOutput (первичным параметром становится id) и
# удаление из вывода id, не соответствующих условию (группы)
def rawOutputToTrueOutput(IDs, seqDB, rawOutput, shortFileNames,
                          maxGroupAbsence, minGroupsWithId):
    groupCount = 1
    for i in range(1, len(shortFileNames)):
        if(shortFileNames[i].split('.')[0] !=
           shortFileNames[i - 1].split('.')[0]):
            groupCount += 1
        i += 1

    # Находим те элементы, которые отсутствуют слишком много раз
    if (maxGroupAbsence != ""):
        for id in IDs:
            i = 0
            groupsWithIdNum = groupCount
            while (i < len(shortFileNames)):
                j = i
                groupAbsence = 0
                while (j < len(shortFileNames)):
                    if(shortFileNames[j].split(".")[0] !=
                       shortFileNames[i].split(".")[0]):
                        break
                    if not (id in rawOutput[shortFileNames[j]]):
                        groupAbsence += 1
                    j = j + 1
                if groupAbsence > int(maxGroupAbsence):
                    groupsWithIdNum -= 1
                i = j
            if (groupsWithIdNum < minGroupsWithId):
                for name in rawOutput:
                    if id in rawOutput[name]:
                        rawOutput[name].pop(id)

    # Конвертируем в trueOutput
    trueOutput = {}
    for id in IDs:
        for name in rawOutput:
            if id in rawOutput[name]:
                if id not in trueOutput:
                    trueOutput[id] = []
                rawOutput[name][id]["Sc_norm"] = (
                    rawOutput[name][id]["Sc_norm"] /
                    rawOutput[name]["sc_norm_summ"])
                rawOutput[name][id]["Psignal_norm"] = (
                    rawOutput[name][id]["Psignal_norm"] /
                    rawOutput[name]["psignal_norm_summ"])
                rawOutput[name][id]["SP_2"] = ((
                    rawOutput[name][id]["Sc_norm"] +
                    rawOutput[name][id]["Psignal_norm"]) / 2)
                rawOutput[name][id]["seq_length"] = seqDB[id]
                trueOutput[id].append(
                    {"fileName": name, "params": rawOutput[name][id]})
    return trueOutput, rawOutput


def outputResultsToFiles(rawOutput, trueOutput, seqDB,
                         seqDescDB, shortFileNames):
    # Создание файлов по параметрам (1 файл - 1 параметр)
    outputFiles = {}
    if not exists("Output"):
        mkdir("Output")
    for param in params:
        outputFiles[param] = (param + ".txt").replace(" ", "_")
    for param in params:
        outputFile = open(OUTPUTPATH + outputFiles[param], "w")
        lines = [[id] for id in trueOutput]
        if param == "Sc_summ":
            lines.insert(0, ["ID", "Length", *shortFileNames])
        else:
            lines.insert(0, ["ID", *shortFileNames])
        for line in lines[1:]:
            if param == "Sc_summ":
                line.append(str(seqDB[line[0]]))
            for fileName in shortFileNames:
                if line[0] in rawOutput[fileName]:
                    line.append(str(rawOutput[fileName][line[0]][param]))
                else:
                    line.append("")

        outputFile.write("\n".join(["\t".join(line) for line in lines]))
        outputFile.close()
    # Создание файла с «общим выводом» (со всеми параметрами)
    with open(OUTPUTPATH + "output.txt", "w") as outputFile:
        outputFile.write("ID\tFileName\t" + "\t".join(outParams) + "\n")
        for key in trueOutput:
            for fileName in shortFileNames:
                for obj in trueOutput[key]:
                    if obj["fileName"] != fileName:
                        continue
                    outputFile.write(
                        ("{}\t{}" + "\t{}" * len(outParams) + "\n").format(
                            key, obj["fileName"],
                            *[obj["params"][val] for val in outParams]
                        )
                    )
                    break
        outputFile.close()

    with open(OUTPUTPATH + "description.txt", 'w') as descFile:
        descFile.write("ID\tDescription\n")
        for key in trueOutput:
            descFile.write(str(key) + "\t" + str(seqDescDB[key]) + "\n")
        descFile.close()


def main():
    if len(sys.argv) == 9:
        inputs = sys.argv[1:]
    else:
        print(sys.argv)
        inputs = [
                    input("Id file name: "),
                    input("ID exclusion file name: "),
                    input("Database file name: "),
                    input("Unused: "),
                    input("Contribution: "),
                    input("Confidence: "),
                    input("Min groups with ID: "),
                    input("Max missing values per group: ")]
    # Получение основных значений из пользовательского ввода
    (idFileName, IDs, IDExclusionList, seqDB, seqDescDB, unusedOp, unused,
     contribOp, contrib, confOp, conf, minGroupsWithId, maxGroupAbsence) = \
        inputParams(inputs[0],
                    inputs[1],
                    inputs[2],
                    inputs[3],
                    inputs[4],
                    inputs[5],
                    inputs[6],
                    inputs[7])

    # Считывание и занесение в оперативную память информации из папки INPUTPATH
    fileNames, db, shortFileNames, rawOutput = readInputFiles()
    # Структурирование данных из файлов, загруженных в оперативную память.
    # Первиный параметр - имя файла
    IDs, rawOutput = getInfoFromFiles(idFileName, IDs, IDExclusionList,
                                      seqDB, db, rawOutput, confOp, conf,
                                      contribOp, contrib, unusedOp, unused)
    trueOutput, rawOutput = rawOutputToTrueOutput(IDs, seqDB, rawOutput,
                                                  shortFileNames,
                                                  maxGroupAbsence,
                                                  minGroupsWithId)
    outputResultsToFiles(rawOutput, trueOutput, seqDB,
                         seqDescDB, shortFileNames)


try:
    main()
except Exception as e:
    input(str(e))
    raise(e)
