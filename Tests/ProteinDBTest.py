import unittest
import shutil
import os
from Classes.ProteinDB import ProteinDB
from Classes.ProteinGroup import ProteinGroup
from Classes.ProteinAccessionWithMaxUnused import ProteinAccessionWithMaxUnused
from PeptideSummaryAnalyzer import ReadSeqDB


class ProteinDBTest(unittest.TestCase):

    proteinDB: ProteinDB
    TESTDIRNAME = "./TEST_TEMP"

    @classmethod
    def setUpClass(cls):
        try:
            os.mkdir(cls.TESTDIRNAME)
        except FileExistsError:
            pass
        files = {
                "1.1_ProteinSummary.txt": [
                    "Unused\tAccession",
                    "4\tA",
                    "0\tB",
                    "2\tC",
                    "1\tD",
                    "0\tE"],
                "1.2_ProteinSummary.txt": [
                    "Unused\tAccession",
                    "2\tA",
                    "4\tB",
                    "0\tC",
                    "1\tD",
                    "2\tE"],
                "1.3_ProteinSummary.txt": [
                    "Unused\tAccession",
                    "2\tB",
                    "1\tG"],
                "1.4_ProteinSummary.txt": [
                    "Unused\tAccession",
                    "2\tF",
                    "0\tG"],
                "EFRA_cont.fasta": [
                    ">A",
                    "VLSLDKOWKSLALSKDLAS",
                    "ASNONOANSONONASONOD",
                    ">B",
                    "AMSMOMAOSMOMDOMAOMA",
                    "PASPMPAMPSMPMPASMPA",
                    "ORKGMKOFLFKGKGLLFLL",
                    ">C",
                    "PMPMPMKJGJJTNMKFGNK",
                    "GNNTJVONONONRONORNO",
                    "RNRINIGNRINGINBIT"
                    ">D",
                    "PMPMPMKJGJJTNMKFGNK",
                    "GNNTJVONONONRONORNO",
                    "RINIGNRINGINBIT"
                    ">E",
                    "PMPMPMKJGJJTNMKFGNK",
                    "GNNTJVONONONRONORNO",
                    "RNARINIGNRINGINBIT",
                    ">F",
                    "AOINODSONAONIONOANO",
                    "PPAONONAOINIDNION",
                    ">G",
                    "ONIAOSnDIONSONODANS",
                    "OASONOANSINDOISOANO"
                ]}

        for filename, lines in files.items():
            with open(os.path.join(cls.TESTDIRNAME, filename), 'w') as f:
                f.write('\n'.join(lines))

        cls.maxDiff = None
        cls.proteinDB = ProteinDB(
                cls.TESTDIRNAME,
                ReadSeqDB(os.path.join(cls.TESTDIRNAME, "EFRA_cont.fasta")),
                None)

    def setUp(self):
        self.proteinDB.ReadDBFromFolder()

    def testGetAccessionsBunch(self):
        self.assertDictEqual(
                self.proteinDB.GetAccessionsBunch(),
                {
                    "A": ProteinAccessionWithMaxUnused("A", 4.0, 2, False),
                    "B": ProteinAccessionWithMaxUnused("B", 4.0, 3, True),
                    "C": ProteinAccessionWithMaxUnused("C", 4.0, 2, False),
                    "D": ProteinAccessionWithMaxUnused("D", 1.0, 2, True),
                    "E": ProteinAccessionWithMaxUnused("E", 2.0, 2, False),
                    "F": ProteinAccessionWithMaxUnused("F", 2.0, 1, False),
                    "G": ProteinAccessionWithMaxUnused("G", 2.0, 2, False),
                    })

    def testReadDBFromFolder(self):
        self.assertDictEqual(
                self.proteinDB.proteinGroupsPerTable,
                {
                    "1.1": [
                        ProteinGroup(4.0, ["A", "B"]),
                        ProteinGroup(2.0, ["C"]),
                        ProteinGroup(1.0, ["D", "E"])],
                    "1.2": [
                        ProteinGroup(2.0, ["A"]),
                        ProteinGroup(4.0, ["B", "C"]),
                        ProteinGroup(1.0, ["D"]),
                        ProteinGroup(2.0, ["E"])],
                    "1.3": [
                        ProteinGroup(2.0, ["B"]),
                        ProteinGroup(1.0, ["G"])],
                    "1.4": [
                        ProteinGroup(2.0, ["F", "G"])]
                })

    def testCalculateRepresentatives(self):
        self.proteinDB.CalculateRepresentatives()
        self.assertDictEqual(
                self.proteinDB.proteinGroupsPerTable,
                {
                    "1.1": [
                        ProteinGroup(4.0, ["A", "B"], "B"),
                        ProteinGroup(2.0, ["C"],      "C"),
                        ProteinGroup(1.0, ["D", "E"], "E")],
                    "1.2": [
                        ProteinGroup(2.0, ["A"],      "A"),
                        ProteinGroup(4.0, ["B", "C"], "B"),
                        ProteinGroup(1.0, ["D"],      "E"),
                        ProteinGroup(2.0, ["E"],      "E")],
                    "1.3": [
                        ProteinGroup(2.0, ["B"],      "B"),
                        ProteinGroup(1.0, ["G"],      "G")],
                    "1.4": [
                        ProteinGroup(2.0, ["F", "G"], "G")]
                })
        self.assertDictEqual(
                self.proteinDB.difficultCases,
                {
                    "1.1": [ProteinGroup(4.0, ["A", "B"], "B")],
                    "1.2": [ProteinGroup(4.0, ["B", "C"], "B")]})

    def testGetAccessionsReplacementsPerTable(self):
        self.proteinDB.CalculateRepresentatives()
        self.assertDictEqual(
            self.proteinDB.GetAccessionsReplacementsPerTable(),
            {
                "1.1": {
                    "A": "B",
                    "B": "B",
                    "C": "C",
                    "D": "E",
                    "E": "E"},
                "1.2": {
                    "A": "A",
                    "B": "B",
                    "C": "B",
                    "D": "E",
                    "E": "E"},
                "1.3": {
                    "B": "B",
                    "G": "G"},
                "1.4": {
                    "F": "G",
                    "G": "G"}
            })

    @classmethod
    def tearDownClass(cls):
        shutil.rmtree(cls.TESTDIRNAME)
        pass
