#!/bin/env python3
import zipfile
from os import listdir
from os.path import join

with zipfile.ZipFile(
    "PeptideSummaryAnalyzer.zip", "w", zipfile.ZIP_STORED
) as archive:
    archive.write("peptide_summary_analyzer.py")
    for filename in listdir("classes"):
        if not filename.startswith("__") and not filename.startswith("."):
            archive.write(join("classes", filename))
