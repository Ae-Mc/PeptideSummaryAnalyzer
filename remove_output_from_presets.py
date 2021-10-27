import os
from shutil import rmtree
from typing import List


to_be_deleted: List[str] = []
for preset in os.listdir("Presets"):
    preset = os.path.join("Presets", preset)
    for element in os.listdir(preset):
        if element == "Output":
            to_be_deleted.append(os.path.join(preset, element))
            print(to_be_deleted[-1] + " will be deleted")

delete = input("Delete (Y/N)? ").upper() == "Y"
if delete:
    for file_path in to_be_deleted:
        rmtree(file_path)
