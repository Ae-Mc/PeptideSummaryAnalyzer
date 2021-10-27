import zipfile
import os


def archive_folder(archive: zipfile.ZipFile, folder_path: str) -> None:
    for element in os.listdir(folder_path):
        element_path = os.path.join(folder_path, element)
        if not (
            element.startswith("__")
            or element.startswith(".")
        ) or element == "__init__.py":
            if os.path.isfile(element_path):
                archive.write(element_path)
            elif os.path.isdir(element_path):
                archive_folder(archive, element_path)


def main():
    with zipfile.ZipFile(
        "PeptideSummaryAnalyzer.zip", "w", zipfile.ZIP_STORED
    ) as archive:
        archive.write("sql.py")
        archive_folder(archive, "classes")


if __name__ == "__main__":
    main()
