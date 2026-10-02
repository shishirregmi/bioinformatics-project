"""Extract Supplement Table 1 from a locally archived publisher DOCX."""
import argparse
import csv
import hashlib
import io
from pathlib import Path
import xml.etree.ElementTree as ET
import zipfile

SOURCE_SHA256 = "3f6d3dd94160a18c418f2eb6911addcd88b7a1dd5664273801a20290f2c7186b"
NS = {"w": "http://schemas.openxmlformats.org/wordprocessingml/2006/main"}


def extract(path):
    content = Path(path).read_bytes()
    if hashlib.sha256(content).hexdigest() != SOURCE_SHA256:
        raise ValueError("Supplement checksum differs from the verified publisher file.")
    document = ET.fromstring(zipfile.ZipFile(io.BytesIO(content)).read("word/document.xml"))
    table = document.find(".//w:tbl", NS)
    rows = [["".join(t.text or "" for t in cell.findall(".//w:t", NS))
             for cell in row.findall("w:tc", NS)] for row in table.findall("w:tr", NS)]
    if rows[0] != ["Gene symbol", "Entrez ID", "Coefficient"] or len(rows) != 35:
        raise ValueError("Unexpected RSS table structure.")
    return [{"gene": gene, "entrez_id": entrez, "weight": float(weight)}
            for gene, entrez, weight in rows[1:]]


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("supplement")
    parser.add_argument("output")
    args = parser.parse_args()
    genes = extract(args.supplement)
    with open(args.output, "w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=["gene", "entrez_id", "weight"])
        writer.writeheader()
        writer.writerows(genes)
