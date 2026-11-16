"""Download the public inputs of the reconstruction into gem/raw_data/.

  MOESM3_cbs51388.xml       template: the CBS 513.88 version of iJB1325 (Brandl et al. 2018, Additional file 3)
  fven_proteins.faa         F. venenatum A3/5 proteome, UniProt UP000245910 (assembly GCA_900007375.1), main and
                            additional FASTA files (13,945 proteins)
  aniger_cbs513.faa         A. niger CBS 513.88 proteome, UniProt UP000006706, with the An locus tag (the template's
                            gene id) as FASTA id; entries without one keep their UniProt accession
  orthofinder_input/*.faa   the 10 proteomes of the OrthoFinder run: the two above, with ids reduced to the
                            accession or locus tag, and eight NCBI RefSeq annotations

The template file is checked by its MD5. UniProt and RefSeq releases change header text and record order over
time, so each proteome is checked by an MD5 over its sorted "id<TAB>sequence" lines. A mismatch means the database
has changed since the version used for iJL1143; the orthology tables in inputs/ were built from that version.

Usage (from gem/):
    python download_data.py                 # everything
    python download_data.py --only template # the template model only (enough for model building and Fig. 1)
"""

from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import io
import re
import sys
import textwrap
import urllib.request
from pathlib import Path

HERE = Path(__file__).resolve().parent
RAW = HERE / "raw_data"
USER_AGENT = "Mozilla/5.0 (compatible; fven-gem data download)"

TEMPLATE = {
    "url": "https://static-content.springer.com/esm/art%3A10.1186%2Fs40694-018-0060-7/MediaObjects/40694_2018_60_MOESM3_ESM.xml",
    "file": "MOESM3_cbs51388.xml",
    "md5": "fc6bccfe766d9c74ee572684a36e5be1",
}
UNIPROT_FTP = "https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/reference_proteomes/Eukaryota"
FVEN_URLS = [f"{UNIPROT_FTP}/UP000245910/UP000245910_56646.fasta.gz",
             f"{UNIPROT_FTP}/UP000245910/UP000245910_56646_additional.fasta.gz"]
ANIGER_TSV = ("https://rest.uniprot.org/uniprotkb/stream?query=proteome:UP000006706"
              "&fields=accession,gene_orf,sequence&format=tsv")
NCBI = "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF"
REFSEQ = {   # OrthoFinder species name -> RefSeq assembly (path below genomes/all/GCF)
    "Afumigatus": "000/002/655/GCF_000002655.1_ASM265v1",                                   # A. fumigatus Af293
    "Anidulans": "000/011/425/GCF_000011425.1_ASM1142v1",                                   # A. nidulans FGSC A4
    "Aoryzae": "000/184/455/GCF_000184455.2_ASM18445v3",                                    # A. oryzae RIB40
    "Aterreus": "000/149/615/GCF_000149615.1_ASM14961v1",                                   # A. terreus NIH2624
    "Ffujikuroi": "900/079/805/GCF_900079805.1_Fusarium_fujikuroi_IMI58289_V2",            # F. fujikuroi IMI 58289
    "Fgraminearum": "000/240/135/GCF_000240135.3_ASM24013v3",                               # F. graminearum PH-1
    "Fverticillioides": "000/149/555/GCF_000149555.1_ASM14955v1",                           # F. verticillioides 7600
    "Treesei": "000/167/675/GCF_000167675.1_v2.0",                                          # T. reesei QM6a
}
# MD5 over the sorted "id<TAB>sequence" lines, and the number of proteins, of the files used for iJL1143
SEQUENCE_MD5 = {
    "fven_proteins.faa": ("21bcbfd629541dc76d071497095b5f5e", 13945),
    "aniger_cbs513.faa": ("a1582a167df11d36e575719832729ba1", 14085),
    "Afumigatus": ("277d829d7cbfce8e44b2ad46e3eb6c32", 9823),
    "Anidulans": ("8c4643a450aec9d7e195e2be03718c54", 10453),
    "Aoryzae": ("9752f5e0810656d6aa6f6c488b67541c", 12074),
    "Aterreus": ("47f9083669125bc7659ce18d237581b0", 10401),
    "Ffujikuroi": ("9baec53e11050e5c58505ddd90545052", 14810),
    "Fgraminearum": ("4d73fc1fe18eec26cccd650a74fffe48", 13312),
    "Fverticillioides": ("70ca91d2ed6a2eecd860e33af1b2eaf8", 20574),
    "Treesei": ("368254f6082c98d535268fdbbaa66aab", 9111),
}
LOCUS = re.compile(r"\bAn\d\dg\d+\b")


def fetch(url: str) -> bytes:
    print(f"  {url}")
    request = urllib.request.Request(url, headers={"User-Agent": USER_AGENT})
    with urllib.request.urlopen(request, timeout=600) as response:
        data = response.read()
    return gzip.decompress(data) if url.endswith(".gz") else data


def records(text: str):
    """(id, sequence) per FASTA record; the id is the UniProt accession for UniProt headers, else the first token."""
    key, seq = None, []
    for line in text.splitlines():
        if line.startswith(">"):
            if key is not None:
                yield key, "".join(seq)
            head = line[1:].strip()
            m = re.match(r"(?:sp|tr)\|([^|]+)\|", head)
            key, seq = (m.group(1) if m else head.split()[0]), []
        elif line.strip():
            seq.append(line.strip())
    if key is not None:
        yield key, "".join(seq)


def fasta(pairs) -> str:
    return "".join(f">{k}\n" + "\n".join(textwrap.wrap(s, 60)) + "\n" for k, s in pairs)


def check(name: str, text: str) -> None:
    lines = sorted(f"{k}\t{s}" for k, s in records(text))
    digest = hashlib.md5("\n".join(lines).encode()).hexdigest()
    expected, n = SEQUENCE_MD5[name]
    if digest == expected:
        print(f"  {name}: {len(lines):,} proteins, sequences identical to the version used for iJL1143")
    else:
        print(f"  WARNING {name}: {len(lines):,} proteins (used: {n:,}); the sequences differ from the version used "
              f"for iJL1143, so a rerun of the orthology step may not reproduce inputs/orthologs_cbs513.tsv")


def template(force: bool) -> None:
    path = RAW / TEMPLATE["file"]
    if force or not path.exists():
        path.write_bytes(fetch(TEMPLATE["url"]))
    digest = hashlib.md5(path.read_bytes()).hexdigest()
    if digest != TEMPLATE["md5"]:
        sys.exit(f"ERROR: {path.name} has MD5 {digest}, expected {TEMPLATE['md5']}")
    print(f"  {path.name}: MD5 matches")


def proteomes(force: bool) -> None:
    of_dir = RAW / "orthofinder_input"
    of_dir.mkdir(parents=True, exist_ok=True)

    path = RAW / "fven_proteins.faa"
    if force or not path.exists():
        path.write_bytes(b"".join(fetch(u) for u in FVEN_URLS))
    text = path.read_text()
    check("fven_proteins.faa", text)
    (of_dir / "Fvenenatum.faa").write_text(fasta(records(text)))

    path = RAW / "aniger_cbs513.faa"
    if force or not path.exists():
        rows = csv.DictReader(io.StringIO(fetch(ANIGER_TSV).decode()), delimiter="\t")
        pairs = []
        for r in rows:
            m = LOCUS.search(r["Gene Names (ORF)"] or "")
            pairs.append((m.group(0) if m else r["Entry"], r["Sequence"]))
        path.write_text(fasta(pairs))
    text = path.read_text()
    check("aniger_cbs513.faa", text)
    (of_dir / "Aniger_CBS513.faa").write_text(text)

    for species, assembly in REFSEQ.items():
        path = of_dir / f"{species}.faa"
        if force or not path.exists():
            name = assembly.rsplit("/", 1)[1]
            path.write_bytes(fetch(f"{NCBI}/{assembly}/{name}_protein.faa.gz"))
        check(species, path.read_text())


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__.split("\n", 1)[0])
    parser.add_argument("--force", action="store_true", help="download again even if a file exists")
    parser.add_argument("--only", choices=("template", "proteomes"))
    args = parser.parse_args()
    RAW.mkdir(parents=True, exist_ok=True)
    if args.only in (None, "template"):
        print("Template model")
        template(args.force)
    if args.only in (None, "proteomes"):
        print("Proteomes")
        proteomes(args.force)


if __name__ == "__main__":
    main()
