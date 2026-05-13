import csv
import shutil
import subprocess
from pathlib import Path


SPECIES = "Thermomyces_lanuginosus" #who are we filtering?

LOOKUP_CSV       = Path(f"/home/jkniblo/Ameoba/isolate_ortholog_sequences/{SPECIES}/lookup_all.csv")
INCENDI_STRUCTS  = Path(f"/home/jkniblo/Ameoba/AF_predicted/incendi_strucs/")
ORTHOLOG_STRUCTS = Path(f"/home/jkniblo/Ameoba/AF_predicted/{SPECIES}_strucs/")
FOLDSEEK_DIR     = Path(f"/home/jkniblo/Ameoba/foldseek/{SPECIES}")
THREADS          = 8

##Get ortholog pairs
pairs      = []
query_ids  = set()
target_ids = set()

with open(LOOKUP_CSV) as f:
    reader = csv.DictReader(f)
    for row in reader:
        q = row["incendi_protein"].strip()
        t = row["ortholog_accession"].strip()
        pairs.append((q, t))
        query_ids.add(q)
        target_ids.add(t)

pair_set = set(pairs)



##Copy structures
#TODO: optimization - make it so it looks in different directories
incendi_db_dir  = FOLDSEEK_DIR / "incendi_structs"
ortholog_db_dir = FOLDSEEK_DIR / "ortholog_structs"
incendi_db_dir.mkdir(parents=True, exist_ok=True)
ortholog_db_dir.mkdir(parents=True, exist_ok=True)


def copy_structures(ids, src_dir, dst_dir):
    copied, missing = 0, []
    for protein_id in ids:
        hits = list(src_dir.glob(f"{protein_id}.pdb"))
        if not hits:
            missing.append(protein_id)
            continue
        for f in hits:
            shutil.copy(f, dst_dir / f.name)
            copied += 1
    return copied, missing


copied_q, missing_q = copy_structures(query_ids,  INCENDI_STRUCTS,  incendi_db_dir)
copied_t, missing_t = copy_structures(target_ids, ORTHOLOG_STRUCTS, ortholog_db_dir)

if missing_q: print("  Missing incendi:",   missing_q[:10])
if missing_t: print("  Missing orthologs:", missing_t[:10])


##Run FoldSeek 
db_q        = FOLDSEEK_DIR / "incendi_db"
db_t        = FOLDSEEK_DIR / "ortholog_db"
raw_results = FOLDSEEK_DIR / "raw_results.tsv"
tmp_dir     = FOLDSEEK_DIR / "tmp"

# clean up tmp to avoid cached results from previous runs
if tmp_dir.exists():
    shutil.rmtree(tmp_dir)

subprocess.run(["foldseek", "createdb", str(incendi_db_dir),  str(db_q)], check=True)
subprocess.run(["foldseek", "createdb", str(ortholog_db_dir), str(db_t)], check=True)

subprocess.run([
    "foldseek", "easy-search",
    str(db_q), str(db_t),
    str(raw_results), str(tmp_dir),
    "--format-output", "query,target,alntmscore,pident,alnlen,evalue,prob",
    "--threads", str(THREADS),
    "--exhaustive-search", "1",
    "--exact-tmscore", "1",
], check=True)


# Keep potential ortholog pairs 
## FoldSeek does an all v all and we don't care about the other comparisons 
out_file = FOLDSEEK_DIR / "filtered_pairs.tsv"
header   = ["incendi_protein", "ortholog_accession",
            "alntmscore", "pident", "alnlen", "evalue", "prob"]

pair_results = {}   # (q, t) → best result row by tmscore

with open(raw_results) as f:
    for line in f:
        parts = line.strip().split()
        if len(parts) < 7:
            continue

        q = parts[0]
        t = parts[1]

        if (q, t) not in pair_set:
            continue

        tmscore = float(parts[2])
        key = (q, t)
        if key not in pair_results or tmscore > pair_results[key][2]:
            pair_results[key] = [q, t] + parts[2:]

with open(out_file, "w", newline="") as f:
    writer = csv.writer(f, delimiter="\t")
    writer.writerow(header)
    for row in sorted(pair_results.values(), key=lambda x: -float(x[2])):
        writer.writerow(row)

# Write approved orthologs (TM-score >= 0.5)
approved_file = FOLDSEEK_DIR / "approved_orthologs.tsv"
approved_rows = [row for row in pair_results.values() if float(row[2]) >= 0.5]
with open(approved_file, "w", newline="") as f:
    writer = csv.writer(f, delimiter="\t")
    writer.writerow(header)
    for row in sorted(approved_rows, key=lambda x: -float(x[2])):
        writer.writerow(row)

##Return how many orthologs have score >0.5
tmscores = [float(r[2]) for r in pair_results.values()]
if tmscores:
    high   = sum(1 for s in tmscores if s >= 0.5)
    low    = sum(1 for s in tmscores if s < 0.5)
    print(f"\nTM-score distribution:")
    print(f"  High   (>=0.5)   : {high}")
    print(f"  Low    (<0.5)    : {low}  (should probably drop these)")
    print(f"\nApproved orthologs written to: {approved_file}")