import re
import csv
from pathlib import Path
from collections import defaultdict

Species = "Thermomyces_lanuginosus"

INCENDI_500_FASTA  = Path("/home/jkniblo/Ameoba/fastas/incendi_500.faa")
ORTHOLOG_FASTA     = Path(f"/home/jkniblo/Ameoba/fastas/{Species}.fasta")
RBH_TSV            = Path(f"/home/jkniblo/Ameoba/reciprocal_best_hits/incendi_500_vs_{Species}_candidates.tsv")
ORTHOFINDER_TSV    = Path(f"/home/jkniblo/Ameoba/orthofinder/orthfinder_output/incendi_v_{Species}/Orthologues/Orthologues_incendi_all/incendi_all__v__{Species}.tsv")
OUT_DIR            = Path(f"/home/jkniblo/Ameoba/isolate_ortholog_sequences/{Species}")

# Output CSVs
LOOKUP_ALL         = OUT_DIR / "lookup_all.csv"
LOOKUP_RBH_ONLY    = OUT_DIR / "lookup_rbh_only.csv"
LOOKUP_OF_ONLY     = OUT_DIR / "lookup_orthofinder_only.csv"
LOOKUP_BOTH        = OUT_DIR / "lookup_both_sources.csv"

# Column names in the OrthoFinder TSV
OF_COL_INCENDI = "incendi_all"
OF_COL_TARGET  = Species


def parse_fasta(fasta_path):
    seqs = {}
    header = seq_lines = None
    with open(fasta_path) as f:
        for line in f:
            line = line.rstrip()
            if line.startswith(">"):
                if header is not None:
                    acc = extract_uniprot_acc(header)
                    seqs[acc] = (header, "".join(seq_lines))
                header    = line[1:]
                seq_lines = []
            elif header is not None:
                seq_lines.append(line)
    if header is not None:
        acc = extract_uniprot_acc(header)
        seqs[acc] = (header, "".join(seq_lines))
    return seqs


def extract_uniprot_acc(header):
    m = re.match(r'(?:sp|tr)\|([A-Z0-9]+)\|', header)
    if m:
        return m.group(1)
    return header.split()[0]


def parse_incendi_ids(fasta_path):
    ids = set()
    with open(fasta_path) as f:
        for line in f:
            if not line.startswith(">"):
                continue
            ids.add(line[1:].strip().split()[0])
    print(f"incendi_500 proteins    : {len(ids)}")
    return ids


def parse_of_accession(entry):
    parts = str(entry).strip().split("|")
    return parts[1] if len(parts) >= 2 else entry.strip()



def collect_rbh_targets(rbh_tsv, incendi_ids):
    accs    = set()
    mapping = defaultdict(set)
    with open(rbh_tsv) as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            query = row["query_id"]
            if query in incendi_ids:
                target = row["target_id"].strip()
                accs.add(target)
                mapping[query].add(target)
    print(f"RBH targets             : {len(accs)}")
    return accs, mapping


def collect_of_targets(of_tsv, col_incendi, col_target, incendi_ids):
    accs    = set()
    mapping = defaultdict(set)
    with open(of_tsv) as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            incendi_raw = str(row.get(col_incendi, "")).strip()
            target_raw  = str(row.get(col_target,  "")).strip()
            if incendi_raw == "nan" or target_raw == "nan":
                continue
            row_incendi_ids = {e.strip() for e in incendi_raw.split(",")}
            matched = row_incendi_ids & incendi_ids
            if not matched:
                continue
            target_accs = set()
            for entry in target_raw.split(","):
                acc = parse_of_accession(entry)
                if acc:
                    target_accs.add(acc)
                    accs.add(acc)
            for inc_id in matched:
                mapping[inc_id] |= target_accs
    print(f"OrthoFinder targets     : {len(accs)}")
    return accs, mapping

def write_lookup(path, mapping, label):
    total = sum(len(v) for v in mapping.values())
    with open(path, "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(["incendi_protein", "ortholog_accession"])
        for inc_id in sorted(mapping):
            for target in sorted(mapping[inc_id]):
                writer.writerow([inc_id, target])


OUT_DIR.mkdir(parents=True, exist_ok=True)

incendi_ids = parse_incendi_ids(INCENDI_500_FASTA) #get incendi ids

rbh_accs, rbh_mapping = collect_rbh_targets(RBH_TSV, incendi_ids) #get RBH targets and mapping

of_accs, of_mapping = collect_of_targets(ORTHOFINDER_TSV, OF_COL_INCENDI, OF_COL_TARGET, incendi_ids) #get OrthoFinder targets and mapping

all_accs = rbh_accs | of_accs #get all unique potential pairs from both sources

##Cateorize potential pairs into 4 groups: 1. combined (everything), 2. RBH only, 3. OrthoFinder only, 4. Both OF and RBH

# combined: union of everything
combined_mapping = defaultdict(set)
for inc_id, targets in rbh_mapping.items():
    combined_mapping[inc_id] |= targets
for inc_id, targets in of_mapping.items():
    combined_mapping[inc_id] |= targets

# rbh_only: pairs in RBH but NOT in OrthoFinder for that same incendi protein
rbh_only_mapping = defaultdict(set)
for inc_id, targets in rbh_mapping.items():
    exclusive = targets - of_mapping.get(inc_id, set())
    if exclusive:
        rbh_only_mapping[inc_id] = exclusive

# of_only: pairs in OrthoFinder but NOT in RBH for that same incendi protein
of_only_mapping = defaultdict(set)
for inc_id, targets in of_mapping.items():
    exclusive = targets - rbh_mapping.get(inc_id, set())
    if exclusive:
        of_only_mapping[inc_id] = exclusive

# both: pairs found by BOTH tools for the same incendi protein
both_mapping = defaultdict(set)
for inc_id, targets in rbh_mapping.items():
    shared = targets & of_mapping.get(inc_id, set())
    if shared:
        both_mapping[inc_id] = shared

# write lookup CSVs for each category
write_lookup(LOOKUP_ALL,      combined_mapping, "All (combined)")
write_lookup(LOOKUP_RBH_ONLY, rbh_only_mapping, "RBH only")
write_lookup(LOOKUP_OF_ONLY,  of_only_mapping,  "OrthoFinder only")
write_lookup(LOOKUP_BOTH,     both_mapping,     "Both sources")

# write one fasta per ortholog
#this is to run multiple AF predictions simultaneously
ortholog_seqs = parse_fasta(ORTHOLOG_FASTA)
for acc in sorted(all_accs):
    entry = ortholog_seqs.get(acc)
    if entry is None:
        continue
    header, seq = entry
    out_file = OUT_DIR / f"{acc}.fasta"
    with open(out_file, "w") as f:
        f.write(f">{header}\n")
        for i in range(0, len(seq), 60):
            f.write(seq[i:i+60] + "\n")
