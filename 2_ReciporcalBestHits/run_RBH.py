import re
import torch
import numpy as np
import faiss
import csv
from pathlib import Path
from collections import defaultdict


BASE_DIR    = Path("/home/jkniblo/Ameoba/embeddings") #Path to embeddings 
QUERY_500   = "incendi_500" #proteins of interest
QUERY_ALL   = "incendi_all" #full proteome of query species 
RUN_NAME    = "incendi_500"
OUT_DIR     = Path("reciprocal_best_hits_out")

LAYERS      = [20, 25] #specify layers to retain from ESM -- see https://genome.cshlp.org/content/early/2024/09/30/gr279127124
K           = 2                      
THRESHOLD   = 0.85  #cosine similarity threshold for candidate hits; adjust based on desired sensitivity/specificity and layer choice                
RUN_TARGETS = ["Thermomyces_lanuginosus"] ##Can also put none, will run on everything in directory.


def extract_protein_id(pt_file):
    full_path = str(pt_file)
    matches = re.findall(r'\|([A-Z0-9]{6,12})\|', full_path)
    if matches:
        return matches[0]
    match = re.search(r'GN=(\S+)', full_path)
    if match:
        return match.group(1)
    return pt_file.stem


def load_proteome_multilayer(species_dir, layers):
    ids  = []
    vecs = defaultdict(list)

    # account for some embeddings being named weird
    pt_files = sorted(Path(species_dir).glob("*.pt"))
    if not pt_files:
        pt_files = sorted(Path(species_dir).rglob("*.pt"))
    if not pt_files:
        pt_files = sorted([f for f in Path(species_dir).iterdir()
                           if f.is_file() and f.suffix == ""])
    if not pt_files:
        raise FileNotFoundError(f"No .pt files found in {species_dir}")

    # load embeddings, hold each layer
    for pt_file in pt_files:
        data = torch.load(pt_file, map_location="cpu")
        if "mean_representations" not in data:
            raise KeyError(f"'mean_representations' missing in {pt_file}. "
                           f"Keys: {list(data.keys())}")
        for layer in layers:
            if layer not in data["mean_representations"]:
                raise KeyError(f"Layer {layer} not found in {pt_file}. "
                               f"Available: {list(data['mean_representations'].keys())}")
            vecs[layer].append(data["mean_representations"][layer].float().numpy())
        ids.append(extract_protein_id(pt_file))

    # return a matrix for each layer/protein
    mats = {layer: np.stack(vecs[layer], axis=0).astype("float32")
            for layer in layers}

    return ids, mats


def normalize(embeddings_vectors):
    # make all vectors same length but keep direction
    # dot product of normalized vectors = cosine similarity
    vectors = embeddings_vectors.copy()
    faiss.normalize_L2(vectors)
    return vectors


def top_k_candidates_single_layer(q500_ids, q500_norm,
                                   qall_ids, qall_norm,
                                   t_ids,    t_norm,
                                   layer):
    dim      = q500_norm.shape[1]
    n_target = len(t_ids)

    # forward: query → ALL target proteins
    idx_target = faiss.IndexFlatIP(dim)
    idx_target.add(t_norm)
    D_all, I_all = idx_target.search(q500_norm, k=n_target)

    # reverse: target → full query proteome, top-K
    idx_qall = faiss.IndexFlatIP(dim)
    idx_qall.add(qall_norm)
    D_rev, I_rev = idx_qall.search(t_norm, k=K)   # keep D_rev for reverse gap

    # for each target protein, which query IDs were in its top-K?
    reverse_topk_ids = {}
    for j in range(len(t_ids)):
        top_query_ids = set()
        for idx in I_rev[j]:
            query_id = qall_ids[idx]
            top_query_ids.add(query_id)
        reverse_topk_ids[j] = top_query_ids

    candidates = []
    for index_query, q_id in enumerate(q500_ids):
        hits = []
        for rank in range(K):
            index_target = int(I_all[index_query, rank])
            score        = float(D_all[index_query, rank])
            if score < THRESHOLD:
                break
            if q_id not in reverse_topk_ids[index_target]:   # reciprocity check
                continue
            hits.append((rank + 1, index_target, score))

        if not hits:
            continue

        # forward gap: how unambiguous was the query's best target hit?
        fwd_gap = round(float(D_all[index_query, 0]) - float(D_all[index_query, 1]), 6) \
                  if n_target >= 2 else 0.0

        for rank, index_target, score in hits:
            # reverse gap: how unambiguous was the target's best query hit?
            rev_gap = round(float(D_rev[index_target, 0]) - float(D_rev[index_target, 1]), 6) \
                      if len(qall_ids) >= 2 else 0.0

            # combined gap: both gaps large = high confidence in both directions
            combined_gap = round((fwd_gap + rev_gap) / 2, 6)

            candidates.append({
                "query_id":    q_id,
                "target_id":   t_ids[index_target],
                "layer":       layer,
                "cosine_sim":  round(score, 6),
                "fwd_rank":    rank,
                "fwd_gap":     fwd_gap,
                "rev_gap":     rev_gap,
                "combined_gap": combined_gap,
            })

    return candidates


def merge_layer_candidates(all_layer_candidates):
    pair_hits = defaultdict(list)
    for c in all_layer_candidates:
        pair_hits[(c["query_id"], c["target_id"])].append(c)

    merged = []
    for (q_id, t_id), hits in pair_hits.items():
        best         = max(hits, key=lambda x: x["cosine_sim"])
        layers_hit   = sorted(set(h["layer"] for h in hits))
        hit_by_layer = {h["layer"]: h["cosine_sim"] for h in hits}

        merged.append({
            "query_id":      q_id,
            "target_id":     t_id,
            "best_layer":    best["layer"],
            "best_cosine":   best["cosine_sim"],
            "fwd_rank":      best["fwd_rank"],
            "fwd_gap":       best["fwd_gap"],
            "rev_gap":       best["rev_gap"],
            "combined_gap":  best["combined_gap"],
            "layer_support": len(layers_hit),
            "layers_hit":    ",".join(str(l) for l in layers_hit),
            **{f"cos_L{l}": round(hit_by_layer.get(l, 0.0), 6) for l in LAYERS},
        })

    return merged


# ─── MAIN ─────────────────────────────────────────────────────────────────────

OUT_DIR.mkdir(parents=True, exist_ok=True)

incendi_dirs = {QUERY_500, QUERY_ALL}
targets = RUN_TARGETS or sorted([
    d.name for d in BASE_DIR.iterdir()
    if d.is_dir() and d.name not in incendi_dirs
])

q500_ids, q500_mats = load_proteome_multilayer(BASE_DIR / QUERY_500, LAYERS)
qall_ids, qall_mats = load_proteome_multilayer(BASE_DIR / QUERY_ALL, LAYERS)

q500_norms = {l: normalize(q500_mats[l]) for l in LAYERS}
qall_norms = {l: normalize(qall_mats[l]) for l in LAYERS}

summary_rows = []

for species in targets:
    try:
        t_ids, t_mats = load_proteome_multilayer(BASE_DIR / species, LAYERS)
    except FileNotFoundError as e:
        print(f"  SKIPPED: {e}")
        continue

    t_norms = {l: normalize(t_mats[l]) for l in LAYERS}

    all_layer_candidates = []
    for layer in LAYERS:
        layer_cands = top_k_candidates_single_layer(
            q500_ids, q500_norms[layer],
            qall_ids, qall_norms[layer],
            t_ids,    t_norms[layer],
            layer,
        )
        all_layer_candidates.extend(layer_cands)

    candidates = merge_layer_candidates(all_layer_candidates)
    queries_with_hit = len(set(c["query_id"] for c in candidates))

    layer_dist = defaultdict(int)
    for c in candidates:
        layer_dist[c["best_layer"]] += 1

    multi = [c for c in candidates if c["layer_support"] > 1]
    candidates_sorted = sorted(candidates,
                               key=lambda x: (-x["layer_support"],
                                              -x["best_cosine"]))

    fieldnames = (["query_id", "target_id",
                   "best_layer", "best_cosine", "fwd_rank",
                   "fwd_gap", "rev_gap", "combined_gap",
                   "layer_support", "layers_hit"]
                  + [f"cos_L{l}" for l in LAYERS])

    out_file = OUT_DIR / f"{RUN_NAME}_vs_{species}_candidates.tsv"
    with open(out_file, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=fieldnames, delimiter="\t")
        w.writeheader()
        w.writerows(candidates_sorted)

    summary_rows.append({
        "species":            species,
        "n_query":            len(q500_ids),
        "n_query_all":        len(qall_ids),
        "n_target":           len(t_ids),
        "n_candidates":       len(candidates),
        "queries_with_hit":   queries_with_hit,
        "pct_query_hit":      round(queries_with_hit / len(q500_ids) * 100, 1),
        "mean_cosine":        round(np.mean([c["best_cosine"] for c in candidates]), 4)
                              if candidates else 0.0,
        "mean_combined_gap":  round(np.mean([c["combined_gap"] for c in candidates]), 4)
                              if candidates else 0.0,
        "mean_layer_support": round(np.mean([c["layer_support"] for c in candidates]), 2)
                              if candidates else 0.0,
        "n_multilayer_hits":  len(multi),
        "pct_multilayer":     round(len(multi) / max(len(candidates), 1) * 100, 1),
        **{f"n_best_L{l}": layer_dist[l] for l in LAYERS},
    })

summary_file = OUT_DIR / "candidates_summary.tsv"
fieldnames_sum = (["species", "n_query", "n_query_all", "n_target",
                   "n_candidates", "queries_with_hit", "pct_query_hit",
                   "mean_cosine", "mean_combined_gap", "mean_layer_support",
                   "n_multilayer_hits", "pct_multilayer"]
                  + [f"n_best_L{l}" for l in LAYERS])

with open(summary_file, "w", newline="") as f:
    w = csv.DictWriter(f, fieldnames=fieldnames_sum, delimiter="\t")
    w.writeheader()
    w.writerows(summary_rows)

