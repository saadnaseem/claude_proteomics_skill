"""Unified proteomics-agent pipeline.
Runs Stages 0, 1, 4, 5, 6, 7, 8 end-to-end given a config dict.
Stages 2 (hypothesis), 3 (deep research), 11 (synthesis) are done by the agent in conversation."""
import argparse
import json
import os
import time
import hashlib
import subprocess
from collections import defaultdict
from datetime import datetime, timezone
from pathlib import Path

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import requests
import requests_cache
import networkx as nx
from scipy.stats import hypergeom
from statsmodels.stats.multitest import multipletests
from tenacity import retry, stop_after_attempt, wait_exponential, retry_if_exception_type

# ===== Args =====
ap = argparse.ArgumentParser()
ap.add_argument("--config", required=True, help="JSON config with input_file, a_id, b_id, a_label, b_label, run_name")
args = ap.parse_args()
cfg = json.loads(Path(args.config).read_text())

INPUT = cfg["input_file"]
A_ID, B_ID = cfg["a_id"], cfg["b_id"]
A_LABEL, B_LABEL = cfg["a_label"], cfg["b_label"]
RUN_NAME = cfg["run_name"]
TAXID = cfg.get("taxid", 160488)
ORGANISM = cfg.get("organism", "Pseudomonas putida KT2440")
KEGG_ORG = cfg.get("kegg_org", "ppu")
DE_P = cfg.get("de_p", 0.05)
DE_FC = cfg.get("de_fc", 1.0)
EXPAND_P = cfg.get("expand_p", 0.05)
RUN_PARENT = Path(cfg.get("run_parent", str(Path.home() / "proteomics_runs")))

# ===== Paths =====
ts = datetime.now().strftime("%Y%m%d_%H%M%S")
RUN_DIR = RUN_PARENT / f"proteomics_run_{ts}_{RUN_NAME}"
FIG = RUN_DIR / "figures"; TBL = RUN_DIR / "tables"; ANN = RUN_DIR / "annotations"
for d in (RUN_DIR, FIG, TBL, ANN): d.mkdir(parents=True, exist_ok=True)
print(f"RUN_DIR={RUN_DIR}")

requests_cache.install_cache(str(ANN / "api_cache"), backend="sqlite",
                              expire_after=86400 * 30, allowable_methods=["GET", "POST"])

A_MEAN, B_MEAN = f"log2_mean_{A_ID}", f"log2_mean_{B_ID}"
A_STD, B_STD = f"log2_std_{A_ID}", f"log2_std_{B_ID}"
COLOR_UP, COLOR_DOWN, COLOR_NS = "#D62728", "#1F77B4", "#AAAAAA"

@retry(stop=stop_after_attempt(4), wait=wait_exponential(multiplier=1, min=1, max=15),
       retry=retry_if_exception_type((requests.HTTPError, requests.ConnectionError, requests.Timeout)))
def safe_get(url, **kwargs):
    r = requests.get(url, timeout=60, **kwargs)
    r.raise_for_status()
    return r

# ===== Stage 0: Load + validate =====
print("\n[0] Validation")
h = hashlib.sha256()
with open(INPUT, "rb") as f:
    for chunk in iter(lambda: f.read(2**20), b""): h.update(chunk)
SHA = h.hexdigest()

xl = pd.ExcelFile(INPUT)
df = pd.read_excel(xl, sheet_name="Full t-test output")
required = ["Protein","Protein.Group","Protein.Names","Protein.Description",
            A_MEAN, B_MEAN, A_STD, B_STD,
            "t-test_stat","p-value","p_adjusted(BH)","log2_Fold_change_A/B"]
missing = [c for c in required if c not in df.columns]
assert not missing, f"Missing: {missing}"
diff = (df["log2_Fold_change_A/B"] - (df[A_MEAN] - df[B_MEAN])).abs()
assert diff.max() < 0.01, f"Sign convention mismatch: {diff.max()}"
df["uniprot_acc"] = df["Protein.Group"].apply(lambda s: str(s).split(";")[0].strip())
print(f"  Proteins: {len(df)}, sign convention OK, schema OK")

# Manifest
skill_rev = subprocess.run(["git","-C",str(Path.home()/".claude/skills/proteomics-agent"),"rev-parse","--short","HEAD"],
                            capture_output=True, text=True).stdout.strip()
manifest = {
    "run_id": RUN_DIR.name, "skill_git_rev": skill_rev, "skill": "proteomics-agent",
    "input_files": [{"path": INPUT, "sha256": SHA, "n_proteins": int(len(df))}],
    "conditions": {
        "A_label": A_LABEL, "A_sample_id": A_ID,
        "B_label": B_LABEL, "B_sample_id": B_ID,
        "fold_change_polarity": f"positive log2FC = higher in {A_LABEL}",
        "fold_change_check": f"PASS (max|diff|={float(diff.max()):.6f})",
    },
    "organism": {"name": ORGANISM, "ncbi_taxid": TAXID, "kegg_org_code": KEGG_ORG},
    "random_seed": 42,
    "created_at": datetime.now(timezone.utc).isoformat(),
    "stages_completed": ["00_intake"],
}
(RUN_DIR/"manifest.json").write_text(json.dumps(manifest, indent=2))

# ===== Stage 1: QC =====
print("\n[1] QC")
log2fc = df["log2_Fold_change_A/B"]; padj = df["p_adjusted(BH)"]
n_total, n_named = len(df), df["Protein.Names"].notna().sum()
nan_int = df[[A_MEAN,B_MEAN,A_STD,B_STD]].isna().any(axis=1).sum()
zero_std = ((df[A_STD]==0)|(df[B_STD]==0)).sum()

fig, axes = plt.subplots(1, 2, figsize=(11, 4))
axes[0].hist(df["p-value"].dropna(), bins=50, color="#666", edgecolor="white")
axes[0].axvline(0.05, color="red", ls="--", lw=1)
axes[0].set_xlabel("raw p-value"); axes[0].set_title("p-value distribution")
axes[1].hist(log2fc, bins=60, color="#666", edgecolor="white")
axes[1].axvline(0, color="black", lw=0.8)
axes[1].set_xlabel(f"log2FC ({A_LABEL} / {B_LABEL})"); axes[1].set_title(f"log2FC distribution range=[{log2fc.min():.2f},{log2fc.max():.2f}]")
plt.tight_layout(); plt.savefig(FIG/"01_qc_distributions.png"); plt.close()

print(f"  n_proteins={n_total}, with_names={n_named}({100*n_named/n_total:.0f}%), nan={nan_int}, std0={zero_std}")
print(f"  log2FC range: [{log2fc.min():.2f}, {log2fc.max():.2f}]")

# ===== Stage 4: DE call (strict + expanded) =====
print("\n[4] DE")
strict_mask_up = (padj < DE_P) & (log2fc > DE_FC)
strict_mask_dn = (padj < DE_P) & (log2fc < -DE_FC)
strict_mask = strict_mask_up | strict_mask_dn
expand_mask_up = (padj < EXPAND_P) & (log2fc > 0)
expand_mask_dn = (padj < EXPAND_P) & (log2fc < 0)
expand_mask = expand_mask_up | expand_mask_dn

print(f"  strict (p_adj<{DE_P}, |log2FC|>{DE_FC}): {strict_mask.sum()} (up={strict_mask_up.sum()}, dn={strict_mask_dn.sum()})")
print(f"  expanded (p_adj<{EXPAND_P}, any FC):     {expand_mask.sum()} (up={expand_mask_up.sum()}, dn={expand_mask_dn.sum()})")

def annotated_de(mask, name):
    out = df[mask].copy()
    out["direction"] = ["UP" if x > 0 else "DN" for x in out["log2_Fold_change_A/B"]]
    out["impact_score"] = out["log2_Fold_change_A/B"].abs() * (-np.log10(out["p_adjusted(BH)"].clip(lower=1e-300)))
    out = out.sort_values("impact_score", ascending=False)
    cols = ["Protein","Protein.Group","Protein.Names","Protein.Description",
            A_MEAN,B_MEAN,"log2_Fold_change_A/B","p-value","p_adjusted(BH)","direction","impact_score","uniprot_acc"]
    out[cols].to_csv(TBL/f"04_de_{name}.csv", index=False)
    return out

de_strict = annotated_de(strict_mask, f"strict_p{DE_P}_fc{DE_FC}")
de_expand = annotated_de(expand_mask, f"expanded_p{EXPAND_P}")

# Volcano
fig, ax = plt.subplots(figsize=(8, 7))
mlog10p = -np.log10(padj.clip(lower=1e-300))
colors = np.where(strict_mask_up, COLOR_UP, np.where(strict_mask_dn, COLOR_DOWN, COLOR_NS))
sizes = np.where(strict_mask, 24, 8)
for c, lab in [(COLOR_NS,"ns"),(COLOR_DOWN,f"DN ({strict_mask_dn.sum()})"),(COLOR_UP,f"UP ({strict_mask_up.sum()})")]:
    m = colors == c
    ax.scatter(log2fc[m], mlog10p[m], s=sizes[m], c=c, alpha=0.7, edgecolors="none", label=lab)
ax.axhline(-np.log10(DE_P), color="grey", ls="--", lw=0.8, alpha=0.5)
for x in (-DE_FC, DE_FC):
    ax.axvline(x, color="grey", ls="--", lw=0.5, alpha=0.5)
ax.set_xlabel(f"log2FC ({A_LABEL} / {B_LABEL})"); ax.set_ylabel("-log10 p_adj")
ax.set_title(f"Volcano: {A_LABEL} vs {B_LABEL}\np_adj<{DE_P}, |log2FC|>{DE_FC}: {strict_mask.sum()} sig")
ax.legend(loc="upper left")
top_lab = de_strict.nlargest(min(12, len(de_strict)), "impact_score")
for _, r in top_lab.iterrows():
    nm = str(r["Protein.Names"]) if pd.notna(r["Protein.Names"]) else r["Protein"]
    ax.annotate(nm, (r["log2_Fold_change_A/B"], -np.log10(max(r["p_adjusted(BH)"],1e-300))),
                fontsize=7, xytext=(3,3), textcoords="offset points")
plt.tight_layout(); plt.savefig(FIG/"04_volcano.png"); plt.close()

# ===== Stage 5: UniProt annotation (full background) =====
print("\n[5] UniProt annotation (full 1021 background)")
fields = "accession,id,gene_names,protein_name,ec,go_id,go_p,go_f,go_c,xref_kegg,xref_string,xref_alphafolddb,reviewed"

def fetch_uniprot(accs, batch_size=100):
    rows = []
    for i in range(0, len(accs), batch_size):
        batch = accs[i:i+batch_size]
        q = " OR ".join([f"accession:{a}" for a in batch])
        r = safe_get("https://rest.uniprot.org/uniprotkb/search",
                     params={"query":q,"format":"json","size":500,"fields":fields})
        for entry in r.json().get("results", []):
            acc = entry.get("primaryAccession","")
            gn_obj = entry.get("genes", [])
            gn = ""
            if gn_obj:
                if "geneName" in gn_obj[0]: gn = gn_obj[0]["geneName"].get("value","")
                elif "orderedLocusNames" in gn_obj[0] and gn_obj[0]["orderedLocusNames"]:
                    gn = gn_obj[0]["orderedLocusNames"][0].get("value","")
            pn = entry.get("proteinDescription",{}).get("recommendedName",{}).get("fullName",{}).get("value","")
            go_p, go_f, go_c, kegg, string_ids, af_ids = [], [], [], [], [], []
            for x in entry.get("uniProtKBCrossReferences", []):
                db = x.get("database")
                if db == "GO":
                    gid = x.get("id","")
                    for prop in x.get("properties", []):
                        if prop.get("key") == "GoTerm":
                            t = prop.get("value","")
                            if t.startswith("P:"): go_p.append(f"{gid}|{t[2:]}")
                            elif t.startswith("F:"): go_f.append(f"{gid}|{t[2:]}")
                            elif t.startswith("C:"): go_c.append(f"{gid}|{t[2:]}")
                elif db == "KEGG": kegg.append(x.get("id",""))
                elif db == "STRING": string_ids.append(x.get("id",""))
                elif db == "AlphaFoldDB": af_ids.append(x.get("id",""))
            rows.append({"acc":acc,"gene":gn,"protein_name":pn,
                         "go_p":go_p,"go_f":go_f,"go_c":go_c,
                         "kegg_ids":kegg,"string_id":string_ids,"alphafold_id":af_ids,
                         "reviewed": entry.get("entryType","").startswith("UniProtKB reviewed")})
        time.sleep(0.15)
    return rows

all_accs = df["uniprot_acc"].tolist()
bg_rows = fetch_uniprot(all_accs, batch_size=100)
print(f"  bg annotated: {len(bg_rows)}/{len(all_accs)}")
ann_df = pd.DataFrame([{
    "uniprot_acc": r["acc"], "gene_name": r["gene"], "protein_name": r["protein_name"],
    "kegg_ppu": ";".join([k for k in r["kegg_ids"] if k.startswith(f"{KEGG_ORG}:")]),
    "kegg_all": ";".join(r["kegg_ids"]),
    "string_id": ";".join(r["string_id"]),
    "alphafold_id": ";".join(r["alphafold_id"]),
    "go_bp": ";".join(r["go_p"]), "go_mf": ";".join(r["go_f"]), "go_cc": ";".join(r["go_c"]),
    "reviewed": r["reviewed"],
} for r in bg_rows])
ann_df.to_csv(TBL/"05_annotations_all.csv", index=False)
de_with_ann = de_expand.merge(ann_df, on="uniprot_acc", how="left")
de_with_ann.to_csv(TBL/"05_de_with_annotation.csv", index=False)
n_sp = ann_df["reviewed"].sum()
print(f"  SwissProt-reviewed in background: {n_sp}/{len(ann_df)}")

# ===== Stage 6: GO over-representation =====
print("\n[6] GO over-representation (custom background)")
acc2go = {"BP":defaultdict(set),"MF":defaultdict(set),"CC":defaultdict(set)}
term2genes = {"BP":defaultdict(set),"MF":defaultdict(set),"CC":defaultdict(set)}
for r in bg_rows:
    a = r["acc"]
    for cat,key in [("BP","go_p"),("MF","go_f"),("CC","go_c")]:
        for entry in r[key]:
            if "|" not in entry: continue
            tid,tname = entry.split("|",1)
            acc2go[cat][a].add((tid,tname))
            term2genes[cat][(tid,tname)].add(a)
bg_set = {a for a in {r["acc"] for r in bg_rows} if any(acc2go[c][a] for c in ["BP","MF","CC"])}

up_set = set(de_expand[de_expand["log2_Fold_change_A/B"]>0]["uniprot_acc"]) & bg_set
dn_set = set(de_expand[de_expand["log2_Fold_change_A/B"]<0]["uniprot_acc"]) & bg_set

def enrich(de_set, t2g, bg_size, min_t=3, max_t=400, min_o=2):
    rows = []
    for (tid,tname), tg in t2g.items():
        tg_bg = tg & bg_set
        K = len(tg_bg)
        if K < min_t or K > max_t: continue
        n = len(de_set); k = len(de_set & tg_bg)
        if k < min_o: continue
        p = hypergeom.sf(k-1, bg_size, K, n)
        rows.append({"term_id":tid,"term_name":tname,"n_overlap":k,"n_term":K,"n_de":n,
                     "p_value":p,"fold":(k/n)/(K/bg_size) if (n>0 and K>0) else 0,
                     "genes":";".join(sorted(de_set & tg_bg))})
    if not rows: return pd.DataFrame()
    out = pd.DataFrame(rows)
    out["p_adj"] = multipletests(out["p_value"], method="fdr_bh")[1]
    return out.sort_values("p_adj")

bg_size = len(bg_set)
all_e = []
for cat in ["BP","MF","CC"]:
    for label, dset in [("up",up_set),("dn",dn_set)]:
        e = enrich(dset, term2genes[cat], bg_size)
        if len(e): e["category"] = cat; e["direction"] = label; all_e.append(e)
if all_e:
    enr = pd.concat(all_e, ignore_index=True)
    enr.to_csv(TBL/"06_enrichment_all.csv", index=False)
    sig = enr[enr["p_adj"]<0.05]
    sig.to_csv(TBL/"06_enrichment_sig.csv", index=False)
    print(f"  GO sig (p_adj<0.05): UP={(sig['direction']=='up').sum()} DN={(sig['direction']=='dn').sum()}")
    for d in ["up","dn"]:
        sub = sig[(sig.direction==d)&(sig.category=="BP")].head(8)
        if len(sub):
            print(f"  {d.upper()} BP top:")
            for _,r in sub.iterrows():
                print(f"    p_adj={r['p_adj']:.2e} fold={r['fold']:.1f}x {r['n_overlap']}/{r['n_term']} {r['term_name'][:55]}")
else:
    print("  No GO enrichments")

# ===== Stage 7: KEGG pathway enrichment =====
print("\n[7] KEGG pathway enrichment")
links_r = safe_get(f"https://rest.kegg.jp/link/pathway/{KEGG_ORG}")
pw2g = defaultdict(set)
for line in links_r.text.strip().split("\n"):
    p = line.split("\t")
    if len(p)==2: pw2g[p[1]].add(p[0])
pn_r = safe_get(f"https://rest.kegg.jp/list/pathway/{KEGG_ORG}")
pname = {}
for line in pn_r.text.strip().split("\n"):
    p = line.split("\t")
    if len(p)==2:
        pid = "path:"+p[0] if not p[0].startswith("path:") else p[0]
        pname[pid] = p[1].split(" - ")[0]

acc2kegg = {}
for r in bg_rows:
    kegg_org = [k for k in r["kegg_ids"] if k.startswith(f"{KEGG_ORG}:")]
    if kegg_org: acc2kegg[r["acc"]] = kegg_org[0]
bg_kegg = set(acc2kegg.values())
up_kegg = {acc2kegg[a] for a in up_set if a in acc2kegg}
dn_kegg = {acc2kegg[a] for a in dn_set if a in acc2kegg}
print(f"  bg KEGG: {len(bg_kegg)}, UP: {len(up_kegg)}, DN: {len(dn_kegg)}")

def enrich_kegg(de_kegg, label):
    rows = []
    for path, pg in pw2g.items():
        pgb = pg & bg_kegg
        K = len(pgb)
        if K < 3 or K > 400: continue
        n = len(de_kegg); k = len(de_kegg & pgb)
        if k < 2: continue
        p = hypergeom.sf(k-1, len(bg_kegg), K, n)
        rows.append({"pathway_id":path.replace("path:",""),"pathway_name":pname.get(path,path),
                     "n_overlap":k,"n_pathway":K,"n_de":n,
                     "p_value":p,"fold":(k/n)/(K/len(bg_kegg)) if (n>0 and K>0) else 0,
                     "direction":label,"genes":";".join(sorted(de_kegg & pgb))})
    if not rows: return pd.DataFrame()
    out = pd.DataFrame(rows)
    out["p_adj"] = multipletests(out["p_value"], method="fdr_bh")[1]
    return out.sort_values("p_adj")

ku = enrich_kegg(up_kegg, "up"); kd = enrich_kegg(dn_kegg, "dn")
all_kegg = pd.concat([ku,kd], ignore_index=True) if (len(ku)+len(kd))>0 else pd.DataFrame()
if len(all_kegg):
    sig_k = all_kegg[all_kegg["p_adj"]<0.05]
    sig_k.to_csv(TBL/"07_kegg_pathway_enrichment.csv", index=False)
    print(f"  KEGG sig: {len(sig_k)}")
    for _, r in sig_k.iterrows():
        print(f"    {r['direction'].upper()} {r['pathway_id']} p_adj={r['p_adj']:.2e} fold={r['fold']:.1f}x {r['n_overlap']}/{r['n_pathway']} {r['pathway_name']}")

# Pathway overlay
overlay = []
for path, pg in pw2g.items():
    pgb = pg & bg_kegg
    nu = len(pgb & up_kegg); nd = len(pgb & dn_kegg)
    if nu+nd >= 2:
        overlay.append({"pathway_id":path.replace("path:",""),"pathway_name":pname.get(path,path),
                        "n_pathway":len(pgb),"n_de_up":nu,"n_de_down":nd,"n_de_total":nu+nd,
                        "up_genes":";".join(sorted(pgb&up_kegg)),"down_genes":";".join(sorted(pgb&dn_kegg))})
overlay_df = pd.DataFrame(overlay).sort_values("n_de_total", ascending=False) if overlay else pd.DataFrame()
if len(overlay_df):
    overlay_df.to_csv(TBL/"07_pathway_overlay.csv", index=False)
    print(f"  Top 6 pathways with ≥2 DE hits:")
    for _, r in overlay_df.head(6).iterrows():
        print(f"    {r['pathway_id']:<10} UP={r['n_de_up']:>2} DN={r['n_de_down']:>2} {r['pathway_name'][:50]}")

# ===== Stage 8: STRING network =====
print("\n[8] STRING network")
de_for_string = de_expand
ids = "%0d".join(de_for_string["uniprot_acc"].tolist())
r = safe_get(f"https://string-db.org/api/json/get_string_ids?identifiers={ids}&species={TAXID}")
mapping = r.json()
acc2sid = {}; sid_info = {}
for m in mapping:
    q, sid = m.get("queryItem",""), m.get("stringId","")
    if q and sid:
        acc2sid[q] = sid
        sid_info[sid] = {"name": m.get("preferredName","")}
ver = safe_get("https://string-db.org/api/json/version").json()[0].get("string_version","?")

string_ids = list(acc2sid.values())
ne = safe_get(f"https://string-db.org/api/json/network?identifiers={'%0d'.join(string_ids)}&species={TAXID}&required_score=700&network_type=physical")
edges = ne.json()
G = nx.Graph()
for _, row in de_for_string.iterrows():
    a = row["uniprot_acc"]; sid = acc2sid.get(a)
    if sid:
        gene = sid_info[sid]["name"] or row.get("Protein.Names") or row["Protein"]
        G.add_node(sid, gene=str(gene), acc=a, log2FC=float(row["log2_Fold_change_A/B"]),
                   p_adj=float(row["p_adjusted(BH)"]),
                   direction=("up" if row["log2_Fold_change_A/B"]>0 else "dn"))
for e in edges:
    a,b = e["stringId_A"], e["stringId_B"]
    if a in G.nodes and b in G.nodes:
        G.add_edge(a, b, score=e.get("score",0))
print(f"  STRING v{ver}: {G.number_of_nodes()} nodes, {G.number_of_edges()} edges")
deg = dict(G.degree())
hubs = sorted([n for n in G.nodes if deg[n]>0], key=lambda n: deg[n], reverse=True)[:15]
hub_rows = []
for n in hubs:
    d = G.nodes[n]
    hub_rows.append({"string_id":n,"gene":d["gene"],"acc":d["acc"],"degree":deg[n],
                     "log2FC":d["log2FC"],"p_adj":d["p_adj"],"direction":d["direction"]})
pd.DataFrame(hub_rows).to_csv(TBL/"08_hubs.csv", index=False)
print(f"  Top 5 hubs: {[(h['gene'],h['degree'],h['direction']) for h in hub_rows[:5]]}")

# Communities
H = G.subgraph(max(nx.connected_components(G), key=len)).copy() if G.number_of_edges()>0 else G
comm_rows = []
if len(H) >= 4:
    try:
        comms = sorted(nx.community.louvain_communities(H, seed=42), key=len, reverse=True)
        for i, c in enumerate(comms):
            if len(c) < 3: continue
            members = sorted([G.nodes[n]["gene"] for n in c])
            nu = sum(1 for n in c if G.nodes[n]["direction"]=="up")
            nd_ = sum(1 for n in c if G.nodes[n]["direction"]=="dn")
            comm_rows.append({"community":i+1,"size":len(c),"n_up":nu,"n_down":nd_,"members":";".join(members)})
        pd.DataFrame(comm_rows).to_csv(TBL/"08_communities.csv", index=False)
        print(f"  Communities ≥3: {len(comm_rows)}")
        for c in comm_rows[:3]:
            print(f"    Comm{c['community']} ({c['size']}n, UP={c['n_up']}, DN={c['n_down']}): {';'.join(c['members'].split(';')[:6])}")
    except Exception as e: print(f"  Louvain failed: {e}")

# Static network viz
H2 = G.subgraph([n for n in G.nodes if deg[n]>0]).copy() if G.number_of_edges()>0 else G
if H2.number_of_nodes() > 0:
    fig, ax = plt.subplots(figsize=(11, 9))
    pos = nx.spring_layout(H2, seed=42, k=0.7, iterations=80)
    nc = ["#D62728" if H2.nodes[n]["direction"]=="up" else "#1F77B4" for n in H2.nodes]
    ns = [80 + 250*abs(H2.nodes[n]["log2FC"]) for n in H2.nodes]
    ew = [(H2.edges[e]["score"]/1000)**2 * 2.5 for e in H2.edges]
    nx.draw_networkx_edges(H2, pos, width=ew, alpha=0.3, edge_color="grey")
    nx.draw_networkx_nodes(H2, pos, node_color=nc, node_size=ns, alpha=0.85, edgecolors="black", linewidths=0.5)
    nx.draw_networkx_labels(H2, pos, labels={n:H2.nodes[n]["gene"] for n in H2.nodes}, font_size=7)
    ax.set_title(f"STRING physical (score≥0.7): {H2.number_of_nodes()} nodes, {H2.number_of_edges()} edges\nRed=UP in {A_LABEL}, Blue=DOWN")
    ax.axis("off"); plt.tight_layout(); plt.savefig(FIG/"08_string_network.png"); plt.close()

# ===== Update manifest =====
mp = RUN_DIR / "manifest.json"
m = json.loads(mp.read_text())
m["stages_completed"] = ["00_intake","01_qc","04_de","05_annotation","06_go_enrichment","07_kegg_enrichment","08_string"]
m["skipped_stages"] = ["09_alphafold (per user)","10_ptm_annotation (per user)"]
m["qc"] = {"n_proteins":int(n_total),"n_named":int(n_named),"nan_rows":int(nan_int),"std_zero":int(zero_std),
           "log2fc_range":[float(log2fc.min()),float(log2fc.max())]}
m["de"] = {"strict":{"p":DE_P,"fc":DE_FC,"n":int(strict_mask.sum()),
                     "up":int(strict_mask_up.sum()),"dn":int(strict_mask_dn.sum())},
           "expanded":{"p":EXPAND_P,"n":int(expand_mask.sum()),
                       "up":int(expand_mask_up.sum()),"dn":int(expand_mask_dn.sum())}}
m["annotation"] = {"bg_n":int(len(ann_df)), "swissprot":int(n_sp)}
m["enrichment"] = {"go_sig_up": int((sig['direction']=='up').sum()) if all_e else 0,
                   "go_sig_dn": int((sig['direction']=='dn').sum()) if all_e else 0,
                   "kegg_sig": int(len(sig_k)) if len(all_kegg) else 0}
m["string"] = {"version":ver,"n_nodes":int(G.number_of_nodes()),"n_edges":int(G.number_of_edges()),
               "n_communities_ge3": len(comm_rows)}
mp.write_text(json.dumps(m, indent=2, default=str))

print(f"\n✅ Pipeline complete: {RUN_DIR}")
print(f"   tables: {len(list(TBL.glob('*')))} files; figures: {len(list(FIG.glob('*')))} files")
