"""Study 4 analysis exactly as specified in SAP v4 (decision rules, alpha 0.01 each)."""
import json, os, sys
import numpy as np
from scipy import stats
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from s4_common import RES, MRS, EXACT_MRS, TAU
ALPHA = 0.01


def cp_lower(k, n, a=ALPHA): return 0.0 if k == 0 else float(stats.beta.ppf(a, k, n - k + 1))
def cp_upper(k, n, a=ALPHA): return 1.0 if k == n else float(stats.beta.ppf(1 - a, k + 1, n - k))


def main():
    M = [json.loads(l) for l in open(os.path.join(RES, "mutants.jsonl"))]
    rows = []
    for m in M:
        if "compile_error" in m:
            rows.append({"mid": m["mid"], "class": "compile_error"}); continue
        runs = [r for f in m["fam"].values() for r in f["runs"].values()]
        nonsilent = any(r["status"] != "OK" for r in runs)
        ref_kill = any(r["status"] == "OK" and r["ref_fail"] for r in runs)
        worst = max([r["worst"] for r in runs if r["status"] == "OK"] + [0.0])
        anyerr = any(r["status"] == "OK" and (r["missing"] or r["misassign"]) for r in runs)
        d = any(r["status"] == "OK" and r["def"]["flag"] for r in runs)
        diff = any(r["status"] == "OK" and r["diff"] > TAU for r in runs)
        mr = {x: any((f["mr"].get(x) is not None and f["mr"][x] > TAU) for f in m["fam"].values()) for x in MRS}
        cls = "nonsilent" if nonsilent else ("silent_killed" if ref_kill else "not_killed")
        rows.append({"mid": m["mid"], "class": cls, "def": d, "diff": diff, "mr": mr,
                     "exact_mr": any(mr[x] for x in EXACT_MRS), "all_mr": any(mr.values()),
                     "worst": worst, "cell_error": anyerr})
    S = [r for r in rows if r["class"] == "silent_killed"]; n = len(S)
    kd = sum(r["def"] for r in S); kdiff = sum(r["diff"] for r in S); kem = sum(r["exact_mr"] for r in S)
    union = [r["def"] or r["all_mr"] for r in S]; kmiss = n - sum(union)
    res = {"n_mutants": len(rows), "classes": {c: sum(r["class"] == c for r in rows) for c in ("compile_error", "nonsilent", "silent_killed", "not_killed")},
           "n_silent_killed": n}
    res["P1"] = {"k": kd, "n": n, "recall": kd / n if n else None, "upper": cp_upper(kd, n) if n else None,
                 "supported": bool(n and cp_upper(kd, n) < 0.95)}
    b = sum(1 for r in S if r["def"] and not r["exact_mr"]); c = sum(1 for r in S if r["exact_mr"] and not r["def"])
    p2 = float(stats.binomtest(b, b + c, 0.5, alternative="greater").pvalue) if b + c else 1.0
    res["P2"] = {"def_only": b, "mr_only": c, "exact_mr_k": kem, "p": p2, "supported": p2 < ALPHA}
    res["P3"] = {"missed": kmiss, "n": n, "rate": kmiss / n if n else None, "lower": cp_lower(kmiss, n) if n else None,
                 "supported": bool(n and cp_lower(kmiss, n) >= 0.05)}
    res["P4"] = {"k": kdiff, "n": n, "recall": kdiff / n if n else None, "lower": cp_lower(kdiff, n) if n else None,
                 "supported": bool(n and cp_lower(kdiff, n) >= 0.90)}
    miss = [r["worst"] for r, u in zip(S, union) if not u]; det = [r["worst"] for r, u in zip(S, union) if u]
    if len(miss) >= 3 and det:
        p5 = float(stats.mannwhitneyu(miss, det, alternative="less").pvalue)
    else:
        p5 = None
    res["P5"] = {"n_missed": len(miss), "n_detected": len(det), "median_missed": float(np.median(miss)) if miss else None,
                 "median_detected": float(np.median(det)) if det else None, "p": p5, "supported": p5 is not None and p5 < ALPHA}
    res["per_mr_recall"] = {x: sum(r["mr"][x] for r in S) for x in MRS}
    json.dump({"summary": res, "rows": rows}, open(os.path.join(RES, "study4_verdicts.json"), "w"), indent=1)
    print(json.dumps(res, indent=1))


if __name__ == "__main__":
    main()
