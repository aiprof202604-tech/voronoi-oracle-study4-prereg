"""Study 4 descriptive counts promised by deviations D-34, D-36 and D-37 (not new analyses):
runs and mutants that were not evaluable or outside the box, mutants stopped early, and how many
silent reference-killed (SRK) mutants the definition checks detect without the D-34/D-36 rule."""
import json, os, sys
from scipy import stats
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from s4_common import RES
V = json.load(open(os.path.join(RES, "study4_verdicts.json")))
M = {json.loads(l)["mid"]: json.loads(l) for l in open(os.path.join(RES, "mutants.jsonl"))}
runs = lambda m: [r for f in M[m].get("fam", {}).values() for r in f["runs"].values()]
S = [r["mid"] for r in V["rows"] if r["class"] == "silent_killed"]
gen = sum(any(r["status"] == "OK" and not r.get("outbox") and not r.get("evalerr") and r["def"]["flag"] for r in runs(m)) for m in S)
w = [r["worst"] for r in V["rows"] if r["class"] == "silent_killed"]
out = {"mutants": len(M), "valid_bases": len(json.load(open(os.path.join(RES, "valid_bases.json")))["valid_bases"]),
       "stopped_early": sum(m.get("stopped_early", False) for m in M.values()),
       "outbox_runs": sum(r.get("outbox", False) for m in M for r in runs(m)),
       "outbox_mutants": sum(any(r.get("outbox", False) for r in runs(m)) for m in M),
       "outbox_srk": sum(any(r.get("outbox", False) for r in runs(m)) for m in S),
       "evalerr_runs": sum(bool(r.get("evalerr")) for m in M for r in runs(m)),
       "srk": len(S), "def_genuine": gen,
       "def_genuine_upper99": float(stats.beta.ppf(0.99, gen + 1, len(S) - gen)) if gen < len(S) else 1.0,
       "srk_min_worst": min(w), "srk_median_worst": sorted(w)[len(w) // 2]}
json.dump(out, open(os.path.join(RES, "s4_report.json"), "w"), indent=1)
print(json.dumps(out, indent=1))
