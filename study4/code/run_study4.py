"""Study 4 runner (SAP v4; evaluation exceptions handled as in deviation D-34). Steps: (1) exact references and Shapely cells for all candidate inputs;
(2) the unmutated subject on all inputs, which selects the valid base inputs (rule in SAP v4);
(3) every first-order mutant on all runs of the valid base inputs.
Resumable (JSONL appends, done keys skipped), parallel (2 workers), progress printed."""
import json, os, sys, pickle, time, types
from multiprocessing import Pool
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from s4_common import *
import s4_common as C
import hp_voronoi as H
import mutate as MU
REFD = os.path.join(RES, "ref"); os.makedirs(REFD, exist_ok=True)


def key(bid, mr): return f"{bid}|{mr}"


def ref_task(k):
    bid, mr = k.rsplit("|", 1)
    P, back, perm = build(bid, mr); bx = vc.clip_box(P)
    fn = os.path.join(REFD, k.replace("|", "__") + ".pkl")
    if not os.path.exists(fn):
        cells, _, _ = vc.reference_cells(P, bx)
        sh = vc.sut_shapely(P, bx)
        pickle.dump({"ref": cells, "shapely": [None if g is None else g.wkb for g in sh]}, open(fn + ".part", "wb"))
        os.replace(fn + ".part", fn)
    return k


def load_ref(k):
    import shapely
    d = pickle.load(open(os.path.join(REFD, k.replace("|", "__") + ".pkl"), "rb"))
    return d["ref"], [None if b is None else shapely.from_wkb(b) for b in d["shapely"]]


def eval_family(fn, bid):
    """Run subject fn on the base input and its five follow-ups; return per-run and per-MR results."""
    out = {"runs": {}, "mr": {}}
    cells_by = {}
    for mr in ["base"] + MRS:
        P, back, perm = build(bid, mr); bx = vc.clip_box(P)
        st, cells = run_subject(fn, P, bx)
        ref, sh = load_ref(key(bid, mr))
        rec = {"status": st}
        if st == "OK":
            # D-34: an exception while EVALUATING a completed run marks that step "not evaluable"
            # and counts as a failure/detection (favours the oracle-free tests).
            rec["evalerr"] = []
            try:
                polys = to_polys(cells)
            except Exception:
                polys = [None] * len(P); rec["evalerr"].append("to_polys")
            try:
                r = vc.compare_to_reference(P, polys, ref, bx)
                rec.update({"missing": r["missing"], "misassign": r["misassign"], "worst": r["worst"]})
                rec["ref_fail"] = bool(r["missing"] or r["misassign"] or r["worst"] > TAU)
            except Exception:
                rec.update({"missing": 0, "misassign": 0, "worst": float("inf"), "ref_fail": True}); rec["evalerr"].append("ref")
            try:
                rec["def"] = def_checks(P, polys, bx)
            except Exception:
                rec["def"] = {"P1": None, "P2": None, "P3": None, "flag": True}; rec["evalerr"].append("def")
            pairs = [(i, i) for i in range(len(P))]
            try:
                rec["diff"] = vc.cells_distance(polys, sh, bx, pairs)
            except Exception:
                rec["diff"] = float("inf"); rec["evalerr"].append("diff")
            cells_by[mr] = (polys, back, perm, P)
        out["runs"][mr] = rec
    if "base" in cells_by:
        A, _, _, P0 = cells_by["base"]; bx0 = vc.clip_box(P0)
        for mr in MRS:
            if mr not in cells_by:
                out["mr"][mr] = None; continue
            B, back, perm, Q = cells_by[mr]
            try:
                Bb = [None if g is None else back(g) for g in B]
            except Exception:
                out["mr"][mr] = float("inf"); out.setdefault("mr_evalerr", []).append(mr); continue
            pairs = [(i, i) for i in range(len(P0))] if perm is None else [(int(perm[k]), k) for k in range(len(Q))]
            try:
                out["mr"][mr] = vc.cells_distance(A, Bb, bx0, pairs)
            except Exception:
                out["mr"][mr] = float("inf"); out.setdefault("mr_evalerr", []).append(mr)
    return out


def orig_task(bid):
    return {"bid": bid, **eval_family(H.voronoi_cells, bid)}


def mutant_task(args):
    mid, src, bids = args
    mod = types.ModuleType("hp_mut")
    try:
        exec(compile(src, "hp_mut", "exec"), mod.__dict__)
    except Exception as e:
        return {"mid": mid, "compile_error": type(e).__name__}
    return {"mid": mid, "fam": {b: eval_family(mod.voronoi_cells, b) for b in bids}}


def valid(fam):
    for mr, r in fam["runs"].items():
        if r["status"] != "OK" or r["ref_fail"] or r["def"]["flag"] or r["diff"] > TAU:
            return False
    return all(d is not None and d <= TAU for d in fam["mr"].values())


def stream(tasks, func, out, label, keyf):
    done = set()
    if os.path.exists(out):
        done = {keyf(json.loads(l)) for l in open(out)}
    todo = [t for t in tasks if (t if isinstance(t, str) else t[0]) not in done]
    print(f"[{label}] {len(tasks)} tasks, {len(todo)} to do", flush=True)
    t0 = time.time()
    with Pool(2) as pool, open(out, "a") as f:
        for i, r in enumerate(pool.imap_unordered(func, todo), 1):
            f.write(json.dumps(r) + "\n"); f.flush()
            if i % max(1, len(todo) // 20) == 0 or i == len(todo):
                el = time.time() - t0
                print(f"[{label}] {i}/{len(todo)} ({100*i/len(todo):.0f}%) {el:.0f}s ETA {el/i*(len(todo)-i):.0f}s", flush=True)


def main():
    keys = [key(b, m) for b in bases() for m in ["base"] + MRS]
    stream(keys, ref_task, os.path.join(RES, "ref_done.jsonl"), "ref", lambda r: r) if False else None
    done = set()
    fdone = os.path.join(RES, "ref_done.txt")
    if os.path.exists(fdone): done = set(open(fdone).read().split())
    todo = [k for k in keys if k not in done]
    print(f"[ref] {len(keys)} inputs, {len(todo)} to do", flush=True)
    with Pool(2) as pool, open(fdone, "a") as f:
        for i, k in enumerate(pool.imap_unordered(ref_task, todo), 1):
            f.write(k + "\n"); f.flush()
            if i % 24 == 0 or i == len(todo): print(f"[ref] {i}/{len(todo)}", flush=True)
    stream(bases(), orig_task, os.path.join(RES, "orig.jsonl"), "orig", lambda r: r["bid"])
    orig = [json.loads(l) for l in open(os.path.join(RES, "orig.jsonl"))]
    vb = sorted(r["bid"] for r in orig if valid(r))
    json.dump({"valid_bases": vb, "n_candidates": len(orig)}, open(os.path.join(RES, "valid_bases.json"), "w"), indent=1)
    print(f"valid base inputs {len(vb)} of {len(orig)}", flush=True)
    M = MU.mutants()
    json.dump([m for m, _ in M], open(os.path.join(RES, "mutant_ids.json"), "w"))
    stream([(m, s, vb) for m, s in M], mutant_task, os.path.join(RES, "mutants.jsonl"), "mutants", lambda r: r["mid"])


if __name__ == "__main__":
    main()
