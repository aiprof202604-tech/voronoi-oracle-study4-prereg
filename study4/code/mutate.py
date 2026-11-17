"""Study 4: deterministic generation of all first-order mutants of hp_voronoi.py.
Operators (applied to every eligible node, in source order):
  AOR  binary arithmetic operator replacement among + - * /           (each other operator)
  ROR  comparison operator replacement among < <= > >= == !=          (each other operator)
  CRP  numeric constant replacement c -> c+1, c-1, 0 (0 skipped if c == 0; duplicates removed)
  LCR  boolean operator replacement and <-> or
  NEG  negation of a whole if/condition test
  SDL  deletion of one statement (replaced by `pass`) except def, return and the docstring
Mutant id: f"{op}:{lineno}:{col}:{end_col}:{k}". Output: list of (id, source) in a fixed order."""
import ast, copy, json, os, sys

SRC = os.path.join(os.path.dirname(os.path.abspath(__file__)), "hp_voronoi.py")
ARITH = [ast.Add, ast.Sub, ast.Mult, ast.Div]
CMP = [ast.Lt, ast.LtE, ast.Gt, ast.GtE, ast.Eq, ast.NotEq]


def _sites(tree):
    sites = []
    for idx, node in enumerate(ast.walk(tree)):
        if isinstance(node, ast.BinOp) and type(node.op) in ARITH:
            for k, op in enumerate(o for o in ARITH if o is not type(node.op)):
                sites.append(("AOR", idx, node, k, ("op", op)))
        elif isinstance(node, ast.Compare) and len(node.ops) == 1 and type(node.ops[0]) in CMP:
            for k, op in enumerate(o for o in CMP if o is not type(node.ops[0])):
                sites.append(("ROR", idx, node, k, ("cmp", op)))
        elif isinstance(node, ast.Constant) and isinstance(node.value, (int, float)) and not isinstance(node.value, bool):
            v = node.value
            new = []
            for w in (v + 1, v - 1, 0 * v):
                if w != v and w not in new:
                    new.append(w)
            for k, w in enumerate(new):
                sites.append(("CRP", idx, node, k, ("const", w)))
        elif isinstance(node, ast.BoolOp):
            sites.append(("LCR", idx, node, 0, ("bool", ast.Or if isinstance(node.op, ast.And) else ast.And)))
        elif isinstance(node, (ast.If, ast.While)):
            sites.append(("NEG", idx, node, 0, ("neg", None)))
    return sites


def _stmts(tree):
    out = []
    for pidx, node in enumerate(ast.walk(tree)):
        for field in ("body", "orelse"):
            body = getattr(node, field, None)
            if not isinstance(body, list):
                continue
            for idx, st in enumerate(body):
                if isinstance(st, (ast.FunctionDef, ast.Return)):
                    continue
                if isinstance(st, ast.Expr) and isinstance(st.value, ast.Constant) and isinstance(st.value.value, str):
                    continue
                out.append((pidx, field, idx, st))
    return out


def mutants():
    """Nodes are located in the deep copy by their index in ast.walk order, which is identical
    for the original tree and its copy (positions alone are ambiguous for nested operators)."""
    src = open(SRC).read()
    base = ast.parse(src)
    res = []
    for op, idx, node, k, (kind, val) in _sites(base):
        t = copy.deepcopy(base)
        m = list(ast.walk(t))[idx]
        if kind == "op": m.op = val()
        elif kind == "cmp": m.ops = [val()]
        elif kind == "const": m.value = val
        elif kind == "bool": m.op = val()
        elif kind == "neg": m.test = ast.UnaryOp(op=ast.Not(), operand=m.test)
        res.append((f"{op}:{node.lineno}:{node.col_offset}:{node.end_col_offset}:{k}", ast.unparse(ast.fix_missing_locations(t))))
    for pidx, field, idx, st in _stmts(base):
        t = copy.deepcopy(base)
        m = list(ast.walk(t))[pidx]
        getattr(m, field)[idx] = ast.Pass()
        res.append((f"SDL:{st.lineno}:{st.col_offset}:{st.end_col_offset}:0", ast.unparse(ast.fix_missing_locations(t))))
    seen, out = set(), []
    orig = ast.unparse(base)
    for mid, s in res:
        if s == orig or s in seen:
            continue
        seen.add(s); out.append((mid, s))
    return out


if __name__ == "__main__":
    print(json.dumps({"n_mutants": len(mutants())}))
