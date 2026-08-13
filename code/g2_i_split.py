"""Split the single-expression g2_i_f.inc (850 KB, written by g2_i.py)
into per-term Fortran statements accumulating into fIv.

flang cannot compile the original single expression in reasonable time
(>20 min, >30 GB RAM at any -O level); the split form compiles in
seconds.  Large parenthesized subexpressions are lifted into tv<N>
temporaries (declared in g2_i_f_split_decl.inc) and depth-0 sums become
`acc = acc + (term)` chains.  The output is verified against the SymPy
integrand at rational test points to double-precision roundoff.
"""
import re

LIMIT = 4000
stmts = []
cnt = [0]


def split_sum(e):
    """Split e at depth-0 binary +/- boundaries."""
    terms, depth, start = [], 0, 0
    for i, ch in enumerate(e):
        if ch == "(":
            depth += 1
        elif ch == ")":
            depth -= 1
        elif ch in "+-" and depth == 0 and i > start:
            prev = e[i - 1]
            if prev in "edED" and i >= 2 and e[i - 2].isdigit():
                continue                      # exponent sign
            j = i - 1
            while j >= 0 and e[j] == " ":
                j -= 1
            if j >= 0 and e[j] in "*/(+-":
                continue                      # unary/operator context
            terms.append(e[start:i].strip())
            start = i
    terms.append(e[start:].strip())
    return terms


def new_temp():
    cnt[0] += 1
    return "tv%d" % cnt[0]


def reduce_expr(e):
    """Return an expression string of len < LIMIT computing e, emitting
    temporary-variable statements as needed."""
    e = e.strip()
    if len(e) < LIMIT:
        return e
    terms = split_sum(e)
    if len(terms) > 1:
        name = new_temp()
        stmts.append("%s = 0" % name)
        for t in terms:
            sign = "+"
            if t[0] in "+-":
                sign, t = t[0], t[1:].strip()
            stmts.append("%s = %s %s (%s)"
                         % (name, name, sign, reduce_expr(t)))
        return name
    # a single product: lift its big parenthesized factors
    out, depth, i, n, last = [], 0, 0, len(e), 0
    gstart = 0
    while i < n:
        if e[i] == "(":
            if depth == 0:
                gstart = i
            depth += 1
        elif e[i] == ")":
            depth -= 1
            if depth == 0:
                inner = e[gstart + 1:i]
                if len(inner) >= LIMIT:
                    out.append(e[last:gstart])
                    rep = reduce_expr(inner)
                    if len(rep) >= LIMIT:
                        nm = new_temp()
                        stmts.append("%s = %s" % (nm, rep))
                        rep = nm
                    out.append("(" + rep + ")")
                    last = i + 1
        i += 1
    out.append(e[last:])
    res = "".join(out)
    if len(res) >= LIMIT:
        nm = new_temp()
        stmts.append("%s = %s" % (nm, res))
        res = nm
    return res


def wrap(stmt):
    """Wrap to free-form Fortran with leading-& continuations (safe to
    break anywhere, including inside tokens)."""
    out, line = [], ""
    for tok in re.split(r"(\s+|\*|\+|/|-|,)", stmt):
        if not tok:
            continue
        if len(line) + len(tok) > 66:
            out.append(line + "&")
            line = "      &" + tok
        else:
            line += tok
    out.append(line)
    return "\n".join(out)


def main():
    src = open("g2_i_f.inc").read()
    body = "\n".join(src.splitlines()[1:])
    body = re.sub(r"&\n\s*", "", body)
    assert body.startswith("fIv = ")
    expr = body[len("fIv = "):].replace("\n", " ")
    stmts.append("fIv = " + reduce_expr(expr))
    print("split into %d statements, %d temps" % (len(stmts), cnt[0]))
    text = "! generated from g2_i_f.inc -- do not edit\n"
    text += "".join(wrap(s) + "\n" for s in stmts)
    open("g2_i_f_split.inc", "w").write(text)
    decl = ("complex(dp) :: "
            + ", ".join("tv%d" % k for k in range(1, cnt[0] + 1)) + "\n")
    open("g2_i_f_split_decl.inc", "w").write(decl)
    # quad-precision twin: all literals are exact small integers written
    # as N.0d0, so promoting them to _qp is a plain rewrite.  Needed
    # because double precision has no significance left where the outer
    # b_hat = (r+v)(1-r-v) vanishes (see g2_i_qmc.f90).
    qtext = re.sub(r"(\d+\.\d+)d0", r"\1_qp", text)
    assert not re.search(r"\dd0", qtext), "unconverted double literal"
    open("g2_i_f_split_qp.inc", "w").write(qtext)
    open("g2_i_f_split_qp_decl.inc", "w").write(
        decl.replace("complex(dp)", "complex(qp)"))


if __name__ == "__main__":
    main()
