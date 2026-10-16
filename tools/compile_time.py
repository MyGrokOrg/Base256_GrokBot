#!/usr/bin/env python3
"""Compile-time and header-size measurement for the base256 headers.

usage: compile_time.py [--cxx g++-14,clang++-19] [-n 8] [--ref REV ...]
                       [--json out.json] [--compare BASELINE.json] [--threshold 0.10]

Compiles each representative TU (below) with -std=c++23 -O2 -c, N times per
compiler, and prints the mean and standard deviation of the compiler's CPU time
(user + sys, robust against other load) and the wall-clock mean, plus the byte
size of each header group. --ref REV measures the headers of a git revision
(extracted with `git archive`) in the same run, interleaved with the working
tree, so revisions can be compared on one machine.

Header groups: "core" is everything <base256/base256.hpp> includes (transitively),
"complex" is what <base256/complex.hpp> adds to it (complex.hpp and
detail/cx_engine.hpp). Each TU belongs to one group. --compare BASELINE.json
measures each group's baseline revision ("groups" in the JSON) next to the working
tree, compares that group's TUs and header bytes with it, and prints a GitHub
`::warning::` for every TU whose mean CPU time grew by more than --threshold
(default 0.10 = 10 %) and for every group whose headers grew by more than that.
It never fails (exit 0) on a regression.
"""
import argparse, json, os, re, statistics, subprocess, sys, tarfile, tempfile, time, io

TUS = {
    "include_only": "#include <base256/base256.hpp>\nint main() { return 0; }\n",
    "arith": ("#include <base256/float256.hpp>\nusing base256::float256;\n"
              "float256 f(float256 a, float256 b, float256 c) { return a * b + c / a - sqrt(b); }\n"
              "int main() { return f(float256(1), float256(2), float256(3)) > float256(0); }\n"),
    "math": ("#include <base256/math.hpp>\nusing base256::float256;\n"
             "float256 f(float256 x) { return base256::exp(x) + base256::log(x) + base256::sin(x) + base256::pow(x, x) + base256::atan(x); }\n"
             "int main() { return f(float256(2)) > float256(0); }\n"),
    "special": ("#include <base256/base256.hpp>\nusing base256::float256;\n"
                "float256 f(float256 x) { return base256::cyl_bessel_j(x, x) + base256::airy_ai(-x) + base256::polygamma(2, -x)"
                " + base256::hurwitz_zeta(x, x) + base256::tgamma_upper(-x, x) + base256::riemann_zeta(x); }\n"
                "int main() { return f(float256(2.5)) > float256(0); }\n"),
    "igamma": ("#include <base256/special_math.hpp>\nusing base256::float256;\n"
               "float256 f(float256 a, float256 x) { return base256::tgamma_upper(a, x); }\n"
               "int main() { return f(float256(-2.5), float256(1)) > float256(0); }\n"),
    "complex": ("#include <base256/complex.hpp>\nusing base256::complex;\nusing base256::float256;\n"
                "complex f(complex z) { return base256::exp(z) * base256::log(z) + base256::sin(z) / base256::sqrt(z)"
                " + base256::pow(z, z) + base256::atanh(z); }\n"
                "int main() { return f(complex(float256(2), float256(1))).real() > float256(0); }\n"),
}
# TUs that need a header newer than some revisions: skipped (n/a) where it is missing.
REQUIRES = {"complex": "base256/complex.hpp"}

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

# Header groups: (top header, group it depends on). A group is the transitive
# include closure of its top header minus the closure of the group it builds on.
GROUPS = {"core": ("base256/base256.hpp", None), "complex": ("base256/complex.hpp", "core")}
TU_GROUP = {"complex": "complex"}  # every other TU: "core"

def closure(inc, top):
    """files reachable from inc/top through #include "..." and <base256/...>"""
    seen, stack = set(), [os.path.normpath(os.path.join(inc, top))]
    while stack:
        f = stack.pop()
        if f in seen or not os.path.exists(f): continue
        seen.add(f)
        for m in re.finditer(r'^\s*#\s*include\s*(?:"([^"]+)"|<(base256/[^>]+)>)', open(f).read(), re.M):
            stack.append(os.path.normpath(os.path.join(os.path.dirname(f), m.group(1)) if m.group(1)
                                          else os.path.join(inc, m.group(2))))
    return seen

def group_files(inc, g):
    top, dep = GROUPS[g]
    return closure(inc, top) - (closure(inc, GROUPS[dep][0]) if dep else set())

def header_bytes(inc):
    """{group: bytes, "all": bytes of every .hpp under base256/}"""
    tot = 0
    for d, _, fs in os.walk(os.path.join(inc, "base256")):
        tot += sum(os.path.getsize(os.path.join(d, f)) for f in fs if f.endswith(".hpp"))
    r = {g: sum(os.path.getsize(f) for f in group_files(inc, g)) for g in GROUPS}
    r["all"] = tot
    return r

def extract(rev, dest):
    data = subprocess.run(["git", "-C", ROOT, "archive", rev, "include"], check=True, capture_output=True).stdout
    tarfile.open(fileobj=io.BytesIO(data)).extractall(dest, filter="data")
    return os.path.join(dest, "include")

def compile_once(cxx, inc, src):
    """(wall, cpu) seconds; cpu = user + sys of the compiler and its children
    (cc1plus), which other load on the machine disturbs much less than wall time."""
    t = time.perf_counter()
    p = subprocess.Popen([cxx, "-std=c++23", "-O2", "-I", inc, "-c", src, "-o", os.devnull])
    _, status, ru = os.wait4(p.pid, 0)
    p.returncode = os.waitstatus_to_exitcode(status)
    if p.returncode != 0: raise SystemExit(f"{cxx} failed on {src}")
    return time.perf_counter() - t, ru.ru_utime + ru.ru_stime

def measure(cxxs, incs, n, tmp, want=None):
    """want(lab, tu): measure this TU for this revision (default: all)"""
    srcs = {}
    for name, code in TUS.items():
        p = os.path.join(tmp, name + ".cpp")
        open(p, "w").write(code); srcs[name] = p
    times = {(lab, cxx, tu): [] for lab in incs for cxx in cxxs for tu in TUS}
    have = {(lab, tu): (tu not in REQUIRES or os.path.exists(os.path.join(inc, REQUIRES[tu])))
            and (want is None or want(lab, tu)) for lab, inc in incs.items() for tu in TUS}
    for cxx in cxxs:  # warm-up (page cache, compiler binary)
        for lab, inc in incs.items(): compile_once(cxx, inc, srcs["include_only"])
    for _ in range(n):
        for cxx in cxxs:
            for tu, src in srcs.items():
                for lab, inc in incs.items():  # interleaved, so drift hits all revisions alike
                    if have[(lab, tu)]: times[(lab, cxx, tu)].append(compile_once(cxx, inc, src))
    res = {}
    for lab, inc in incs.items():
        r = {"header_bytes": header_bytes(inc), "n": n, "times": {}}
        for cxx in cxxs:
            for tu in TUS:
                if not have[(lab, tu)]: continue
                w = [t[0] for t in times[(lab, cxx, tu)]]
                v = [t[1] for t in times[(lab, cxx, tu)]]
                sd = lambda u: statistics.stdev(u) if len(u) > 1 else 0.0
                r["times"][f"{cxx} {tu}"] = {"mean": round(statistics.mean(v), 4), "sd": round(sd(v), 4),
                                             "wall": round(statistics.mean(w), 4), "wall_sd": round(sd(w), 4)}
        res[lab] = r
    return res

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--cxx", default="g++-14,clang++-19")
    ap.add_argument("-n", type=int, default=8)
    ap.add_argument("--ref", action="append", default=[])
    ap.add_argument("--json")
    ap.add_argument("--compare")
    ap.add_argument("--threshold", type=float, default=0.10)
    a = ap.parse_args()
    cxxs = a.cxx.split(",")
    base = json.load(open(a.compare)) if a.compare else None
    gref = {g: base["groups"][g]["ref"] for g in GROUPS} if base else {}
    refs = list(dict.fromkeys(list(a.ref) + list(gref.values())))
    # a baseline revision that is not also an explicit --ref only compiles its own group's TUs
    want = lambda lab, tu: (lab == "worktree" or lab in a.ref or gref.get(TU_GROUP.get(tu, "core")) == lab)
    with tempfile.TemporaryDirectory() as tmp:
        incs = {"worktree": os.path.join(ROOT, "include")}
        for r in refs: incs[r] = extract(r, os.path.join(tmp, "rev-" + r))
        res = measure(cxxs, incs, a.n, tmp, want)
    print(f"compile time, -std=c++23 -O2 -c, CPU seconds (user+sys), mean ± sd of {a.n}; wall-clock means in [brackets]")
    labs = list(res)
    print("| TU | " + " | ".join(labs) + " |"); print("|---" * (len(labs) + 1) + "|")
    for g in list(GROUPS) + ["all"]:
        name = f"{g} header bytes" if g != "all" else "all header bytes"
        print(f"| {name} | " + " | ".join(f"{res[l]['header_bytes'][g]:,}" if res[l]['header_bytes'][g] else "n/a" for l in labs) + " |")
    def cell(l, key):
        t = res[l]["times"].get(key)
        return "n/a" if t is None else f"{t['mean']:.3f} ± {t['sd']:.3f} [{t['wall']:.3f}]"
    for key in res["worktree"]["times"]:
        print(f"| {key} | " + " | ".join(cell(l, key) for l in labs) + " |")
    if a.json: json.dump(res, open(a.json, "w"), indent=1)
    if base:
        w = res["worktree"]
        warned = 0
        for g in GROUPS:
            ref = gref[g]; b = res[ref]
            print(f"\n{g} group vs {ref} ({base['groups'][g].get('version', '')}) on this machine, threshold +{a.threshold:.0%} (warn only)")
            for key, v in w["times"].items():
                if TU_GROUP.get(key.split(" ", 1)[1], "core") != g: continue
                if key not in b["times"]:
                    print(f"  {key:28s} {v['mean']:.3f} s  (new TU, not in {ref})")
                    continue
                r = v["mean"] / b["times"][key]["mean"]
                flag = r > 1 + a.threshold
                print(f"  {key:28s} {b['times'][key]['mean']:.3f} -> {v['mean']:.3f} s  ({r - 1:+.1%}){'  REGRESSION' if flag else ''}")
                if flag:
                    warned += 1
                    print(f"::warning title=compile-time regression ({g})::{key}: {b['times'][key]['mean']:.3f} s -> {v['mean']:.3f} s "
                          f"({r - 1:+.1%}, threshold +{a.threshold:.0%}, baseline {ref})")
            hb, hw = b["header_bytes"][g], w["header_bytes"][g]
            if not hb:
                print(f"  {g} header bytes {hw:,} (group not in {ref})")
                continue
            hr = hw / hb
            print(f"  {g} header bytes {hb:,} -> {hw:,} ({hr - 1:+.1%})")
            if hr > 1 + a.threshold:
                warned += 1
                print(f"::warning title=header size regression ({g})::{hb:,} -> {hw:,} bytes "
                      f"({hr - 1:+.1%}, threshold +{a.threshold:.0%}, baseline {ref})")
        print(f"{warned} warning(s)")

if __name__ == "__main__":
    main()
