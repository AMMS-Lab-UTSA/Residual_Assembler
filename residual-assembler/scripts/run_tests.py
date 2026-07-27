#!/usr/bin/env python3
"""Run the standalone test scripts and summarize (offline by default).

Each tests/**/test_*.py is a self-contained script whose exit code is its result
(0 = pass). Tests that need OTILib or a real Abaqus run fail here with a
recognizable dependency message; this runner marks those SKIP rather than FAIL
unless you opt in with RESASM_RUN_OTILIB=1 / RESASM_RUN_ABAQUS=1.

    python scripts/run_tests.py                # everything; dep-gated -> SKIP
    python scripts/run_tests.py --core         # fast offline smoke set only
    python scripts/run_tests.py <file> [...]   # specific test file(s)
"""
import glob
import os
import subprocess
import sys
import time

RA = os.path.abspath(os.path.join(os.path.dirname(__file__), os.pardir))

# a fast, fully-offline smoke set (no OTILib, no Abaqus)
CORE = [
    "tests/framework/test_solid3d_kernel.py",
    "tests/framework/test_job_validation.py",
    "tests/framework/test_interface_versions.py",
]

# substrings that mean "this failed only because an external dependency or an
# optional data fixture is absent" (not a real code failure)
DEP_MARKERS = ("otilib", "OTILib", "libotilib", "No module named 'oti", "No module named 'pyoti",
               ".odb", "gfortran", "ABA_PARAM", "cannot find -l",
               "OTI library", "libpython", "Abaqus command",
               "fixture mesh not vendored", "SKIP:")
RUN_OTILIB = os.environ.get("RESASM_RUN_OTILIB") == "1"
RUN_ABAQUS = os.environ.get("RESASM_RUN_ABAQUS") == "1"


def classify(rc, out):
    if rc == 0:
        return "PASS"
    if not (RUN_OTILIB and RUN_ABAQUS) and any(m in out for m in DEP_MARKERS):
        return "SKIP"
    return "FAIL"


def main(argv):
    core = "--core" in argv
    files = [a for a in argv if not a.startswith("--")]
    if files:
        tests = [f if os.path.isabs(f) else os.path.join(RA, f) for f in files]
    elif core:
        tests = [os.path.join(RA, f) for f in CORE]
    else:
        tests = sorted(glob.glob(os.path.join(RA, "tests", "**", "test_*.py"), recursive=True))

    results = {}
    print("running %d test script(s)%s\n" % (len(tests), "  (--core)" if core else ""))
    for t in tests:
        rel = os.path.relpath(t, RA)
        t0 = time.time()
        try:
            p = subprocess.run([sys.executable, t], cwd=RA, capture_output=True,
                               text=True, timeout=420)
            out = p.stdout + p.stderr
            cls = classify(p.returncode, out)
        except subprocess.TimeoutExpired:
            cls = "FAIL"
        results[rel] = cls
        print("  %-5s  %-54s %5.1fs" % (cls, rel, time.time() - t0))

    n = {"PASS": 0, "FAIL": 0, "SKIP": 0}
    for c in results.values():
        n[c] += 1
    print("=" * 74)
    print("  %d passed, %d skipped (missing deps), %d failed  of %d"
          % (n["PASS"], n["SKIP"], n["FAIL"], len(results)))
    if n["FAIL"]:
        print("  FAILED: " + ", ".join(k for k, c in results.items() if c == "FAIL"))
    return 1 if n["FAIL"] else 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
