#!/usr/bin/env python3
"""Offline reproducibility check for the pinned external source submodules.

Runs with no network access. Verifies that what a fresh clone would get is
self-consistent and that the license-tier policy is actually enforced by
configuration rather than only by documentation:

  1. .gitmodules exists.
  2. Every gitlink recorded in HEAD has a .gitmodules mapping (path + url).
  3. Every .gitmodules entry corresponds to a real gitlink (no dead mappings).
  4. copyleft/ and license-unknown/ entries are `update = none`, so a default
     `git submodule update --init` cannot fetch restricted source.
  5. Permissive entries are NOT `update = none` (they must stay fetchable).
  6. Pinned SHAs match the values this project verified upstream.
  7. Any populated submodule sits exactly on its pinned commit.
  8. Restricted tiers are empty on disk.
  9. A populated permissive submodule still carries its upstream license file.

Exit code 0 = all checks pass, 1 = at least one failed.

    python scripts/verify_source_submodules.py [-v]
"""

import os
import subprocess
import sys

REPO_ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), os.pardir))

# SHAs confirmed present in the named upstream repositories. Changing a value
# here is a deliberate re-pin and must be justified in the commit message.
EXPECTED_PINS = {
    "sources/permissive/ngrilli_Oxford_Crystal_Plasticity":
        "85102ed35dc4592edc5fd4aaec542fa63e7162fd",
    "sources/permissive/bibekanandadatta_Abaqus-UEL-Elasticity":
        "9187e54a4069ed0e6196a29e97e92185d14724a5",
    "sources/permissive/bibekanandadatta_Abaqus-UEL-Hyperelasticity":
        "ba5bf018b71e0717d3a3136a351637a27e6ff0ed",
    "sources/permissive/jgomezc1_ABAQUS-US":
        "54181407aa7aa23055e33d354d0b2a3abc266365",
    "sources/copyleft/ICAMS_Crystal_Plasticity_UMAT":
        "469464b646a6f250d996828559342b81a1dfbbd7",
    "sources/license-unknown/TarletonGroup_CrystalPlasticity":
        "2f0909472f4cdf1b0b71da6ba97db900dacc6f05",
}

RESTRICTED_TIERS = ("copyleft", "license-unknown")
LICENSE_NAMES = ("LICENSE", "LICENSE.txt", "LICENSE.md", "COPYING", "COPYING.txt")


def git(*args):
    return subprocess.run(["git", "-C", REPO_ROOT] + list(args),
                          capture_output=True, text=True).stdout.strip()


def gitlinks():
    """path -> pinned sha, for every gitlink in HEAD."""
    out = {}
    for line in git("ls-tree", "-r", "HEAD").splitlines():
        meta, _, path = line.partition("\t")
        parts = meta.split()
        if len(parts) >= 3 and parts[1] == "commit":
            out[path] = parts[2]
    return out


def gitmodules_entries():
    """name -> {path, url, update} parsed from .gitmodules via git config."""
    path = os.path.join(REPO_ROOT, ".gitmodules")
    if not os.path.exists(path):
        return None
    raw = git("config", "-f", ".gitmodules", "--list")
    entries = {}
    for line in raw.splitlines():
        key, _, value = line.partition("=")
        if not key.startswith("submodule."):
            continue
        rest = key[len("submodule."):]
        name, _, field = rest.rpartition(".")
        entries.setdefault(name, {})[field] = value
    return entries


def tier_of(path):
    parts = path.split("/")
    return parts[1] if len(parts) > 2 and parts[0] == "sources" else None


def main():
    verbose = "-v" in sys.argv
    failures = []
    notes = []

    def check(ok, label, detail=""):
        (notes if ok else failures).append((label, detail))
        if verbose or not ok:
            print("  %s %s%s" % ("PASS" if ok else "FAIL", label,
                                 ("  -- " + detail) if detail else ""))

    print("External source submodule reproducibility check")
    print("repo: %s\n" % REPO_ROOT)

    entries = gitmodules_entries()
    if entries is None:
        print("  FAIL .gitmodules is missing -- a fresh clone cannot initialize "
              "any pinned external source.")
        return 1
    check(True, ".gitmodules present")

    links = gitlinks()
    check(bool(links), "HEAD records gitlinks", "%d found" % len(links))

    mapped_paths = {e.get("path") for e in entries.values()}

    # 2 + 3: mapping completeness, both directions.
    for path, sha in sorted(links.items()):
        check(path in mapped_paths, "mapped in .gitmodules: %s" % path)
    for name, e in sorted(entries.items()):
        p = e.get("path")
        check(p in links, "mapping points at a real gitlink: %s" % name,
              "" if p in links else "path %r is not a gitlink in HEAD" % p)
        check(bool(e.get("url")), "mapping has a url: %s" % name)

    # 4 + 5: the tier policy must be enforced by config.
    for name, e in sorted(entries.items()):
        p = e.get("path") or ""
        tier = tier_of(p)
        upd = e.get("update")
        if tier in RESTRICTED_TIERS:
            check(upd == "none",
                  "restricted tier is update=none: %s" % p,
                  "" if upd == "none" else
                  "update=%r would let a default init fetch restricted source" % upd)
        elif tier == "permissive":
            check(upd != "none", "permissive tier stays fetchable: %s" % p)

    # 6: pinned SHAs unchanged.
    for path, expected in sorted(EXPECTED_PINS.items()):
        actual = links.get(path)
        check(actual == expected, "pinned SHA preserved: %s" % path,
              "" if actual == expected else "expected %s got %s" % (expected, actual))

    # 7 + 8 + 9: on-disk state.
    for path, sha in sorted(links.items()):
        abs_path = os.path.join(REPO_ROOT, path)
        populated = os.path.isdir(abs_path) and bool(os.listdir(abs_path))
        tier = tier_of(path)

        if tier in RESTRICTED_TIERS:
            check(not populated, "restricted tier not vendored: %s" % path,
                  "" if not populated else "directory has content; it must stay empty")
            continue

        if not populated:
            if verbose:
                print("  .... not initialized (optional): %s" % path)
            continue

        head = subprocess.run(["git", "-C", abs_path, "rev-parse", "HEAD"],
                              capture_output=True, text=True).stdout.strip()
        check(head == sha, "populated submodule on pinned commit: %s" % path,
              "" if head == sha else "worktree at %s but tree pins %s" % (head, sha))
        has_license = any(os.path.exists(os.path.join(abs_path, n)) for n in LICENSE_NAMES)
        check(has_license, "upstream license retained: %s" % path,
              "" if has_license else "no LICENSE/COPYING file found")

    print("\n%d checks passed, %d failed" % (len(notes), len(failures)))
    if failures:
        print("\nFailures:")
        for label, detail in failures:
            print("  - %s%s" % (label, ("  -- " + detail) if detail else ""))
        return 1
    print("Fresh-clone reproducibility contract holds.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
