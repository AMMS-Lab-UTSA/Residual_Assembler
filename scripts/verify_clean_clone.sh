#!/usr/bin/env bash
#
# Fresh-clone reproducibility procedure.
#
# Clones this repository into a throwaway directory exactly as a new
# collaborator would, then proves the documented setup path actually works:
#
#   1. `git submodule status` succeeds (it used to die on a missing mapping).
#   2. The offline verifier passes on the clone.
#   3. A bare `git submodule update --init` does NOT fetch copyleft or
#      license-unknown source.
#   4. `./scripts/init_permissive_sources.sh` lands the test dependency on its
#      pinned commit.
#   5. `pytest -q` reports the expected counts.
#   6. The parent worktree is still clean afterwards.
#
# Usage:
#   ./scripts/verify_clean_clone.sh [<python>] [<source-repo>] [<ref>]
#
# Defaults to the shared venv interpreter, to this checkout as the clone source
# (so it can be run offline against local objects), and to whatever ref this
# working copy currently has checked out. Pass the GitHub URL to verify what
# collaborators actually receive, and a ref to check a specific branch or tag.
#
# The ref is deliberately NOT hard-coded to a feature branch: this procedure has
# to keep working after that branch is merged or released.
#
# Environment:
#   EXPECTED_PYTEST   baseline line pytest must report (default "68 passed, 9 skipped")

set -euo pipefail

REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
PYTHON="${1:-$HOME/softwarex_work/.venv/bin/python}"
SOURCE="${2:-$REPO_ROOT}"
if [ -n "${3:-}" ]; then
  REF="$3"
else
  REF="$(git -C "$REPO_ROOT" symbolic-ref --quiet --short HEAD 2>/dev/null || true)"
  if [ -z "$REF" ]; then
    # Detached HEAD: fall back to the exact commit so the check still means something.
    REF="$(git -C "$REPO_ROOT" rev-parse HEAD)"
  fi
fi
EXPECTED_PYTEST="${EXPECTED_PYTEST:-68 passed, 9 skipped}"

WORK="$(mktemp -d)"
cleanup() { rm -rf "$WORK"; }
trap cleanup EXIT

fail() { echo; echo "FAILED: $*" >&2; exit 1; }

echo "clean-clone verification"
echo "  source : $SOURCE"
echo "  ref    : $REF"
echo "  python : $PYTHON"
echo "  workdir: $WORK"
echo

echo "[1/7] cloning"
if ! git clone -q --branch "$REF" "$SOURCE" "$WORK/clone" 2>/dev/null; then
  # Not a branch or tag name (e.g. a raw commit): clone then check it out.
  git clone -q "$SOURCE" "$WORK/clone"
  git -C "$WORK/clone" checkout -q --detach "$REF"
fi
cd "$WORK/clone"
echo "      at $(git rev-parse --short HEAD)"

echo "[2/7] git submodule status"
if ! git submodule status > "$WORK/status.txt" 2>"$WORK/status.err"; then
  cat "$WORK/status.err" >&2
  fail "git submodule status returned non-zero on a fresh clone"
fi
sed 's/^/      /' "$WORK/status.txt"

echo "[3/7] offline verifier"
"$PYTHON" scripts/verify_source_submodules.py || fail "verify_source_submodules.py"

echo "[4/7] bare 'git submodule update --init' must skip restricted tiers"
git submodule update --init > "$WORK/init.log" 2>&1 || true
for p in sources/copyleft/* sources/license-unknown/*; do
  [ -d "$p" ] || continue
  git ls-tree HEAD "$p" | grep -q '^160000' || continue
  n=$(ls -A "$p" 2>/dev/null | wc -l)
  [ "$n" -eq 0 ] || fail "restricted source was fetched into $p ($n entries)"
  echo "      empty as required: $p"
done

echo "[5/7] bootstrap permissive sources"
./scripts/init_permissive_sources.sh > "$WORK/boot.log" 2>&1 \
  || { cat "$WORK/boot.log"; fail "init_permissive_sources.sh"; }
sub="sources/permissive/ngrilli_Oxford_Crystal_Plasticity"
want="$(git ls-tree HEAD "$sub" | awk '{print $3}')"
got="$(git -C "$sub" rev-parse HEAD)"
[ "$want" = "$got" ] || fail "$sub at $got, expected pinned $want"
echo "      $sub @ $got"

echo "[6/7] strict mode must refuse to run green without the sources"
git submodule deinit -qf "$sub" >/dev/null 2>&1
set +e
strict_out="$(REQUIRE_EXTERNAL_TEST_SOURCES=1 "$PYTHON" -m pytest -q 2>&1 | tail -2)"
set -e
echo "$strict_out" | sed 's/^/      /'
echo "$strict_out" | grep -qE '[0-9]+ (failed|error)' \
  || fail "strict mode stayed green with the required external source removed"
./scripts/init_permissive_sources.sh > "$WORK/reboot.log" 2>&1 \
  || { cat "$WORK/reboot.log"; fail "re-bootstrap after strict-mode check"; }

echo "[7/7] pytest and final cleanliness"
set +e
out="$(REQUIRE_EXTERNAL_TEST_SOURCES=1 "$PYTHON" -m pytest -q 2>&1 | tail -3)"
set -e
echo "$out" | sed 's/^/      /'
echo "$out" | grep -q "$EXPECTED_PYTEST" \
  || fail "expected '$EXPECTED_PYTEST' from pytest on a fresh clone"

dirty="$(git status --porcelain)"
[ -z "$dirty" ] || { echo "$dirty"; fail "parent worktree is dirty after setup"; }
echo "      parent worktree clean"

echo
echo "PASS: fresh clone initializes, tests at '$EXPECTED_PYTEST', worktree clean."
