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
#   ./scripts/verify_clean_clone.sh [<python>] [<source-repo>]
#
# Defaults to the shared venv interpreter and to this checkout as the clone
# source (so it can be run offline against local objects). Pass the GitHub URL
# to verify what collaborators actually receive.

set -euo pipefail

REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
PYTHON="${1:-$HOME/softwarex_work/.venv/bin/python}"
SOURCE="${2:-$REPO_ROOT}"
BRANCH="feature/umat-oti-residual-bridge"
EXPECTED_PYTEST="68 passed, 9 skipped"

WORK="$(mktemp -d)"
cleanup() { rm -rf "$WORK"; }
trap cleanup EXIT

fail() { echo; echo "FAILED: $*" >&2; exit 1; }

echo "clean-clone verification"
echo "  source : $SOURCE"
echo "  branch : $BRANCH"
echo "  python : $PYTHON"
echo "  workdir: $WORK"
echo

echo "[1/6] cloning"
git clone -q --branch "$BRANCH" "$SOURCE" "$WORK/clone"
cd "$WORK/clone"

echo "[2/6] git submodule status"
if ! git submodule status > "$WORK/status.txt" 2>"$WORK/status.err"; then
  cat "$WORK/status.err" >&2
  fail "git submodule status returned non-zero on a fresh clone"
fi
sed 's/^/      /' "$WORK/status.txt"

echo "[3/6] offline verifier"
"$PYTHON" scripts/verify_source_submodules.py || fail "verify_source_submodules.py"

echo "[4/6] bare 'git submodule update --init' must skip restricted tiers"
git submodule update --init > "$WORK/init.log" 2>&1 || true
for p in sources/copyleft/* sources/license-unknown/*; do
  [ -d "$p" ] || continue
  git ls-tree HEAD "$p" | grep -q '^160000' || continue
  n=$(ls -A "$p" 2>/dev/null | wc -l)
  [ "$n" -eq 0 ] || fail "restricted source was fetched into $p ($n entries)"
  echo "      empty as required: $p"
done

echo "[5/6] bootstrap permissive sources"
./scripts/init_permissive_sources.sh > "$WORK/boot.log" 2>&1 \
  || { cat "$WORK/boot.log"; fail "init_permissive_sources.sh"; }
sub="sources/permissive/ngrilli_Oxford_Crystal_Plasticity"
want="$(git ls-tree HEAD "$sub" | awk '{print $3}')"
got="$(git -C "$sub" rev-parse HEAD)"
[ "$want" = "$got" ] || fail "$sub at $got, expected pinned $want"
echo "      $sub @ $got"

echo "[6/6] pytest and final cleanliness"
set +e
out="$("$PYTHON" -m pytest -q 2>&1 | tail -3)"
set -e
echo "$out" | sed 's/^/      /'
echo "$out" | grep -q "$EXPECTED_PYTEST" \
  || fail "expected '$EXPECTED_PYTEST' from pytest on a fresh clone"

dirty="$(git status --porcelain)"
[ -z "$dirty" ] || { echo "$dirty"; fail "parent worktree is dirty after setup"; }
echo "      parent worktree clean"

echo
echo "PASS: fresh clone initializes, tests at '$EXPECTED_PYTEST', worktree clean."
