#!/usr/bin/env bash
#
# Initialize the permissive (MIT / BSD-3) external source submodules that the
# offline test suite needs. Nothing else is fetched.
#
#   ./scripts/init_permissive_sources.sh                  # test-required only
#   ./scripts/init_permissive_sources.sh --all-permissive # every permissive tier entry
#   ./scripts/init_permissive_sources.sh --check          # report state, fetch nothing
#
# This script will never initialize sources/copyleft/** or
# sources/license-unknown/**. Those are reference-only: GPL/AGPL or
# no-license-granted. They are pinned in the tree for provenance and marked
# `update = none` in .gitmodules so even `git submodule update --init` skips
# them. Fetching one has to be a deliberate act, not a side effect of setup.
#
# Exit codes: 0 ok, 1 usage/precondition failure, 2 a submodule did not land on
# its pinned commit.

set -euo pipefail

REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$REPO_ROOT"

# Submodules the offline pytest suite actually reads. Keep this list minimal:
# every entry here is a mandatory network dependency for a fresh clone.
REQUIRED=(
  "sources/permissive/ngrilli_Oxford_Crystal_Plasticity"
)

# Every permissive-tier entry. Safe to fetch (MIT / BSD-3) but not needed offline.
ALL_PERMISSIVE=(
  "sources/permissive/bibekanandadatta_Abaqus-UEL-Elasticity"
  "sources/permissive/bibekanandadatta_Abaqus-UEL-Hyperelasticity"
  "sources/permissive/jgomezc1_ABAQUS-US"
  "sources/permissive/ngrilli_Oxford_Crystal_Plasticity"
)

MODE="required"
case "${1:-}" in
  ""|--required-only) MODE="required" ;;
  --all-permissive)   MODE="all" ;;
  --check)            MODE="check" ;;
  -h|--help)          sed -n '2,20p' "$0" | sed 's/^# \{0,1\}//'; exit 0 ;;
  *) echo "unknown option: $1 (try --help)" >&2; exit 1 ;;
esac

if [ ! -f .gitmodules ]; then
  echo "ERROR: .gitmodules is missing. This repository pins external sources as" >&2
  echo "       submodules; without .gitmodules they cannot be initialized." >&2
  exit 1
fi

if [ "$MODE" = "required" ]; then TARGETS=("${REQUIRED[@]}"); else TARGETS=("${ALL_PERMISSIVE[@]}"); fi

# Guard against anyone editing the lists above to smuggle in a restricted tier.
for path in "${TARGETS[@]}"; do
  case "$path" in
    sources/permissive/*) ;;
    *) echo "REFUSING: $path is not under sources/permissive/." >&2
       echo "          Only permissive-tier sources may be initialized by this script." >&2
       exit 1 ;;
  esac
done

pinned_sha() { git ls-tree HEAD "$1" | awk '{print $3}'; }
actual_sha() { git -C "$1" rev-parse HEAD 2>/dev/null || echo "-"; }

if [ "$MODE" = "check" ]; then
  echo "External source submodule state (nothing fetched):"
  git submodule status
  echo
  echo "Restricted tiers (must stay empty unless deliberately fetched):"
  for p in sources/copyleft/* sources/license-unknown/*; do
    [ -d "$p" ] || continue
    git ls-tree HEAD "$p" | grep -q '^160000' || continue
    n=$(ls -A "$p" 2>/dev/null | wc -l)
    echo "  $p -> $n entries"
  done
  exit 0
fi

echo "Initializing permissive external sources (mode: $MODE)"
for path in "${TARGETS[@]}"; do
  want="$(pinned_sha "$path")"
  if [ -z "$want" ]; then
    echo "ERROR: $path is not a pinned gitlink in HEAD." >&2
    exit 1
  fi
  if [ "$(actual_sha "$path")" = "$want" ]; then
    echo "  ok (already at pinned commit)  $path"
    continue
  fi
  echo "  fetching  $path"
  git submodule update --init -- "$path"
  got="$(actual_sha "$path")"
  if [ "$got" != "$want" ]; then
    echo "ERROR: $path landed on $got but the tree pins $want." >&2
    echo "       Refusing to continue with an unpinned external source." >&2
    exit 2
  fi
done

echo
echo "Permissive sources ready. Restricted tiers were not touched:"
for p in sources/copyleft/* sources/license-unknown/*; do
  [ -d "$p" ] || continue
  git ls-tree HEAD "$p" | grep -q '^160000' || continue
  n=$(ls -A "$p" 2>/dev/null | wc -l)
  echo "  $p -> $n entries (expected 0)"
done
