#!/bin/bash
# Default-path identity: build tools/mode_digest.cpp against the headers of REV
# and of the working tree (same compiler and flags) and require identical output.
# usage: tools/mode_digest.sh REV [CXX=g++-14] [arith_cases] [fn_args]
set -euo pipefail
rev=$1; cxx=${2:-g++-14}; n=${3:-200000}; k=${4:-24}
root=$(cd "$(dirname "$0")/.." && pwd)
tmp=$(mktemp -d); trap 'rm -rf "$tmp"' EXIT
git -C "$root" archive "$rev" include | tar -x -C "$tmp"
mkdir -p "$tmp/tests"; cp "$root/tests/fn_lists.hpp" "$tmp/tests/"
"$cxx" -std=c++23 -O2 -I "$tmp/include" "$root/tools/mode_digest.cpp" -o "$tmp/old"
"$cxx" -std=c++23 -O2 -I "$root/include" "$root/tools/mode_digest.cpp" -o "$tmp/new"
"$tmp/old" "$n" "$k" > "$tmp/old.txt" & "$tmp/new" "$n" "$k" > "$tmp/new.txt"; wait
lines=$(wc -l < "$tmp/new.txt")
if cmp -s "$tmp/old.txt" "$tmp/new.txt"; then
    echo "identical: $lines digests, $rev vs working tree, $cxx, arith_cases=$n fn_args=$k"
else
    diff "$tmp/old.txt" "$tmp/new.txt" || true
    echo "DIFFERENT: $rev vs working tree"; exit 1
fi
