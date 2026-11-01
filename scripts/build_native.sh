#!/usr/bin/env bash
set -euo pipefail

root="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
source_dir="$root/src/sctransform_cpu/native"
build_dir="${SCT_NATIVE_BUILD_DIR:-$root/build/native-linux}"
output="$source_dir/libsct_port_r45.so"

if [[ "$(uname -s)" != "Linux" || "$(uname -m)" != "x86_64" ]]; then
    echo "the pinned native build currently supports Linux x86_64 only" >&2
    exit 2
fi

cc="${CC:-gcc}"
cxx="${CXX:-g++}"
fc="${FC:-gfortran}"
for tool in "$cc" "$cxx" "$fc" sha256sum; do
    command -v "$tool" >/dev/null || {
        echo "required build tool is unavailable: $tool" >&2
        exit 2
    }
done

mkdir -p "$build_dir"
common=(-O2 -Wall -fPIC -mfpmath=sse -msse2 -mstackrealign)

"$fc" "${common[@]}" -c "$source_dir/portsrc.f" \
    -o "$build_dir/portsrc.o"
for name in port_standalone r45_mathlib_subset glmgampoi_compat; do
    "$cc" "${common[@]}" -std=gnu2x -c "$source_dir/$name.c" \
        -o "$build_dir/$name.o"
done
"$cxx" "${common[@]}" -std=gnu++17 -c "$source_dir/rcpp_compat.cpp" \
    -o "$build_dir/rcpp_compat.o"

"$fc" -shared -o "$output" \
    "$build_dir/portsrc.o" \
    "$build_dir/port_standalone.o" \
    "$build_dir/r45_mathlib_subset.o" \
    "$build_dir/glmgampoi_compat.o" \
    "$build_dir/rcpp_compat.o" \
    -static-libgfortran -static-libgcc -lm \
    -Wl,-z,relro,-z,now

file "$output"
sha256sum "$output"
