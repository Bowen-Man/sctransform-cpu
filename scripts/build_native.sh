#!/usr/bin/env bash
set -euo pipefail

root="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
source_dir="$root/src/sctransform_cpu/native"
build_dir="${SCT_NATIVE_BUILD_DIR:-$root/build/native-linux}"
output="$source_dir/libsct_port_r45.so"
dependency_dir="$root/build/native-deps"
armadillo_include="${SCT_ARMADILLO_INCLUDE:-$dependency_dir/RcppArmadillo/inst/include}"

if [[ "$(uname -s)" != "Linux" || "$(uname -m)" != "x86_64" ]]; then
    echo "the pinned native build currently supports Linux x86_64 only" >&2
    exit 2
fi

cc="${CC:-gcc}"
cxx="${CXX:-g++}"
fc="${FC:-gfortran}"
for tool in "$cc" "$cxx" "$fc" curl tar strip sha256sum; do
    command -v "$tool" >/dev/null || {
        echo "required build tool is unavailable: $tool" >&2
        exit 2
    }
done

if [[ ! -d "$armadillo_include" ]]; then
    mkdir -p "$dependency_dir"
    archive="$dependency_dir/RcppArmadillo_15.2.3-1.tar.gz"
    curl -fL --retry 3 \
        https://cran.r-project.org/src/contrib/Archive/RcppArmadillo/RcppArmadillo_15.2.3-1.tar.gz \
        -o "$archive"
    echo "5bbb4c32001439e4d2090edaff3e0b5a00249aaa41f89c94a5c3c1a6548656f3  $archive" \
        | sha256sum -c -
    tar -xf "$archive" -C "$dependency_dir"
fi

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
"$cxx" "${common[@]}" -std=gnu++17 -I"$armadillo_include" \
    -c "$source_dir/generic_design.cpp" \
    -o "$build_dir/generic_design.o"

"$fc" -shared -o "$output" \
    "$build_dir/portsrc.o" \
    "$build_dir/port_standalone.o" \
    "$build_dir/r45_mathlib_subset.o" \
    "$build_dir/glmgampoi_compat.o" \
    "$build_dir/rcpp_compat.o" \
    "$build_dir/generic_design.o" \
    -static-libgcc \
    -Wl,-Bstatic -lstdc++ -llapack -lblas -Wl,-Bdynamic -lm \
    -Wl,-z,relro,-z,now

strip --strip-unneeded "$output"

file "$output"
sha256sum "$output"
