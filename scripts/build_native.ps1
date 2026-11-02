param(
    [string]$RtoolsHome = 'D:\Program Files\rtools45',
    [string]$ArmadilloInclude = ''
)

$ErrorActionPreference = 'Stop'
$root = Split-Path -Parent $PSScriptRoot
$source = Join-Path $root 'src\sctransform_cpu\native'
$toolBin = Join-Path $RtoolsHome 'x86_64-w64-mingw32.static.posix\bin'
$toolUsr = Join-Path $RtoolsHome 'usr\bin'
$gcc = Join-Path $toolBin 'gcc.exe'
$gxx = Join-Path $toolBin 'g++.exe'
$gfortran = Join-Path $toolBin 'gfortran.exe'
$strip = Join-Path $toolBin 'strip.exe'
$build = Join-Path $root 'build\native'
$output = Join-Path $source 'sct_port_r45.dll'
$dependencyRoot = Join-Path $root 'build\native-deps'
if (-not $ArmadilloInclude) {
    $ArmadilloInclude = Join-Path $dependencyRoot 'RcppArmadillo\inst\include'
}
if (-not (Test-Path -LiteralPath $ArmadilloInclude)) {
    New-Item -ItemType Directory -Force -Path $dependencyRoot | Out-Null
    $archive = Join-Path $dependencyRoot 'RcppArmadillo_15.2.3-1.tar.gz'
    Invoke-WebRequest `
        'https://cran.r-project.org/src/contrib/Archive/RcppArmadillo/RcppArmadillo_15.2.3-1.tar.gz' `
        -OutFile $archive
    $hash = (Get-FileHash -Algorithm SHA256 -LiteralPath $archive).Hash
    if ($hash -ne '5BBB4C32001439E4D2090EDAFF3E0B5A00249AAA41F89C94A5C3C1A6548656F3') {
        throw 'RcppArmadillo source hash mismatch'
    }
    tar -xf $archive -C $dependencyRoot
}

foreach ($required in @($gcc, $gfortran, $strip, $ArmadilloInclude, (Join-Path $source 'portsrc.f'))) {
    if (-not (Test-Path -LiteralPath $required)) {
        throw "required build input is missing: $required"
    }
}
New-Item -ItemType Directory -Force -Path $build | Out-Null
$oldPath = $env:PATH
try {
    $env:PATH = "$toolBin;$toolUsr;$oldPath"
    $cflags = @('-O2', '-Wall', '-std=gnu2x', '-mfpmath=sse', '-msse2', '-mstackrealign')
    & $gfortran -O2 -mfpmath=sse -msse2 -mstackrealign -c `
        (Join-Path $source 'portsrc.f') -o (Join-Path $build 'portsrc.o')
    if ($LASTEXITCODE) { throw 'portsrc.f compilation failed' }
    foreach ($name in @('port_standalone.c', 'r45_mathlib_subset.c', 'glmgampoi_compat.c')) {
        & $gcc @cflags -c (Join-Path $source $name) `
            -o (Join-Path $build ($name -replace '\.c$', '.o'))
        if ($LASTEXITCODE) { throw "$name compilation failed" }
    }
    & $gxx -O2 -Wall -std=gnu++17 -mfpmath=sse -msse2 -mstackrealign -c `
        (Join-Path $source 'rcpp_compat.cpp') -o (Join-Path $build 'rcpp_compat.o')
    if ($LASTEXITCODE) { throw 'rcpp_compat.cpp compilation failed' }
    & $gxx -O2 -Wall -std=gnu++17 -mfpmath=sse -msse2 -mstackrealign `
        "-I$ArmadilloInclude" -c `
        (Join-Path $source 'generic_design.cpp') -o (Join-Path $build 'generic_design.o')
    if ($LASTEXITCODE) { throw 'generic_design.cpp compilation failed' }
    & $gfortran -shared -o $output `
        (Join-Path $build 'portsrc.o') `
        (Join-Path $build 'port_standalone.o') `
        (Join-Path $build 'r45_mathlib_subset.o') `
        (Join-Path $build 'glmgampoi_compat.o') `
        (Join-Path $build 'rcpp_compat.o') `
        (Join-Path $build 'generic_design.o') `
        -static-libgfortran -static-libgcc -static-libstdc++ -lstdc++ `
        -llapack -lblas
    if ($LASTEXITCODE) { throw 'link failed' }
    & $strip --strip-unneeded $output
    if ($LASTEXITCODE) { throw 'strip failed' }
} finally {
    $env:PATH = $oldPath
}

Get-FileHash -Algorithm SHA256 -LiteralPath $output
