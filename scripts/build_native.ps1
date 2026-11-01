param(
    [string]$RtoolsHome = 'D:\Program Files\rtools45'
)

$ErrorActionPreference = 'Stop'
$root = Split-Path -Parent $PSScriptRoot
$source = Join-Path $root 'src\sctransform_cpu\native'
$toolBin = Join-Path $RtoolsHome 'x86_64-w64-mingw32.static.posix\bin'
$toolUsr = Join-Path $RtoolsHome 'usr\bin'
$gcc = Join-Path $toolBin 'gcc.exe'
$gfortran = Join-Path $toolBin 'gfortran.exe'
$build = Join-Path $root 'build\native'
$output = Join-Path $source 'sct_port_r45.dll'

foreach ($required in @($gcc, $gfortran, (Join-Path $source 'portsrc.f'))) {
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
    & $gfortran -shared -o $output `
        (Join-Path $build 'portsrc.o') `
        (Join-Path $build 'port_standalone.o') `
        (Join-Path $build 'r45_mathlib_subset.o') `
        (Join-Path $build 'glmgampoi_compat.o') `
        -static-libgfortran -static-libgcc
    if ($LASTEXITCODE) { throw 'link failed' }
} finally {
    $env:PATH = $oldPath
}

Get-FileHash -Algorithm SHA256 -LiteralPath $output
