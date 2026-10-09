#!/bin/sh
# Build the TaiChi reduced-capacity collision attack.
# Picks sensible flags for Apple clang, Homebrew gcc, or Linux g++/clang.
set -e
SRC=taichi_real_attack.cpp
OUT=taichi_real
STD="-O3 -std=c++17 -pthread"

uname_s=$(uname -s)
if [ "$uname_s" = "Darwin" ]; then
  # Apple clang does not accept -march=native; -mcpu=native works on recent clang,
  # otherwise fall back to a plain build (still fully vectorised by the compiler).
  if clang++ $STD -mcpu=native -x c++ -c -o /dev/null - </dev/null 2>/dev/null; then
    CXX="${CXX:-clang++}"; ARCH="-mcpu=native"
  else
    CXX="${CXX:-clang++}"; ARCH=""
  fi
else
  CXX="${CXX:-g++}"; ARCH="-march=native -mprefer-vector-width=512"
fi

echo "building with: $CXX $STD $ARCH"
$CXX $STD $ARCH "$SRC" -o "$OUT"
echo "built ./$OUT"
echo "reference ground truth (optional): cc -O2 -std=c11 ref.c -o ref"
