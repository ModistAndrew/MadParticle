#!/bin/bash
TARGET="tower"
BLEND="${TARGET}.blend"
BLENDER_PATH="$HOME/zjx/blender/blender"

for f in $(ls -v ../output/${TARGET}_*.obj); do
    [ -e "$f" ] || continue
    base=$(basename "$f" .obj)
    frame="${base##*_}"
    out="../renders/${TARGET}_${frame}.png"
    $BLENDER_PATH -b "$BLEND" --python batch_render.py -- "$f" "$out"
done