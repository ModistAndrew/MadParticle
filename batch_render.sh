#!/bin/bash

if [ $# -eq 0 ]; then
    echo "Usage: $0 <TARGET>"
    exit 1
fi

TARGET="$1"
BLEND="assets/${TARGET}.blend"
BLENDER_PATH="$HOME/zjx/blender/blender"

for f in $(ls -v output/${TARGET}_*.obj); do
    [ -e "$f" ] || continue
    base=$(basename "$f" .obj)
    frame="${base##*_}"
    out="renders/${TARGET}_${frame}.png"
    $BLENDER_PATH -b "$BLEND" --python batch_render.py -- "$f" "$out"
done