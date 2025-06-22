#!/bin/bash

if [ $# -lt 1 ]; then
    echo "Usage: $0 <TARGET> [FRAME_NUMBER]"
    exit 1
fi

TARGET="$1"
FRAME="$2"
BLEND="assets/${TARGET}.blend"
BLENDER_PATH="$HOME/zjx/blender/blender"

for f in $(ls -v output/${TARGET}_*.obj); do
    [ -e "$f" ] || continue
    base=$(basename "$f" .obj)
    frame="${base##*_}"

    if [ -n "$FRAME" ] && [ "$frame" -ne "$FRAME" ]; then
        continue
    fi

    out="renders_spec/${TARGET}_${frame}.png"
    $BLENDER_PATH -b "$BLEND" --python batch_render.py -- "$f" "$out"
done