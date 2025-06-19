#!/bin/bash

if [ $# -eq 0 ]; then
    echo "Usage: $0 <TARGET>"
    exit 1
fi

TARGET="$1"
ffmpeg -framerate 60 -i renders/${TARGET}_%d.png -c:v libx264 -pix_fmt yuv420p -crf 18 ${TARGET}.mp4