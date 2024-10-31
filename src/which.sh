#!/bin/bash
make clean
rm prof.nvvp
make

echo " "
echo " "
echo " "
echo " "
# echo "UNCOMMENTED"
# echo " "

nvprof --unified-memory-profiling per-process-device --track-memory-allocations on -o prof.nvvp ./main.exe
# nvprof --print-gpu-trace ./main.exe
# ./main.exe

# rm Frames/wgpuk*

echo " "
echo " "
echo " "
echo " "