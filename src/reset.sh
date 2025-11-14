#!/bin/bash

rm ./out*
rm ./slurm*
mv ./fields/*.001 ./
mv ./fields/*.000 ./
rm ./fields/*
mv ./w.0* ./fields/
mv ./p.0* ./fields/
mv ./files/seed.000 ./
rm ./files/*
mv ./seed.000 ./files/
rm ./curframe.dat
cp ../src/comm/curframe.dat curframe.dat

