#!/bin/bash

# rm ./out*
# rm ./slurm*
# rm ./fields/*
# mv ./files/seed.000 ./seed.000
# rm ./files/*
# mv ./seed.000 ./files/seed.000
# rm ./curframe.dat
# cp ../src/startframe.dat curframe.dat

rm -r 2run
rm -r 3run
rm -r 4run
rm -r 5run
rm -r 6run
rm -r 7run
rm -r 8run
rm -r 9run

cp -vr 1run/ 2run
cp -vr 1run/ 3run
cp -vr 1run/ 4run
cp -vr 1run/ 5run
cp -vr 1run/ 6run
cp -vr 1run/ 7run
cp -vr 1run/ 8run
cp -vr 1run/ 9run

source ~/my_venv/bin/activate

echo "sourced python"

echo "1run"
cd ./1run/
python3 gen_noise.py
sbatch jobscript

echo "2run"
cd ../2run/
python3 gen_noise.py
sbatch jobscript

echo "3run"
cd ../3run/
python3 gen_noise.py
sbatch jobscript

echo "4run"
cd ../4run/
python3 gen_noise.py
sbatch jobscript

echo "5run"
cd ../5run/
python3 gen_noise.py
sbatch jobscript

echo "6run"
cd ../6run/
python3 gen_noise.py
sbatch jobscript

echo "7run"
cd ../7run/
python3 gen_noise.py
sbatch jobscript

echo "8run"
cd ../8run/
python3 gen_noise.py
sbatch jobscript

echo "9run"
cd ../9run/
python3 gen_noise.py
sbatch jobscript

deactivate

echo "Finished"