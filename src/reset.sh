#!/bin/bash

rm ./out*
rm ./slurm*
mv ./fields/*.001 ./
rm ./fields/*
rm ./fields/*
mv ./w.001 ./fields/
# mv ./p.001 ./fields/
mv ./files/seed.000 ./
rm ./files/*
mv ./seed.000 ./files/
rm ./curframe.dat
cp ../src/comm/curframe.dat curframe.dat


# rm ./out*
# rm ./slurm*
# rm ./fields/*
# mv ./files/seed.000 ./seed.000
# rm ./files/*
# mv ./seed.000 ./files/seed.000
# rm ./curframe.dat
# cp ../src/comm/curframe.dat curframe.dat

# scp -r ./valadao/works/2025/lag-lyap/franco/a0.16/ . ftle1.* valadao@turbo.to.infn.it:/scratch1/victor/PerFranco/


#!/bin/bash

# ssh -keygen -t rsa -b 4096 
# cat .ssh/id_rsa.pub > .ssh/authorized_keys
# 
# # login node
# ssh -L 9999:localhost:9999 vvaladao@login01-ext.leonardo.cineca.it
# 
# # compute node
# # ssh -L 9999:localhost:9999 vvaladao@loginXX-ext.leonardo.cineca.it ssh -L 9999:localhost:9999 -N lrdnXXXX
# 
# alias connect='source ~/my_env/bin/activate'