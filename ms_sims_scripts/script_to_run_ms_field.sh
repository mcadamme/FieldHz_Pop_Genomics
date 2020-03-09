#script to run msms - simulations of field populations
#03052020

cd /home/megan/ms_sims/FieldPops

mkdir Ne_fivethou_NoSize Ne_twopointfivethou_NoSize Ne_tenthou_NoSize Ne_twentythou_NoSize

#simulating present pop based on No = 5000, u = 2.9e-9 (from Keightly et al. 2015; used by Anderson et al. 2018), and r (4cM/Mb from Martin et al. 2019).  We take temporal samples (-eA), assuming that our field population has undergone 75 generations between sampling points (5 generations/year based on Pan et al. 2016 for 15 years).
/home/megan/src/msdir/ms 24 20000 -t 2.32 -r 32 40000 -eA 0.00375 1 24 > ./Ne_fivethou_NoSize/Ne_fivethou_NoSize.txt

#allowing present population size at current time to vary - 2500
/home/megan/src/msdir/ms 24 20000 -t 1.16 -r 16 40000 -eA 0.0075 1 24 > ./Ne_twopointfivethou_NoSize/Ne_twopointfivethou_NoSize.txt

#allowing present population size at current time to vary - 10000
/home/megan/src/msdir/ms 24 20000 -t 4.64 -r 64 40000 -eA 0.0019 1 24 > ./Ne_tenthou_NoSize.txt

#allowing present population size at current time to vary - 20000
/home/megan/src/msdir/ms 24 20000 -t 9.28 -r 128 40000 -eA 0.0009 1 24 > ./Ne_twentythou_NoSize.txt

