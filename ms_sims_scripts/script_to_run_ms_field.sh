#script to run msms - simulations of field populations
#03052020

cd /home/megan/ms_sims/FieldPops

mkdir N0_21thou N0_63thou N0_105thou N0_140thou

#simulating 40,000 kb fragments, where present pop based on No = 21,000, u = 2.9e-9 (from Keightly et al. 2015; used by Anderson et al. 2018), and r (4cM/Mb from Martin et al. 2019).  We take temporal samples (-eA), assuming that our field population has undergone 75 generations between sampling points (6 generations/year based on Pan et al. 2016 for 15 years).
/home/megan/src/msdir/ms 24 20000 -t 9.7 -r 134 40000 -eA 0.00089 1 24 > ./N0_21thou/N0_21thou.txt

#allowing present population size at current time to vary - 63,000
/home/megan/src/msdir/ms 24 20000 -t 29.2 -r 403 40000 -eA 0.00030 1 24 > ./N0_63thou/N0_63thou.txt

#allowing present population size at current time to vary - 105,000
/home/megan/src/msdir/ms 24 20000 -t 48.7 -r 672 40000 -eA 0.00018 1 24 > ./N0_105thou/N0_105thou.txt

#allowing present population size at current time to vary - 140,000
/home/megan/src/msdir/ms 24 20000 -t 64.9 -r 896 40000 -eA 0.00013 1 24 > ./N0_140thou/N0_140thou.txt

