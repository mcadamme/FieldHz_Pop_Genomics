#script to run msms - simulations of field populations
#03052020; updated 3/26/2021

cd /home/megan/ms_sims/FieldPops

mkdir N0_21thou N0_63thou N0_105thou N0_140thou

#simulating 40kb fragments, where present pop based on No = 21,000, u = 2.9e-9 (from Keightly et al. 2015; used by Anderson et al. 2018), and r (4cM/Mb from Martin et al. 2019).  We take temporal samples (-eA), assuming that our field population has undergone 75 generations between sampling points (6 generations/year based on Pan et al. 2016 for 15 years).
/home/megan/src/msdir/ms 24 20000 -t 9.7 -r 134 40000 -eA 0.00089 1 24 > ./N0_21thou/N0_21thou.txt

#allowing present population size at current time to vary - 63,000
/home/megan/src/msdir/ms 24 20000 -t 29.2 -r 403 40000 -eA 0.00030 1 24 > ./N0_63thou/N0_63thou.txt

#allowing present population size at current time to vary - 105,000
/home/megan/src/msdir/ms 24 20000 -t 48.7 -r 672 40000 -eA 0.00018 1 24 > ./N0_105thou/N0_105thou.txt

#allowing present population size at current time to vary - 140,000
/home/megan/src/msdir/ms 24 20000 -t 64.9 -r 896 40000 -eA 0.00013 1 24 > ./N0_140thou/N0_140thou.txt


#Now simulating 10kb fragments, but using exactly the same parameter settings as above.
mkdir N0_21thou_10kb N0_63thou_10kb N0_105thou_10kb N0_140thou_10kb

/home/megan/src/msdir/ms 24 20000 -t 2.4 -r 33.6 10000 -eA 0.00089 1 24 > ./N0_21thou_10kb/N0_21thou.txt

#allowing present population size at current time to vary - 63,000
/home/megan/src/msdir/ms 24 20000 -t 7.3 -r 100.8 10000 -eA 0.00030 1 24 > ./N0_63thou_10kb/N0_63thou.txt

#allowing present population size at current time to vary - 105,000
/home/megan/src/msdir/ms 24 20000 -t 12.2 -r 168 10000 -eA 0.00018 1 24 > ./N0_105thou_10kb/N0_105thou.txt

#allowing present population size at current time to vary - 140,000
/home/megan/src/msdir/ms 24 20000 -t 16.2 -r 224 10000 -eA 0.00013 1 24 > ./N0_140thou_10kb/N0_140thou.txt

