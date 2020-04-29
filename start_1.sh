#! /bin/bash
nvals=(250 500 1000 5000 10000 20000)
setup="simu1"
obsnoisevals=(0 0.5 1)

repidxs=($(seq 1 1 50))

for n in "${nvals[@]}"
do
      for obsnoise in "${obsnoisevals[@]}"
      do
	      for repidx in "${repidxs[@]}"
	      do


	      OMP_NUM_THREADS=1 Rscript main.R $setup $n $obsnoise $repidx 2>&1 | tee 'logs/'$setup'-'$n'-'$obsnoise'-'$repidx'.log' &
	      echo "OMP_NUM_THREADS=1 Rscript main.R $setup $n $obsnoise $repidx 2>&1 | tee 'logs/'$setup'-'$n'-'$obsnoise'-'$repidx'.log' &"
      	done
     done
done
