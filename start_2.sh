#! /bin/bash
nvals=(250 500 1000 5000 10000 20000 30000)
setup="simu2"
obsnoisevals=(0 0.5)
betas=(0.5 1)
sigmas=(1 3)

repidxs=($(seq 1 1 50))

for n in "${nvals[@]}"
do
  for obsnoise in "${obsnoisevals[@]}"
  do
    for beta in "${betas[@]}"
    do
      for sigma in "${sigmas[@]}"
      do
  	    for repidx in "${repidxs[@]}"
  	    do

	     OMP_NUM_THREADS=1 Rscript main.R $setup $n $obsnoise $beta $sigma $repidx 2>&1 | tee 'logs/'$setup'-'$n'-'$obsnoise'-'$beta'-'$sigma'-'$repidx'.log' &
	     echo "OMP_NUM_THREADS=1 Rscript main.R $setup $n $obsnoise $beta $sigma $repidx 2>&1 | tee 'logs/'$setup'-'$n'-'$obsnoise'-'$beta'-'$sigma'-'$repidx'.log' &"
            done
         done
      done
    done
done 
