#! /bin/bash
setup="simu2"
obsnoisevals=(0 0.5)
betas=(0.5 1)
sigmas=(1 3)


  for obsnoise in "${obsnoisevals[@]}"
  do
    for beta in "${betas[@]}"
    do
      for sigma in "${sigmas[@]}"
      do


	     Rscript run_oracle.R $setup $obsnoise $beta $sigma 2>&1 | tee 'logs/'$setup'-'$obsnoise'-'$beta'-'$sigma'.log' &
	     echo "Rscript run_oracle.R $setup $obsnoise $beta $sigma 2>&1 | tee 'logs/'$setup'-'$obsnoise'-'$beta'-'$sigma'.log' &"
      	done
     done
  done
