#! /bin/bash
setup="simu1"
obsnoisevals=(0 0.5 1)


      for obsnoise in "${obsnoisevals[@]}"
      do

	      Rscript run_oracle.R $setup $obsnoise 2>&1 | tee 'logs/'$setup'-'$obsnoise'.log' &
	      echo "Rscript run_oracle.R $setup $obsnoise 2>&1 | tee 'logs/'$setup'-'$obsnoise'.log' &"
     done
