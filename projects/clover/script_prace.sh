ensembles="beta1.726/cA211ab.12.48  beta1.726/cA211ab.30.32  beta1.726/cA211ab.40.24  beta1.726/cA211ab.53.24 beta1.778/cB211ab.072.64 beta1.778/cB211ab.14.64 beta1.778/cB211ab.25.48  beta1.836/cC211ab.06.80 "
ensembles="$ensembles beta1.778/cB211ab.25.32  beta1.778/cB211ab.25.24  "
#ensembles="beta1.778/cB211ab.25.32  beta1.778/cB211ab.25.24"


pdf=no
sampling=jack

make

for e in $ensembles
do
#	t=`head -n 1 ../../$e/analysis/main/plateaux_masses.txt | awk '{print $1"   "$2}'`
#        echo $e"    "  $t   	
     ./form_factors read_plateaux -p ../../$e/analysis/main $sampling $pdf
 # ./petros_correlators read_plateaux -p ../../$e/analysis/main jack $pdf
done
#  ./form_factors_out_max_twist read_plateaux -p ../../beta1.726/cA211ab.12.48_no_rew/analysis/main $sampling $pdf
#  ./petros_correlators read_plateaux -p ../../beta1.726/cA211ab.12.48/analysis/main jack $pdf
