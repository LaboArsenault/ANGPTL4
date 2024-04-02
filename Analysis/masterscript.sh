#!/bin/bash
myArray=(drugpipeline.sh 1a_create_out2.R 1b_resper_tgsd2.R 1d_cophescan.R 2a_ukb_pheno.R 2b_ukb_lof.R 2c_ukb_diseases.R 2d_create_lofpheno.R 3a_hyprcoloc.R 4a_plotallfigures.R 5a_writetables.R 6_writetext.R)
for str in ${myArray[@]}; do
chmod u+x ./$str
done
echo 'Initializing drugpipeline.sh' && ./drugpipeline.sh && echo 'Initializing 1a_create_out2.R' && ./1a_create_out2.R && echo 'Initializing 1b_resper_tgsd2.R' && ./1b_resper_tgsd2.R && echo 'Initializing 1d_cophescan.R' && ./1d_cophescan.R && echo 'Initializing 2a_ukb_pheno.R' && ./2a_ukb_pheno.R && echo 'Initializing 2b_ukb_lof.R' && ./2b_ukb_lof.R && echo 'Initializing 2c_ukb_diseases.R' && ./2c_ukb_diseases.R && echo 'Initializing 2d_create_lofpheno.R' && ./2d_create_lofpheno.R && echo 'Initializing 3a_hyprcoloc.R' && ./3a_hyprcoloc.R && echo 'Initializing 4a_plotallfigures.R' && ./4a_plotallfigures.R && echo 'Initializing 5a_writetables.R' && ./5a_writetables.R && echo 'Initializing 6_writetext.R' && ./6_writetext.R && echo 'The master script finished without errors'
