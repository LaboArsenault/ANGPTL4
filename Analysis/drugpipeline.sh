#!/bin/bash
start_time=$(date +%s)
myArray=(/mnt/sda/gagelo01/Projects/Pipelines/Robust_efficient_drug_pipeline/Analysis/1aa_create_dt_gene_region.R /mnt/sda/gagelo01/Projects/Pipelines/Robust_efficient_drug_pipeline/Analysis/1b_create_inst.R /mnt/sda/gagelo01/Projects/Pipelines/Robust_efficient_drug_pipeline/Analysis/1c_create_out.R /mnt/sda/gagelo01/Projects/Pipelines/Robust_efficient_drug_pipeline/Analysis/1dd_vep.R /mnt/sda/gagelo01/Projects/Pipelines/Robust_efficient_drug_pipeline/Analysis/1g_multicisandpan.R)
for str in ${myArray[@]}; do
chmod u+x $str
done
echo 'Initializing /mnt/sda/gagelo01/Projects/Pipelines/Robust_efficient_drug_pipeline/Analysis/1aa_create_dt_gene_region.R' && /mnt/sda/gagelo01/Projects/Pipelines/Robust_efficient_drug_pipeline/Analysis/1aa_create_dt_gene_region.R tosource.R && echo 'Initializing /mnt/sda/gagelo01/Projects/Pipelines/Robust_efficient_drug_pipeline/Analysis/1b_create_inst.R' && /mnt/sda/gagelo01/Projects/Pipelines/Robust_efficient_drug_pipeline/Analysis/1b_create_inst.R tosource.R && echo 'Initializing /mnt/sda/gagelo01/Projects/Pipelines/Robust_efficient_drug_pipeline/Analysis/1c_create_out.R' && /mnt/sda/gagelo01/Projects/Pipelines/Robust_efficient_drug_pipeline/Analysis/1c_create_out.R tosource.R && echo 'Initializing /mnt/sda/gagelo01/Projects/Pipelines/Robust_efficient_drug_pipeline/Analysis/1dd_vep.R' && /mnt/sda/gagelo01/Projects/Pipelines/Robust_efficient_drug_pipeline/Analysis/1dd_vep.R tosource.R && echo 'Initializing /mnt/sda/gagelo01/Projects/Pipelines/Robust_efficient_drug_pipeline/Analysis/1g_multicisandpan.R' && /mnt/sda/gagelo01/Projects/Pipelines/Robust_efficient_drug_pipeline/Analysis/1g_multicisandpan.R tosource.R && echo 'The master script finished without errors'
end_time=$(date +%s)
elapsed_time=$((end_time - start_time))
echo "The script took $elapsed_time seconds to run."
