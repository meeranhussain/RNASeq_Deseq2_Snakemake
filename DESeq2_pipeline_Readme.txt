# STEPS TO RUN RNA_SEQ SNAKEMAKE FILE
STEP-1  Make  a project folder with Project_ID
STEP-2  Copy Snakefile, Deseq2_final.R, create_combinations.R, config.yaml and Master_file.txt file into the project folder
STEP-3  Create a sub-folder "1_Data"
STEP-4  Copy sample files to 1_Data        Note: Use underscore(_) instead Hyphen(-) Ex: replace'-' with '_' eg. Tumor-1_R1.fq.gz --> Tumor_1_R1.fq.gz
STEP-5  Create "Master_file.txt" to specify the combinations and replicates. Please refer to example file for better clarity
STEP-6  Use config file to add additional information
STEP-7  Open Terminal in Project folder
STEP-8  Type command "snakemake --configfile=config.yaml --cores 5" # cores can be specified based on availability

