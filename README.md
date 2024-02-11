# RNA_seq-analysis
## Steps to Run RNA_SEQ Snakemake File

### Step 1: Make a Project Folder with Project_ID
Create a project folder and give it a meaningful Project_ID.

### Step 2: Copy Files into Project Folder
Copy the following files into the project folder:
- `Snakefile`
- `Deseq2_final.R`
- `create_combinations.R`
- `config.yaml`
- `Master_file.txt`

### Step 3: Create a Sub-folder "1_Data"
Inside the project folder, create a sub-folder named `1_Data`.

### Step 4: Copy Sample Files to 1_Data
Copy the sample files into the `1_Data` folder. Make sure to use underscores (_) instead of hyphens (-) in file names. For example, replace '-' with '_' (e.g., `Tumor-1_R1.fq.gz` --> `Tumor_1_R1.fq.gz`).

### Step 5: Create "Master_file.txt"
Create a file named `Master_file.txt` in the project folder. This file should specify the combinations and replicates. Refer to the example file provided for better clarity.

### Step 6: Use Config File to Add Additional Information
Utilize the `config.yaml` file to add any additional information required for the workflow.

### Step 7: Open Terminal in Project Folder
Navigate to the project folder in your terminal/command prompt.

## Step 8: Run Snakemake
Type the following command in the terminal:
```bash
snakemake --configfile=config.yaml --cores 5

