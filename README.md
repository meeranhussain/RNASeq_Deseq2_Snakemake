# RNA_seq-analysis
This workflow is for differential gene expression study for the samples with replicates
## FOR STAND ALONE RUN (INDIVIDUAL COMMANDS)
**1. Perform quality check on fastq file using FASTQC/MULTIQC**

FASTQC:
```bash
fastqc *.fastq -o <output_directory>
```
-o : Directory to save output files (file must be created)
"*.fastq" represents to select all files with the ".fastq" extension in the working directory

MULTIQC:
```bash
multiqc <fastqc_results_directory>/
```

**2. Quality control using Trimgalore**
```bash
trim_galore -q 20 --paired --fastqc --cores <number_of_threads> <input_R1_fq.gz> <input_R2.fq.gz> -o <output_directory>
```

**3. Indexing reference file**
```bash
STAR --runMode genomeGenerate --genomeDir <index_dir_name> --genomeFastaFiles <path to ".fasta" file> --sjdbGTFfile <path to ".gtf" file> --sjdbOverhang 100 --runThreadN 10
```

**4. Alignment using STAR**
```bash
STAR --genomeDir <index_dir_name> --runThreadN <number_of_threads> --outSAMtype BAM SortedByCoordinate --readFilesCommand zcat --readFilesIn <input_R1.fastq.gz> <input_R2.fastq.gz>  --outFileNamePrefix <output_filename>
```

**5. Read Quantification using Feature count**
```bash
featureCounts -p -T <number_of_threads> --verbose -t exon -g gene_id -a <path to ".gtf" file> -o <out_count_file_name> <List of BAM files as Input files>
```

# FOR SNAKEMAKE RUN
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
Copy the sample files into the `1_Data` folder. If in case you want to use characters in sample name make sure to use underscores (_) instead of hyphens (-) in file names. For example, replace '-' with '_' (e.g., `Tumor-1_R1.fq.gz` --> `Tumor_1_R1.fq.gz`).

### Step 5: Create "Master_file.txt"
Create a file named `Master_file.txt` in the project folder. This file should specify the combinations and replicates. Refer to the example file provided for better clarity.

### Step 6: Use Config File to Add Additional Information
Utilize the `config.yaml` file to add any additional information required for the workflow.

#### Config.yaml Content for RNA_SEQ Snakemake Workflow (Example file)

```yaml
#### Enter organism name (Scientific name)
org: "Homo sapiens"

#### Enter Kegg organism code
org_code: "hsa"

#### Specify Number of threads
threads: "40"

#### Specify Combinations using "+" between combinations
combinations: "control_Tumor + Tumor_control"

#### Path to indexed reference folder (Reference indexing command provided below)
reference: "</Path/to/indexed/reference/folder>"
```
##### Genome indexing using STAR
```bash
STAR --runMode genomeGenerate --genomeDir {index_dir_name} --genomeFastaFiles {path to ".fasta" file} --sjdbGTFfile {path to ".gtf" file} --sjdbOverhang 100 --runThreadN 10
```
### Step 7: Open Terminal in Project Folder
Navigate to the project folder in your terminal/command prompt.

## Step 8: Run Snakemake
Type the following command in the terminal:
```bash
snakemake --configfile=config.yaml --cores 5
```
