import os
from os.path import join
from pyfiglet import Figlet
f = Figlet(font='slant')
print (f.renderText('RNA-Seq Reference based STAR Pipeline'))
dev = Figlet(font='digital')
print (dev.renderText('Developed : Meeran Hussain'))
#################### Get Project Details ######################################
#user_name = input("Enter your name:")
#print ("Hello",user_name, "Please enter following details:")
org = config["org"]
#print ("Organism:",org)
org_code = config["org_code"]
#print ("Organism code:",org_code)
###############################################################################

######################## Reference Version #################################
reference = config["reference"]
####### Fetch path of specified organism Index ################################
STAR_INDEX = reference
SAMPLE_DIR = "1_Data"
SAMPLES, = glob_wildcards(SAMPLE_DIR + "/{sample}_R1.fq.gz")
###############################################################################

################ Fetch GTF for the mentioned organism ########################
# folder path
dir = ("/mnt/data1/EGI/01_RNASeq/Automation_testing/Reference_file", reference)
dir_path = "/".join(dir)

# list to store files
res = []
# Iterate directory
for file in os.listdir(dir_path):
    # check only text files
    if file.endswith('.gtf'):
        res.append(file)
print(res[0])

gtf_name = res[0]

gtf_path = (dir_path, gtf_name)

gtf = "/".join(gtf_path)
###############################################################################
R1 = '{sample}_R1.fq.gz'
R2 = '{sample}_R2.fq.gz'


rule all:
    input:
        expand("3_Analysis/1_STAR/{sample}/Aligned.sortedByCoord.out.bam", sample=SAMPLES),
        '3_Analysis/2_Deseq2/counts.txt',
        '2_Combinations',
        '3_Analysis/2_Deseq2/DGE_files'

rule alignment:
    input:
        r1 = join(SAMPLE_DIR, R1),
        r2 = join(SAMPLE_DIR, R2)
    params:
        STAR_INDEX = STAR_INDEX,
        #threads = config["threads"]
    threads:
        15
    output:
        "3_Analysis/1_STAR/{sample}/Aligned.sortedByCoord.out.bam"
    message:
        "--- Mapping STAR---"
    shell:"""
        mkdir -p 3_Analysis
        mkdir -p 3_Analysis/1_STAR
        mkdir -p 3_Analysis/2_Deseq2
        mkdir -p 3_Analysis/1_STAR/{wildcards.sample}
        /mnt/data1/ShareApps/STAR-2.7.10a/source/STAR --genomeDir {params.STAR_INDEX}  --runThreadN {threads} --outSAMtype BAM  SortedByCoordinate --readFilesCommand zcat  --readFilesIn {input.r1} {input.r2}  --outFileNamePrefix 3_Analysis/1_STAR/{wildcards.sample}/
    """

rule feature_count:
    input:
        expand("3_Analysis/1_STAR/{sample}/Aligned.sortedByCoord.out.bam", sample=SAMPLES),
        gtf = gtf
    output:
        '3_Analysis/2_Deseq2/counts.txt'
    params:
        inpstring = 'Aligned.sortedByCoord.out.bam'
    threads:
        15
    message:
        "---- Running Feature Count----"
    shell:
        'inputfiles=$(find . -name {params.inpstring} | xargs echo) &&'
        'readlink -m $inputfiles &&'
        '/mnt/data1/ShareApps/subread-2.0.3-source/bin/featureCounts -p -T {threads} --verbose -t exon -g gene_id -a {input.gtf} -o {output} $inputfiles'

rule create_metatable:
    input:
        master_file = "Master_file.txt"
    output:
        directory("2_Combinations")
    params:
        comb = config["combinations"]
    message:
        "--- Creating Combination files ---"
    shell:"""
        mkdir -p {output}
        Rscript /mnt/data1/ShareApps/Automation_scripts/Reference_Replicates_Deseq2/create_combinations.R {input.master_file} {params.comb}
    """

rule Deseq2_Analysis:
    input:
        counts = '3_Analysis/2_Deseq2/counts.txt'
    output:
        directory("3_Analysis/2_Deseq2/DGE_files")
    params:
        org = config["org"],
        org_code = config["org_code"]
    message:
        "---- Running Deseq2 Analysis----"
    shell:"""
        mkdir -p {output}
        mkdir -p 3_Analysis/2_Deseq2/Plots
        Rscript /mnt/data1/ShareApps/Automation_scripts/Reference_Replicates_Deseq2/Deseq2_final.R {input.counts} '{params.org}' {params.org_code}
    """
