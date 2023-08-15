# VirolSDC

This is an HBV integration site detection tool.

## Installation

This tool requires a virtual environment on conda.<br>
Download the source code from github.
```
git clone https://github.com/QTazimi/ViroISDC.git
cd ViroISDC
```
Create a python virtual environment and install the required packages. 
If your device is cuda available, you can choose to use tensorflow with gpu.
```
conda create -n virolsdc python=3.6
conda activate virolsdc
pip install -r requirements.txt
```
## Prepare data

## Download

First, use prefetch to download the SRA sample data in batches, 
and then use fastq-dump to convert the sra data into fastq format. Examples are as follows:
```
prefetch --option-file SraAccList.txt
fastq-dump  sample.sra --split-3 --gzip -O /home/data/sample/
```
## Preprocess

Before preprocessing the data, please ensure that you have installed bwa, trimmomatic, and samtools.
Additionally, you will need the human reference genome sequence GRCh37_latest_genomic.fna as well as the reference sequence for Hepatitis B virus (hbv.fasta). 
Below is an example of data preprocessing.
```
Trimmomatic-0.39/trimmomatic-0.39.jar PE -threads 2 -phred33 sample_1.fastq.gz sample_2.fastq.gz LEADING:5 SLIDINGWINDOW:5:20 MINLEN:20 -baseout  sample_trim.fastq.gz
bwa mem -t 4 -R "@RG\tID:hbv\tPL:illumina\tSM:hbv" hbv.fasta sample_trim_1P.fastq.gz   sample_trim_1P.fastq.gz \
            | samtools view -bS - > sample.hbv.bam
samtools view -bF 4 sample.hbv.bam | samtools fastq -1 sample.hbvmap.r1.fq -2 sample.hbvmap.r2.fq -s sample_single.fq -
bwa mem -t 2 -R "@RG\tID:human\tPL:illumina\tSM:human" GRCh37_latest_genomic.fna sample.hbvmap.r1.fq sample.hbvmap.r2.fq >sample.hbvhuman.sam
python data_process/sam_filter.py sample.hbvhuman.sam sample.fil.sam 
samtools view -bS sample.fil.sam | samtools sort -o samplefil_sort.bam -
samtools index samplefil_sort.bam
```
The output file (samplefile_sort.bam) contains preprocessed data that can be directly used for integration site detection. 
Furthermore, the aforementioned data preprocessing steps can be scripted into a shell script for batch processing.

## Usage

## Quick Start

```
python ViroISDC/get_integration.py \
--task_name=hbv \
--do_predict=true \
--vocab_file= VirolSDC/config/word2index.txt \
--bert_config_file= VirolSDC/config/bert_config.json \
--init_checkpoint=/home/HBV/maxpooltest/class_output/model.ckpt-68000 \
--region_path= VirolSDC/file/region.csv \
--fasta_path=/home/hbv/ref/GRCh38/GRCh37_latest_genomic.fna \
--bam_path=samplefil_sort.bam \
--threads=30 \
--filter_depth=5 \
--output_dir=VirolSDC/result/
```
The results of the tool will be generated in the VirolSDC/result/ directory, where the integrations.csv file contains the integration sites of the samples. 
If further filtering is required, you can use the data_process/filter_tool.py for filtering. The usage is as follows:
```
python data_process/filter_tool.py -op VirolSDC/result/integrations.csv -sp VirolSDC/result/integrations_filter.csv
```
The integrations_filter.csv contains the filtered results.

## Example

[gen_filbam.sh](https://github.com/QTazimi/ViroISDC/blob/main/bash/gen_filbam.sh) is an example script for data preprocessing.<br>
[runme.sh](https://github.com/QTazimi/ViroISDC/blob/main/bash/runme.sh) is an example script for running the integration site detection tool.
