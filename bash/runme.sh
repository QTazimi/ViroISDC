#!/bin/bash

echo '--start--'

path=/home/HBV/SRA335342/
ls $path | while read file
do
  echo $file
  mkdir -p ${path}${file}/VirolSDC1

#  samtools sort ${path}${file}/bio_bert/${file}.hbvhuman.sam -o ${path}${file}/VirolSDC1/${file}_sorted.hbvhuman.bam
#  samtools index ${path}${file}/VirolSDC1/${file}_sorted.hbvhuman.bam

  python /home/program/ViroISDC/get_integration.py \
  --task_name=hbv \
  --do_predict=true \
  --vocab_file=/home/hbv/program/config/word2index.txt \
  --bert_config_file=/home/hbv/program/config/bert_config.json \
  --init_checkpoint=/home/HBV/maxpooltest/class_output/model.ckpt-68000 \
  --region_path= /home/program/VirolSDC/file/region.csv \
  --fasta_path=/home/hbv/ref/GRCh38/GRCh37_latest_genomic.fna \
  --bam_path=${path}${file}/bio_bert/${file}fil_sort.bam \
  --threads=30 \
  --filter_depth=5 \
  --output_dir=${path}${file}/VirolSDC1/

  python /home/program/ViroISDC/data_process/filter_tool.py -op ${path}${file}/VirolSDC1/integrations.csv -sp ${path}${file}/VirolSDC1/integrations_filter.csv


done

echo '--end--'
