#!/usr/bin/zsh

echo '--start--'

start_time=`date +%s`  #����ű����еĿ�ʼʱ��

tmp_fifofile="/tmp/$$.fifo"
mkfifo $tmp_fifofile   # �½�һ��FIFO���͵��ļ�
exec 6<>$tmp_fifofile  # ��FD6ָ��FIFO����
rm $tmp_fifofile  #ɾҲ���ԣ�

thread_num=15  # ��������߳���

path=/home/HBV/SRA335342/

#�����߳��������������Ƹ���
#��ʵ�Ͼ�����fd6�з�����$thread_num���س���
for ((i=0;i<${thread_num};i++));do
    echo
done >&6
ls $path | while read file
do
    read -u6
    {
    
        echo $file
        mkdir -p /home/HBV/SRA335342/${file}/bio_bert
        #java -jar /home/package/Trimmomatic-0.39/trimmomatic-0.39.jar PE -threads 2 -phred33 /home/HBV/SRA335342/${file}/${file}_1.fastq.gz /home/HBV/SRA335342/${file}/${file}_2.fastq.gz LEADING:5 SLIDINGWINDOW:5:20 MINLEN:20 -baseout  /home/HBV/SRA335342/${file}/${file}_trim.fastq.gz

        bwa mem -t 4 -R "@RG\tID:hbv\tPL:illumina\tSM:hbv" /home/HBV_wgsim/ref/10791.fasta /home/HBV/SRA335342/${file}/${file}_trim_1P.fastq.gz   /home/HBV/SRA335342/${file}/${file}_trim_2P.fastq.gz \
            | samtools view -bS - > /home/HBV/SRA335342/${file}/bio_bert/${file}.hbv.bam
        
        
        samtools view -bF 4 /home/HBV/SRA335342/${file}/bio_bert/${file}.hbv.bam | samtools fastq -1 /home/HBV/SRA335342/${file}/bio_bert/${file}.hbvmap.r1.fq -2 /home/HBV/SRA335342/${file}/bio_bert/${file}.hbvmap.r2.fq -s /home/HBV/SRA335342/${file}/bio_bert/${file}_single.fq -
        
        
        bwa mem -t 2 -R "@RG\tID:human\tPL:illumina\tSM:human" /home/hbv/ref/GRCh38/GRCh37_latest_genomic.fna /home/HBV/SRA335342/${file}/bio_bert/${file}.hbvmap.r1.fq /home/HBV/SRA335342/${file}/bio_bert/${file}.hbvmap.r2.fq > /home/HBV/SRA335342/${file}/bio_bert/${file}.hbvhuman.sam
        
        
        python /home/HBV_wgsim/lincode/tools/sam_filter.py /home/HBV/SRA335342/${file}/bio_bert/${file}.hbvhuman.sam /home/HBV/SRA335342/${file}/bio_bert/${file}.fil.sam 
        
        samtools view -bS /home/HBV/SRA335342/${file}/bio_bert/${file}.fil.sam | samtools sort -o /home/HBV/SRA335342/${file}/bio_bert/${file}fil_sort.bam -
        samtools index /home/HBV/SRA335342/${file}/bio_bert/${file}fil_sort.bam
        
                       
        echo >&6
    } &
      
done


wait # Ҫ��wait���ȴ������߳̽���

stop_time=`date +%s` # ����ű����еĽ���ʱ��

echo "TIME:`expr $stop_time - $start_time`" # ����ű�����ʱ��

exec 6>&- # �ر�FD6
echo '--done--'  # ��ʾ�ű����н���

