# -*- coing:utf-8 -*-

import gzip
from pysam import VariantFile, FastaFile, AlignmentFile
import pysam
import csv


def ref_atcg(fasta_file, chr_id, start, end):
    # print('chr_id', chr_id, start, end)
    for rec in fasta_file.fetch(chr_id, start, end):
        if rec == 'N':
            # print("ref_base is none")
            return 'N'
        return rec


def cigar_info(poss, mapstr_list):
    match_num1, softclip_num1, hardclip_num1, other_num1 = 0, 0, 0, 0
    match_num2, softclip_num2, hardclip_num2, other_num2 = 0, 0, 0, 0
    match_num3, softclip_num3, hardclip_num3, other_num3 = 0, 0, 0, 0
    pos1_cigar_str, pos2_cigar_str, pos3_cigar_str = '', '', ''
    misrate = 0.8
    for pos, mapstr in zip(poss, mapstr_list):
        pos1 = int(pos) - 1
        if pos1 >= 0 and pos1 <= len(mapstr) - 1:
            maptype1 = mapstr[pos1]
            if maptype1 == '0':
                match_num1 += 1
                pos1_cigar_str = pos1_cigar_str + 'M'
            elif maptype1 == '4':
                softclip_num1 += 1
                pos1_cigar_str = pos1_cigar_str + 'S'
            elif maptype1 == '5':
                hardclip_num1 += 1
                pos1_cigar_str = pos1_cigar_str + 'H'
            else:
                other_num1 += 1
                pos1_cigar_str = pos1_cigar_str + 'O'
        pos2 = int(pos)
        if pos2 >= 0 and pos2 <= len(mapstr) - 1:
            maptype2 = mapstr[pos2]
            if maptype2 == '0':
                match_num2 += 1
                pos2_cigar_str = pos2_cigar_str + 'M'
            elif maptype2 == '4':
                softclip_num2 += 1
                pos2_cigar_str = pos2_cigar_str + 'S'
            elif maptype2 == '5':
                hardclip_num2 += 1
                pos2_cigar_str = pos2_cigar_str + 'H'
            else:
                other_num2 += 1
                pos2_cigar_str = pos2_cigar_str + 'O'
        pos3 = int(pos) + 1
        if pos3 >= 0 and pos3 <= len(mapstr) - 1:
            maptype3 = mapstr[pos3]
            if maptype3 == '0':
                match_num3 += 1
                pos3_cigar_str = pos3_cigar_str + 'M'
            elif maptype3 == '4':
                softclip_num3 += 1
                pos3_cigar_str = pos3_cigar_str + 'S'
            elif maptype3 == '5':
                hardclip_num3 += 1
                pos3_cigar_str = pos3_cigar_str + 'H'
            else:
                other_num3 += 1
                pos3_cigar_str = pos3_cigar_str + 'O'

    total_depth1 = match_num1 + softclip_num1 + hardclip_num1 + other_num1
    total_depth2 = match_num2 + softclip_num2 + hardclip_num2 + other_num2
    total_depth3 = match_num3 + softclip_num3 + hardclip_num3 + other_num3
    cigar_str = [pos1_cigar_str, pos2_cigar_str, pos3_cigar_str]
    if total_depth1 == 0 or total_depth2 == 0 or total_depth3 == 0:
        break_point_infor = ['2']
        return break_point_infor
    else:
        if (float(hardclip_num1) + float(softclip_num1)) / float(total_depth1) > misrate and float(match_num2) / float(
                total_depth2) > misrate and float(match_num3) / float(total_depth3) > misrate:
            # print('break_pointL:')
            break_point_infor = ['1', cigar_str]
            return break_point_infor
        elif float(match_num1) / float(total_depth1) > misrate and float(match_num2) / float(
                total_depth2) > misrate and (float(hardclip_num3) + float(softclip_num3)) / float(
            total_depth3) > misrate:
            # print('break_pointR:')
            break_point_infor = ['1', cigar_str]
            return break_point_infor
        # elif float(match_num1) / float(total_depth1) > misrate and float(match_num2) / float(total_depth2) > misrate \
        #         and float(match_num3) / float(total_depth3) > misrate:
        #     # print('not break_point')
        #     break_point_infor = ['3', cigar_str]
        #     return break_point_infor
        else:
            break_point_infor = ['2']
            return break_point_infor


def get_pos_cigar(chr, pos, bam_file):
    for rec in bam_file.pileup(chr, pos, pos + 1, stepper='nofilter', ignore_overlaps=True):
        if rec.pos == pos:
            poss = rec.get_query_positions()
            # print(poss)
            align_list = [tmp.alignment for tmp in rec.pileups]
            # cigar
            # print('cigar')
            cigartuples = [read.cigartuples for read in align_list]
            # print(cigartuples)

            mapstr_list = []
            for item in cigartuples:
                stritem = ''
                for tupleitem in item:
                    stritem = stritem + str(tupleitem[0]) * tupleitem[1]
                mapstr_list.append(stritem)
            # print(mapstr_list)
            pileup_infor = [poss, mapstr_list]
            return pileup_infor


def get_breakpoint(bam_path, chrr, region_left, region_right, fasta_path):
    print('start')
    bam_file = AlignmentFile(bam_path, 'rb')
    fasta_file = FastaFile(fasta_path)
    # break_point = 117337870 - 1
    rows_data = []
    for break_point in range(region_left, region_right):  # chrr_len
        break_point = break_point - 1
        ref1 = ref_atcg(fasta_file, chrr, break_point - 1, break_point)
        ref2 = ref_atcg(fasta_file, chrr, break_point, break_point + 1)
        ref3 = ref_atcg(fasta_file, chrr, break_point + 1, break_point + 2)
        if ref1 == 'N' or ref2 == 'N' or ref3 == 'N':
            continue
        pileup_infor = get_pos_cigar(chrr, break_point, bam_file)
        if pileup_infor == None:
            continue
        poss, mapstr_list = pileup_infor[0], pileup_infor[1]
        if len(poss) == 0 or len(mapstr_list) == 0:
            continue
        break_point_infor = cigar_info(poss, mapstr_list)
        flag = break_point_infor[0]
        #返回可疑窗口的数据
        if flag == '1':
            rows_data.append([chrr, break_point + 1, break_point_infor[1][0], break_point_infor[1][1], break_point_infor[1][2]])

    return rows_data


#############################
from multiprocessing import Pool
import os, time
import fcntl


def long_time_task(bam_path, chrr, region_left, region_right, save_path, fasta_path):
    print('Run task %s (%s)...' % (chrr, os.getpid()))
    start = time.time()

    print('chrr', chrr, 'chrr_left:', region_left, 'chrr_right:', region_right)
    row_datas = get_breakpoint(bam_path, chrr, region_left, region_right, fasta_path)
    if len(row_datas) != 0:
        with open(save_path, "a", encoding='utf-8', newline='\n') as f:
            fcntl.flock(f.fileno(), fcntl.LOCK_EX)
            writer = csv.writer(f)
            writer.writerows(row_datas)

    # time.sleep(1)
    end = time.time()
    print('Task %s runs %0.2f seconds.' % (chrr, (end - start)))
    return 1


def total_process(fasta_path, bam_path, region_path, save_path, threads):
    bam_file = AlignmentFile(bam_path, 'rb')
    bam_list = bam_file.get_index_statistics()
    ll = []
    for item in bam_list:
        if item.total != 0:
            ll.append(item.contig)
    bam_file.close()
    print(ll)

    with open(region_path, 'r') as f:
        reader = csv.reader(f, delimiter=",", quotechar=None)
        lines = []
        # for line in reader:
        for (i, line) in enumerate(reader):
            if line[0] == 'Chr':
                continue
            else:
                lines.append(line)

    f = open(save_path, 'w', encoding='utf-8', newline='\n')
    writer = csv.writer(f)
    writer.writerow(
        ['Chr', 'break_point', 'pos1_cigar_str', 'pos2_cigar_str', 'pos3_cigar_str'])
    f.close()

    print('Parent process %s.' % os.getpid())
    p = Pool(threads)
    results = []
    for line in lines:
        chrr = line[0]
        region_left = int(line[1]) - 20
        region_right = int(line[2]) + 20
        if region_left <= 0:
            region_left = 2
        # print(chrr, region_left, region_right)
        if chrr in ll:
            # time.sleep(1)
            result = p.apply_async(long_time_task, args=(bam_path, chrr, region_left, region_right, save_path, fasta_path))
            results.append(result)
    data = [r.get() for r in results]
    # print(len(data))
    print('Waiting for all subprocesses done...')
    p.close()
    p.join()
    print('All subprocesses done.')


#########################################################


import argparse

def args_func():
    parser = argparse.ArgumentParser(description="Use this function to get breakpoints about chrr")
    parser.add_argument('--fasta', '-fa', default="/home/hbv/ref/GRCh38/GRCh37_latest_genomic.fna", help='fasta filename')
    parser.add_argument('--bam', '-b', default="/home/program/bio_bert/file/33_sort.fil.bam", help='bam filename,for example /home/cap4/cap4_sort.bam')
    parser.add_argument('--region_path', '-rf', default="/home/program/bio_bert/file/a_f.csv", help='region_path,for example /home/cap4/cap4_filter.csv')
    parser.add_argument('--savefile_path', '-sp', default="/home/program/bio_bert/file/test.csv", help='savefile_path,for example /home/cap4/cap4.csv')
    parser.add_argument('--threads', '-t', default=20, help='Number of thread')
    args = parser.parse_args()#['-h']
    return args


def main():
    args = args_func()
    fasta_path = args.fasta
    bam_path = args.bam
    save_path = args.savefile_path
    region_path = args.region_path
    threads = args.threads
    if fasta_path == None or bam_path == None or save_path == None or region_path == None:
        print('Parameter cannot be empty！')
    else:
        print('Parameter:', fasta_path, bam_path, region_path, save_path)
        total_process(fasta_path, bam_path, region_path, save_path, threads)



if __name__ == '__main__':
    start = time.time()
    main()
    end = time.time()
    print('Task %s runs %0.2f seconds.' % ('all', (end - start)))