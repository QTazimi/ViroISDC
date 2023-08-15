# -*- coing:utf-8 -*-
from pysam import VariantFile, FastaFile, AlignmentFile
import pysam
import csv
import os

def get_samread_region(bam_path, save_path):
    bf = pysam.AlignmentFile(bam_path, 'rb')
    # if os.path.isfile(save_path):
    #     pass
    # else:
    #     os.makedir(save_path)
    f = open(save_path, 'w', encoding='utf-8', newline='\n')
    writer = csv.writer(f)
    writer.writerow(['Chr', 'left', 'right', 'length'])
    f.close()
    row_datas = []
    for r in bf:
        if r.reference_name == None:
            continue
        row_datas.append([r.reference_name, r.pos, r.pos + r.query_alignment_end - r.query_alignment_start,
                          r.query_alignment_end - r.query_alignment_start])
        # print(r.reference_name, r.pos, r.pos + r.query_alignment_end - r.query_alignment_start, r.query_alignment_end - r.query_alignment_start)
        # print(r.reference_name, r.pos, r.mapq, r.isize, r.reference_end)
    row_datas.sort(key=lambda range_item: range_item[0])
    with open(save_path, "a", encoding='utf-8', newline='\n') as f:
        # print(row_datas)
        writer = csv.writer(f)
        writer.writerows(row_datas)

def get_union(rangea, rangeb):
    rangea_lower, rangea_upper = rangea
    rangeb_lower, rangeb_upper = rangeb
    if rangea_upper + 1 == rangeb_lower:
        return [1, rangea_lower, rangeb_upper]
    if rangea_upper < rangeb_lower:
        return [2, rangea, rangeb]
    if rangea_lower <= rangeb_lower and rangea_upper >= rangeb_upper:
        return [1, rangea_lower, rangea_upper]
    elif rangea_lower <= rangeb_lower and rangea_upper < rangeb_upper:
        return [1, rangea_lower, rangeb_upper]

def get_nset_union(section_lists):
    section_lists.sort(key=lambda range_item: range_item[0])
    i, j = 0, len(section_lists) - 1
    if j < 0:
      return []
    while True:
        if i == j:
            break
        rangea, rangeb = section_lists.pop(i),  section_lists.pop(i)
        union_item = get_union(rangea, rangeb)
        if union_item[0] == 1:
            section_lists.insert(i, [union_item[1], union_item[2]])
        elif union_item[0] == 2:
            section_lists.insert(i, union_item[1])
            section_lists.insert(i + 1, union_item[2])
            i = i + 1
        j = len(section_lists) - 1
    return section_lists

def get_merge_region(orginal_path, save_path):
    row_datas = {}
    with open(orginal_path, 'r') as f:
        reader = csv.reader(f, delimiter=",", quotechar=None)
        lines = []
        # for line in reader:NC_000001.10
        temp_chrr = 'chr1'
        for (i, line) in enumerate(reader):
            if i == 0:
                continue
            elif i == 1:
                temp_chrr = line[0]
                lines.append([int(line[1]) + 1, int(line[2]) + 1])
            elif line[0] == temp_chrr:
                lines.append([int(line[1]) + 1, int(line[2]) + 1])
            else:
                row_datas[temp_chrr] = get_nset_union(lines)
                temp_chrr = line[0]
                lines = []
                lines.append([int(line[1]) + 1, int(line[2]) + 1])
        row_datas[temp_chrr] = get_nset_union(lines)
    # x = 0
    for k, v in row_datas.items():
        # xx = [item[1] - item[0] for item in row_datas[k]]
        # x = x + sum(xx)
        row_datas[k] = [item for item in row_datas[k] if item.insert(0, k) == None]

    # for k, v in row_datas.items():
    #     print(k, v)
    # print(x)
    with open(save_path, "w", encoding='utf-8', newline='\n') as f:
        writer = csv.writer(f)
        writer.writerow(['Chr', 'start', 'end'])
        for k, v in row_datas.items():
             writer.writerows(row_datas[k])

# get_merge_region("/mnt/xiaohuang/qiaolei/data/HBV/region/testregion_raw.csv", "/mnt/xiaohuang/qiaolei/data/HBV/region/testregion.csv")

import argparse

def args_func():
    parser = argparse.ArgumentParser(description="Use this function to get regions about reads")
    parser.add_argument('--bam', '-b', default='/home/program/bio_bert/file/33_sort.fil.bam', help='bam filename,for example /home/cap4/cap4_sort.bam')
    parser.add_argument('--savefile_path', '-sp', default='/home/program/bio_bert/file/a.csv', help='savefile_path,for example /home/cap4/cap4_raw.csv')
    parser.add_argument('--savefilter_path', '-rf', default='/home/program/bio_bert/file/a_f.csv', help='savefile_path,for example /home/cap4/cap4_filter.csv')
    args = parser.parse_args()#['-h']
    return args


def main():
    args = args_func()
    bam_path = args.bam
    save_path = args.savefile_path
    save_filter_path = args.savefilter_path
    if bam_path == None or save_path == None or save_filter_path == None:
        print('Parameter cannot be empty?')
    else:
        print('Parameter:', bam_path, save_path, save_filter_path)
        get_samread_region(bam_path, save_path)
        get_merge_region(save_path, save_filter_path)

if __name__ == '__main__':
    print('start!')
    main()
    print('end!')