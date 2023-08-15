# -*- coing:utf-8 -*-
import csv


def get_sorted(row_datas):
    tempdir = {}
    for line in row_datas:
        if tempdir.get(line[0]) == None:
            tempdir[line[0]] = [line[0:]]
        else:
            tempdir[line[0]].append(line[0:])
    res = []
    for key, value in tempdir.items():
        value.sort(key=lambda s: int(s[1]))
        res.extend(value)
    return res

#####过滤比较接近的断点
def filter_breakpoint(original_path, save_path):
    print('start!')
    row_datas = [] #先收集过滤
    temp_row_datas = []
    with open(original_path, "r", encoding='utf-8', newline='\n') as f:
        reader = csv.reader(f, delimiter=",", quotechar=None)
        for (i, line) in enumerate(reader):
            if i == 0 or int(line[3]) < 5:
                continue
            else:
                temp_row_datas.append(line)

    temp_row_datas = get_sorted(temp_row_datas)
    # print('qqq',temp_row_datas)
    temp_line = temp_row_datas[0]
    for line in temp_row_datas[1:]:
        if float(line[2]) < 0.8:
            print('skip1', line)
            continue
        elif line[0] != temp_line[0] or int(line[1]) > int(temp_line[1]) + 30:
            row_datas.append(temp_line)
            temp_line = line
        else:
            if int(line[3]) > int(temp_line[3]):
                print('skip2', temp_line)
                temp_line = line
                # row_datas.append(line)
            else:
                print('skip3', line)
                # row_datas.append(temp_line)
    if temp_line not in row_datas:
        row_datas.append(temp_line)

    f = open(save_path, 'w', encoding='utf-8', newline='\n')
    writer = csv.writer(f)
    writer.writerow(['Chr', 'Breakpoint', 'Reliability', 'Depth'])
    #['Chr', 'break_point', 'type', 'match', 'softclip', 'hardxlip', 'other', 'total', 'match', 'softclip', 'hardxlip', 'other', 'total']
    f.close()
    with open(save_path, "a", encoding='utf-8', newline='\n') as f:
        writer = csv.writer(f)
        writer.writerows(row_datas)

    print('ALL Done!')



import argparse

def args_func():
    parser = argparse.ArgumentParser(description="Use this function to filter breakpoints about chrr")
    parser.add_argument('--original_path', '-op', help='original path about breakpoint,for example /home/cap4/cap4_filter.csv')
    parser.add_argument('--save_path', '-sp', help='filter breakpoint path,for example /home/cap4/cap4.csv')
    args = parser.parse_args()#['-h']
    return args


def main():
    args = args_func()
    original_path = args.original_path
    save_path = args.save_path
    if original_path == None or save_path == None:
        print('Parameter cannot be empty！')
    else:
        print('Parameter:', original_path, save_path)
        filter_breakpoint(original_path, save_path)


if __name__ == '__main__':
    main()


# original_path = "/mnt/xiaohuang/qiaolei/data/HBV/region/filter/cap4.csv"
# save_path = "/mnt/xiaohuang/qiaolei/data/HBV/region/filter/cap4_filter.csv"
# filter_breakpoint(original_path, save_path)