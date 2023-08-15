import sys

f=open(sys.argv[1],'r')
g=open(sys.argv[2],'w+')

errls = []
sqnlist = []

for each in f:
    try:
        if '@' not in each.split('\t')[0]:
            cigar = each.split('\t')[5]
            if cigar.count('S') > 1 or cigar.count('H') > 1:
                continue
            elif cigar.count('S') == 1:
                if cigar[-1] == 'M':
                    Mnum = cigar.split('S')[1].strip('M')
                    Snum = cigar.split('S')[0]
                else:
                    Mnum = cigar.split('M')[0]
                    Snum = cigar.split('M')[1].strip('S')
                if int(Mnum) < 30 or int(Snum) < 11:
                    continue
                else:
                    seqname = each.split('\t')[0]
                    if seqname not in sqnlist:
                        sqnlist.append(seqname)
                        g.write(each)
                    else:
                        continue
            elif cigar.count('H') == 1:
                if cigar[-1] == 'M':
                    Mnum = cigar.split('H')[1].strip('M')
                    Snum = cigar.split('H')[0]
                else:
                    Mnum = cigar.split('M')[0]
                    Snum = cigar.split('M')[1].strip('H')
                if int(Mnum) < 30 or int(Snum) < 11:
                    continue
                else:
                    seqname = each.split('\t')[0]
                    if seqname not in sqnlist:
                        sqnlist.append(seqname)
                        g.write(each)
                    else:
                        continue    
            elif '*' in cigar:
                continue
            elif cigar.count('S') == 0 and cigar.count('H') == 0:
                continue
            else:
                seqname = each.split('\t')[0]
                if seqname not in sqnlist:
                    sqnlist.append(seqname)
                    g.write(each)
                else:
                    continue      
        else:
            g.write(each)
    except ValueError:
        errls.append(each.split('\t')[5])

print(errls, len(errls))

f.close()
g.close()