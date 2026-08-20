import os
import sys
import csv
import copy

input=[]
for s in sys.argv[1:]:
    input.append(s)

genes=[]
info=[]
heptamers=[]
nonamers=[]
for i in input:
    with open(i,"r") as read:
        reader=csv.reader(read)
        header=next(reader)
        for row in reader:
            if row[:4] not in genes:
                genes.append(row[:4])
            else:
                continue
            inf=copy.deepcopy(row)
            inf.append(i)
            info.append(inf)
            if row[9] not in heptamers:
                heptamers.append(row[9])
            if row[10] not in nonamers:
                nonamers.append(row[10])
        read.close()

for h in heptamers:
    num=0
    neg=0
    pos=0
    dupe=0
    not_dupe=0
    direct=0
    one=0
    two=0
    species=0
    holo=0
    holo_list=[]
    species_list=[]
    for i in info:
        if h==i[9]:
            if i[0] not in holo_list:
                holo_list.append(i[0])
                not_dupe+=1
                holo+=1
            else:
                if dupe==0:
                    dupe+=2
                    not_dupe-=1
                else:
                    dupe+=1
            if i[0].split("/")[1] not in species_list:
                species_list.append(i[0].split("/")[1])
                species+=1
            num+=1
            if i[4]=="-":
                neg+=1
            else:
                pos+=1
            if "directdown" in i[11]:
                direct+=1
            elif "onedown" in i[11]:
                one+=1
            elif "twodown" in i[11]:
                two+=1
    print(h," ",num," ",neg," ",pos," ",not_dupe," ",dupe," ",direct," ",one," ",two," ",holo," ",species)

for h in nonamers:
    num=0
    neg=0
    pos=0
    dupe=0
    not_dupe=0
    direct=0
    one=0
    two=0
    species=0
    holo=0
    holo_list=[]
    species_list=[]
    for i in info:
        if h==i[10]:
            if i[0] not in holo_list:
                holo_list.append(i[0])
                not_dupe+=1
                holo+=1
            else:
                if dupe==0:
                    dupe+=2
                    not_dupe-=1
                else:
                    dupe+=1
            if i[0].split("/")[1] not in species_list:
                species_list.append(i[0].split("/")[1])
                species+=1
            num+=1
            if i[4]=="-":
                neg+=1
            else:
                pos+=1
            if "directdown" in i[11]:
                direct+=1
            elif "onedown" in i[11]:
                one+=1
            elif "twodown" in i[11]:
                two+=1
    print(h," ",num," ",neg," ",pos," ",not_dupe," ",dupe," ",direct," ",one," ",two," ",holo," ",species)

with open("combined_rss_zones.csv","w",newline="") as com:
    writer=csv.writer(com)
    writer.writerow(['Source','GeneType','Contig','Pos','Strand','Sequence','Productive','Locus','Extracted Section (gene +-50bp)','Heptamer','Nonamer',"Number of bp Downstream (RSS)"])
    for io in info:
        if "directdown" in io[-1]:
            io[-1]="0"
        elif "onedown" in io[-1]:
            io[-1]="1"
        elif "twodown" in io[-1]:
            io[-1]="2"
        writer.writerow(io)
    com.close()
    