import os
import sys
import csv
import copy
from collections import defaultdict
from Bio.Blast.Applications import NcbiblastnCommandline

def get_largest_contig(csv_file):
    contig_counts = defaultdict(int)
    
    with open(csv_file, 'r') as f:
        reader = csv.DictReader(f, delimiter='\t')
        for row in reader:
            contig = row['Contig']
            contig_counts[contig] += 1
    
    largest_contig = max(contig_counts.items(), key=lambda x: x[1])
    return largest_contig

with open("../input_data/contig_list.csv") as read:
    reader=csv.reader(read)
    header=next(reader)
    contigs=[]
    for row in reader:
        contigs.append(row)
    read.close()

location=sys.argv[1]
locus=sys.argv[2]
limit=sys.argv[3]
order=""
species=""
hap=""
if limit=="-o":
    order=sys.argv[4]
if limit=="-s":
    species=sys.argv[4]
if limit=="-h":
    hap=sys.argv[4]

with open("target_genes.csv","w",newline="") as write_file:
    writer=csv.writer(write_file)
    writer.writerow(["Source","GeneType","Contig","Pos","Strand","Sequence","Productive","Locus"])
    write_file.close()

for f in os.listdir(location):
    if os.path.isdir(location+"/"+f):
        if order!="" and order!=f:
            continue 
        for f1 in os.listdir(location+"/"+f):
            if os.path.isdir(location+"/"+f+"/"+f1):
                if species!="" and species!=f1:
                    continue
                for f2 in os.listdir(location+"/"+f+"/"+f1):
                    if os.path.isdir(location+"/"+f+"/"+f1+"/"+f2):
                        if hap!="" and hap!=f2:
                            continue
                        for f3 in os.listdir(location+"/"+f+"/"+f1+"/"+f2):
                            if f3=="combined_genes_"+locus+".txt":
                                igh_gene_list=location+"/"+f+"/"+f1+"/"+f2+"/"+f3
                                with open(igh_gene_list,"r") as check:
                                    reader=csv.reader(check,delimiter="\t")
                                    header=next(reader)
                                    num=0
                                    for row in reader:
                                        num+=1
                                    if num==0:
                                        move_on=False
                                    else:
                                        move_on=True
                                    check.close()
                                if move_on==True:
                                    with open(igh_gene_list,"r") as g_igh:
                                        reader=csv.reader(g_igh,delimiter="\t")
                                        header=next(reader)
                                        with open(locus+"_gene_fasta.fasta","w",newline="") as igh_f:
                                            igh_f.close()
                                        genes=[]
                                        for row in reader:
                                            for con in contigs:
                                                if str(row[1])==str(con[0]) and str(f1+"/"+f2) in str(con[1]):
                                                    ap=copy.deepcopy(row)
                                                    ap.insert(0,f+"/"+f1+"/"+f2)
                                                    genes.append(ap)
                                                    igh_gene=">"+f+"/"+f1+"/"+f2+"|"+str(row[1])+"|"+str(row[2])+"\n"+str(row[4])+"\n"
                                                    with open(locus+"_gene_fasta.fasta","a",newline="") as igh_f:
                                                        igh_f.write(igh_gene)
                                                        igh_f.close()
                                        with open("target_genes.csv","a",newline="") as write_file:
                                            writer=csv.writer(write_file)
                                            for g in genes:
                                                writer.writerow(g)
