import os
import csv
import sys
import shutil
import numpy as np
import math
import copy
import statistics
import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator

data_filter = sys.argv[1] # -none / filter string (ex: Landfowl)  :  locus (ex: IGH)
locus=data_filter.split(":")[1]
data_dir=sys.argv[2] # directory to pull data from
gene_data=sys.argv[3] # -csv (provide gene_list.csv in input_data)  |  -fold : path to IgDetective data folder with (order --> species --> haplotype) subdir structure
output_folder=sys.argv[4] #  -final_out : name of final output folder

data=[]
for folder in os.listdir(data_dir):
    if os.path.isdir(data_dir+"/"+folder):
        with open(data_dir+"/"+folder+"/d_genes.csv","r") as read:
            reader=csv.reader(read)
            header=next(reader)
            found=False
            for d in data:
                if d[0]==locus:
                    found=True
            if found==False:
                data.append([locus,[]])
            for d in data:
                if d[0]==locus:
                    d[1].append([int(folder.split("_")[-1].split("-")[0]),int(folder.split("_")[-1].split("-")[1]),[]])
                    for row in reader:
                        d[1][-1][2].append(row)
            read.close()

for d in data:
    print(d[0])
    d[1].sort(key=lambda x: x[0])
    old_num=0
    new_num=0
    hep=''
    for d2 in d[1]:
        if d2[0]!=hep:
            hep=d2[0]
            if new_num!=0:
                if old_num==0:
                    d[1][:new_num] = sorted(d[1][:new_num], key=lambda x: x[1])
                else:
                    d[1][old_num:new_num] = sorted(d[1][old_num:new_num], key=lambda x: x[1])
                old_num=new_num
        new_num+=1
    d[1][old_num:new_num] = sorted(d[1][old_num:new_num], key=lambda x: x[1])

    
    #d[1][:5] = sorted(d[1][:5], key=lambda x: x[1])
    #d[1][5:10] = sorted(d[1][5:10], key=lambda x: x[1])
    #d[1][10:] = sorted(d[1][10:], key=lambda x: x[1])
    for d1 in d[1]:
        print(d1[0], d1[1], len(d1[2]))

contigs=[]
with open("../input_data/contig_list.csv","r") as read:
    reader=csv.reader(read)
    header=next(reader)
    for row in reader:
        contigs.append(row)
    read.close()

all_haplotypes=[]
if gene_data.startswith("-fold:"):
    gene_fold=gene_data.split(":")[1]
    for c in os.listdir(gene_fold):
        if os.path.isdir(gene_fold+"/"+c):
            for c1 in os.listdir(gene_fold+"/"+c):
                if os.path.isdir(gene_fold+"/"+c+"/"+c1):
                    for c2 in os.listdir(gene_fold+"/"+c+"/"+c1):
                        if os.path.isdir(gene_fold+"/"+c+"/"+c1+"/"+c2):
                            if data_filter.startswith("-none") or data_filter.split(":")[0].replace("-","") in (c+"/"+c1+"/"+c2):
                                if (c+"/"+c1+"/"+c2) not in all_haplotypes:
                                    all_haplotypes.append(c+"/"+c1+"/"+c2)
elif gene_data.startswith("-csv"):
    with open("../input_data/gene_list.csv","r") as read:
        reader=csv.reader(read)
        header=next(reader)
        for row in reader:
            if row[0] not in all_haplotypes:
                all_haplotypes.append(row[0])
        read.close()
print(all_haplotypes)

def data_create(name_data):
    output_data=[]
    for haplo in name_data:
        hap=[haplo,[]]
        for d in data:
            locus=d[0]
            found=False
            for i in hap[1]:
                if locus==i[0]:
                    found=True
                    break
            if found==False:
                hap[1].append([locus,[]])
        for h in hap[1]:
            for con in contigs:
                if con[2]==h[0] and haplo in con[1]:
                    h[1].append([con[0],[]])
        for d in data:
            for d1 in d[1]:
                hep_thresh=d1[0]
                non_thresh=d1[1]
                all_genes=d1[2]

                for gene in all_genes:
                    for i in hap[1]:
                        for c in i[1]:
                            if d[0]==i[0]:
                                found=False
                                for i1 in c[1]:
                                    if hep_thresh==i1[0] and non_thresh==i1[1]:
                                        found=True
                                if found==False:
                                    c[1].append([hep_thresh,non_thresh,[]])
                                
                                if hap[0] in gene[0] and c[0]==gene[2]:
                                    for i1 in c[1]:
                                        if hep_thresh==i1[0] and non_thresh==i1[1]:
                                            i1[2].append(gene)
        output_data.append(hap)
    return output_data

all_haplotype_data = data_create(all_haplotypes)
#all_species_data = data_create(all_species)
#all_order_data = data_create(all_orders)

first_iter=[]
output_data=[]
total_genes=[]

out_folder=os.path.dirname(__file__)+"/"+output_folder.split(":")[1]
if os.path.isdir(out_folder)==False:
    os.mkdir(out_folder)
else:
    shutil.rmtree(out_folder)
    os.mkdir(out_folder)

cluster_threshold=250
for hap in all_haplotype_data:
    print(hap[0])

    for v in hap[1]:

        if v[0]!=data_filter.split(":")[1]:
            continue

        if v[0] not in first_iter:
            first_iter.append(v[0])    
            with open(out_folder+"/d_genes.csv","w") as write:
                writer=csv.writer(write)
                writer.writerow(["Source","GeneType","Contig","Pos","Strand","Sequence","Productive","Locus","Upstream Heptamer","Upstream Nonamer","Downstream Heptamer","Downstream Nonamer","Whole Section","Location Relative to V-Cluster","Heptamer Threshold","Nonamer Threshold"])
                write.close()

        print(v[0])
        for cont in v[1]:
            contig=cont[0]
            print(contig)
            
            for c in cont[1]:
                print(c[0],c[1],len(c[2]))

            baseline=[]
            cluster_genes=[]

            for q in cont[1]:
                gene_num=len(q[2])

                if gene_num!=0:
                    if baseline==[]:
                        for gene in q[2]:
                            if len(gene)==14:
                                gene.append(q[0])
                                gene.append(q[1])
                            baseline.append(gene)
                            cluster_genes.append(gene)
                        break
            print(baseline)
            
            def build_genes(contigs,baseline,cluster_genes):
                for q in contigs:
                    gene_num=len(q[2])
                    if gene_num!=0:
                        for gene in q[2]:
                            if gene[4]=="-":
                                gene_pos=int(gene[3])-int(len(gene[5])/2)
                            elif gene[4]=="+":
                                gene_pos=int(gene[3])+int(len(gene[5])/2)
                            
                            if str(gene[4])=="-":
                                test_gene = copy.deepcopy(gene)
                                test_gene[3] = int(gene[3])-len(gene[5])
                            else:
                                test_gene = copy.deepcopy(gene)
                            for b in baseline:
                                if b[4]=="-":
                                    b_pos=int(b[3])-int(len(b[5])/2)
                                elif b[4]=="+":
                                    b_pos=int(b[3])+int(len(b[5])/2)
                                if abs(gene_pos-b_pos)<=cluster_threshold:
                                    there_already=False
                                    for b1 in cluster_genes:
                                        if str(b1[4])=="-":
                                            test_b = copy.deepcopy(b1)
                                            test_b[3] = int(b1[3])-len(b1[5])
                                        else:
                                            test_b = copy.deepcopy(b1)

                                        if test_b[:5]==test_gene[:5]:
                                            if len(gene[5])>len(b1[5]):
                                                if len(gene)==14:
                                                    gene.append(q[0])
                                                    gene.append(q[1])
                                                b1[:]=gene
                                            there_already=True
                                        
                                        if str(b1[3])==str(gene[3]) and str(b1[4])==str(gene[4]):
                                            if len(gene[5])>len(b1[5]):
                                                if len(gene)==14:
                                                    gene.append(q[0])
                                                    gene.append(q[1])
                                                b1[:]=gene
                                            there_already=True
                                    if there_already==False:
                                        if len(gene)==14:
                                            gene.append(q[0])
                                            gene.append(q[1])
                                        cluster_genes.append(gene)
                
                for cl in cluster_genes:
                    if str(cl[4])=="-":
                        test_cl = copy.deepcopy(cl)
                        test_cl[3] = int(cl[3])-len(cl[5])
                    else:
                        test_cl = copy.deepcopy(cl)
                    there_already=False

                    for b1 in baseline:
                        if str(b1[4])=="-":
                            test_b = copy.deepcopy(b1)
                            test_b[3] = int(b1[3])-len(b1[5])
                        else:
                            test_b = copy.deepcopy(b1)
                        
                        if test_b[:5]==test_cl[:5]:
                            if len(cl[5])>len(b1[5]):
                                b1[:]=cl
                            there_already=True
                    if there_already==False:
                        baseline.append(cl)
                

            end=False
            while end==False:
                old_cluster=copy.deepcopy(cluster_genes)
                build_genes(cont[1],baseline,cluster_genes)
                print("\n",len(cluster_genes),"\n")
                if cluster_genes==old_cluster:
                    end=True
            
            '''keep_cluster_genes=[]
            for cg in cluster_genes:
                found=False
                for cg1 in cluster_genes:
                    if cg[3]==cg1[3] and cg[4]==cg1[4] and cg[5]!=cg1[5]:
                        found=True
                        if len(cg1[5])>len(cg[5]):
                            if cg1 not in keep_cluster_genes:
                                keep_cluster_genes.append(cg1)
                        else:
                            if cg not in keep_cluster_genes:
                                keep_cluster_genes.append(cg)
                if found==False:
                    keep_cluster_genes.append(cg)
            cluster_genes=keep_cluster_genes'''

            with open(out_folder+"/d_genes.csv","a") as write:
                writer=csv.writer(write)
                cluster_genes.sort(key=lambda x: int(x[3]))
                for g in cluster_genes:
                    if g[13].startswith("/"):
                        g[13]==g[13].replace("/","")
                    writer.writerow(g)
                    total_genes.append(g)
                write.close()
            
            if [hap[0],v[0],cont[0],len(cluster_genes),"",""] not in output_data:
                print([hap[0],v[0],cont[0],len(cluster_genes),"",""])
                output_data.append([hap[0],v[0],cont[0],len(cluster_genes),"",""])

with open(out_folder+"/thresh_results.csv","w",newline="") as write:
    writer=csv.writer(write)
    writer.writerow(["Source","Locus","Contig","Number of Genes","Cluster Score","Number of Clusters"])    
    for outd in output_data:
        print(outd)
        writer.writerow(outd)
    write.close()