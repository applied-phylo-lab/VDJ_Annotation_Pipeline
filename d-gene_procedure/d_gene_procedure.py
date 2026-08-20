import os
import sys
import subprocess
import csv

hep_thresh_range=sys.argv[1] # ex) -h:1-3
non_thresh_range=sys.argv[2] # ex) -n:2-6
thresh_exclude=sys.argv[3] # ex) -e:3-6,3-5 (-all to use all threshes)
thresh_type=sys.argv[4] # ex) -base:IGH
both_same=sys.argv[5] # ex) -both_same:2-9 (-none to not use)
motif=sys.argv[6] # ex) motif:1-2:prop (-none to not use)
size=sys.argv[7] # ex) -size:2 (-none to not use)
output=sys.argv[8] # ex) -outdel:all_data/all_birds (-none to not save to a folder)
if output!="-none":
    if len(output.split(":")[1].split("/"))!=2:
        sys.exit("Please enter a name for the main output dir and all child dirs (main/child)")
    output=output+"_"+thresh_type.split(":")[1].lower()
data_location=sys.argv[9] # ex) -csv (provide gene_list.csv in input_data)  |  -fold : path to IgDetective data folder with (order --> species --> haplotype) subdir structure
output_folder=sys.argv[10] # ex) -final_out : name of final output folder
figure_thresh=sys.argv[11] # ex) -fig_thresh:hap-0-0
filtering=sys.argv[12] # ex) -f:bStrDea1_pri (-n for none)

if filtering=="-n": 
    threshold_arg = "-none:"
else:
    threshold_arg = "-"+str(filtering.split(":")[1])+":"

threshold_arg = threshold_arg+str(thresh_type.split(":")[1])


d_search_runner = "d_search_runner.py"
threshold_chooser = "threshold_chooser.py"
d_gene_figures = "d_gene_figures.py"

subprocess.run(["python",d_search_runner,hep_thresh_range,non_thresh_range,thresh_exclude,thresh_type,both_same,motif,size,output,data_location,filtering])
subprocess.run(["python",threshold_chooser,threshold_arg,output.split(":")[1].split("/")[0],data_location,output_folder])
subprocess.run(["python",d_gene_figures,output_folder.split(":")[1]+"/d_genes.csv",thresh_type.split(":")[1],figure_thresh,"rss:size:tree"])
