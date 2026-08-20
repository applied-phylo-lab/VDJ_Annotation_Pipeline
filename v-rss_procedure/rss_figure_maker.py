import os
import sys
import csv
import itertools
import numpy as np
from scipy.stats import entropy
import matplotlib.pyplot as plt
from matplotlib.textpath import TextPath
from matplotlib.font_manager import FontProperties
from matplotlib.patches import PathPatch
from matplotlib.transforms import Affine2D

data_dir=sys.argv[1]
if sys.argv[2]!="-none" and sys.argv[2]!="-def":
    heptamer_threshold=""
    nonamer_threshold=""
    for thr in sys.argv[2].split(":"):
        if "h" in thr:
            heptamer_threshold=int(thr.replace("h",""))
        elif "n" in thr:
            nonamer_threshold=int(thr.replace("n",""))        
    if heptamer_threshold=="" or nonamer_threshold=="":
        sys.exit("Please enter valid RSS analysis threshold argument")

rss_data=[]
with open(data_dir+"/combined_rss_zones.csv","r") as read:
    reader=csv.reader(read)
    header=next(reader)
    for row in reader:
        rss_data.append(row)
    read.close()

gene_data=[]
with open("/local/storage/kav67/clean_birds/gene_list.csv") as read:
    reader=csv.reader(read)
    header=next(reader)
    for row in reader:
        gene_data.append(row)
    read.close()

order_list=[]
for r in rss_data:
    found=False
    for o in order_list:
        if r[0].split("/")[0]==o[0]:
            found=True
            o[1].append(r)
    if found==False:
        order_list.append([r[0].split("/")[0],[r]])

species_list=[]
for r in rss_data:
    found=False
    for s in species_list:
        if r[0].split("/")[1]==s[0]:
            found=True
            s[1].append(r)
    if found==False:
        species_list.append([r[0].split("/")[1],[r]])
    
hap_list=[]
for r in rss_data:
    found=False
    for h in hap_list:
        if r[0]==h[0]:
            found=True
            h[1].append(r)
    if found==False:
        hap_list.append([r[0],[r]])


#total genes vs rss genes scatter plots

order_names=["Cormorants","Cranes","Doves","Eagles","Falcons","Hummingbirds","Ibises","Landfowl","MiscBirds","Owls","Parrots","Plovers","Songbirds","Suboscines","Waterfowl","Woodpeckers"] #names of all bird order folders
clean_birds="/local/storage/kav67/clean_birds"

total_hap_list=[]
for f in os.listdir(clean_birds):
    if f in order_names and os.path.isdir(clean_birds+"/"+f):        
        for f1 in os.listdir(clean_birds+"/"+f):
            if os.path.isdir(clean_birds+"/"+f+"/"+f1) and f1!="patchworkplot":
                for f2 in os.listdir(clean_birds+"/"+f+"/"+f1):
                    if os.path.isdir(clean_birds+"/"+f+"/"+f1+"/"+f2):
                        total_hap_list.append([f+"/"+f1+"/"+f2,["IGH",[],[]],["IGL",[],[]],["TRA",[],[]],["TRB",[],[]],["TRD",[],[]],["TRG",[],[]]])
          
for t in total_hap_list:
    for rss in rss_data:
        if rss[0]==t[0]:
            for locus in t[1:]:
                if locus[0]==rss[7]:
                    locus[1].append(rss)

    for g in gene_data:
        if g[0]==t[0]:
            for locus in t[1:]:
                if locus[0]==g[7]:
                    locus[2].append(g)

loci_plot_data=[]
for t in total_hap_list:
    loci_plot_data.append([t[0],["IGH",0,0],["IGL",0,0],["TRA",0,0],["TRB",0,0],["TRD",0,0],["TRG",0,0]])
    for dat in t[1:]:
        for l in loci_plot_data[-1]:
            if dat[0]==l[0]:
                l[1]=len(dat[1])
                l[2]=len(dat[2])

LOCI = ["IGH", "IGL", "TRA", "TRB", "TRD", "TRG"]

#heptamer haplotype bar chart + meme figure

def stacked_sequence_logo(
    data,
    weight_by="gene_count",
    bases=("A", "C", "G", "T"),
    col_width=0.82,
    total_height=1.0,
    colors=None,
    figsize=None,
    ax=None,
    output_path=None,
):
    """
    Build a classic stacked sequence logo, weighted by a chosen metadata
    column, with every position's stack summing to the same height.
 
    Parameters
    ----------
    data : list of list
        Rows of ``[rss_seq, n_genes, n_orders, n_species, n_haplotypes,
        n_contigs, n_pos_genes, n_neg_genes]``. All ``rss_seq`` must share
        the same length.
    weight_by : {"gene_count", "spec_count", "hap_count", "contig_count"}
        Which metadata column weights each sequence's contribution to the
        per-position base frequencies.
    bases : tuple of str
        Symbol alphabet to tally (default DNA ``A/C/G/T``).
    col_width : float
        Width of each position's column (columns are centered on integers
        0..L-1); keep < 1 to leave a small gap between positions.
    total_height : float
        Height every stack sums to (since proportions sum to 1, every
        column reaches exactly this height).
    colors : dict, optional
        Override base -> color mapping.
    figsize : (float, float), optional
    ax : matplotlib Axes, optional
        Draw into an existing axes instead of creating a new figure.
    output_path : str, optional
        If given, the figure is saved to this path as an SVG (regardless
        of the file extension you pass).
 
    Returns
    -------
    fig, ax, proportions
        ``proportions`` is an (L, len(bases)) array of the weighted base
        proportions used (each row sums to 1), for inspection / reuse.
    """

    _WEIGHT_COLUMNS = {
        "gene_count": 1,    # # of all genes
        "spec_count": 3,    # # of species
        "hap_count": 4,     # # of haplotypes
        "contig_count": 5,  # # of contigs
    }
    
    _BASE_COLORS = {'A': '#F8CD9C', 'C': '#172869', 'G': '#088BBE', 'T': '#EA7580'}

    if weight_by not in _WEIGHT_COLUMNS:
        raise ValueError(
            f"weight_by must be one of {list(_WEIGHT_COLUMNS)}, got {weight_by!r}"
        )
    col = _WEIGHT_COLUMNS[weight_by]
    color_map = {**_BASE_COLORS, **(colors or {})}
 
    seqs = [row[0].upper() for row in data]
    weights = np.array([float(row[col]) for row in data], dtype=float)
 
    lengths = {len(s) for s in seqs}
    if len(lengths) != 1:
        raise ValueError(f"All sequences must be equal length, got lengths {sorted(lengths)}")
    seq_len = lengths.pop()
 
    if np.any(weights < 0):
        raise ValueError("weights must be non-negative")
 
    # -------------------------------------------------- weighted counts
    base_idx = {b: i for i, b in enumerate(bases)}
    counts = np.zeros((seq_len, len(bases)))
    for seq, w in zip(seqs, weights):
        if w == 0:
            continue
        for pos, ch in enumerate(seq):
            i = base_idx.get(ch)
            if i is not None:
                counts[pos, i] += w
 
    totals = counts.sum(axis=1, keepdims=True)
    totals_safe = np.where(totals == 0, 1, totals)
    props = counts / totals_safe  # (seq_len, n_bases), each row sums to 1
 
    # -------------------------------------------------- figure setup
    if ax is None:
        if figsize is None:
            figsize = (max(4, seq_len * 0.55), 4.5)
        fig = plt.figure(figsize=figsize)
        # axes fill the entire figure (no margins) so the letters always
        # scale directly with figsize, rather than being shrunk by a
        # fixed-size layout padding
        ax = fig.add_axes([0, 0, 1, 1])
    else:
        fig = ax.figure
 
    fp = FontProperties(weight="bold")
 
    for pos in range(seq_len):
        p = props[pos]
        # stack smallest -> largest, bottom -> top, so the most frequent
        # base ends up on top (standard sequence-logo convention)
        order = np.argsort(p)
        y0 = 0.0
        for i in order:
            frac = p[i]
            if frac <= 1e-6:
                continue
            base = bases[i]
            letter_h = frac * total_height
 
            tp = TextPath((0, 0), base, size=100, prop=fp)
            verts = tp.vertices
            xmin, ymin = verts.min(axis=0)
            xmax, ymax = verts.max(axis=0)
            w_ = max(xmax - xmin, 1e-6)
            h_ = max(ymax - ymin, 1e-6)
 
            scale_x = col_width / w_
            scale_y = letter_h / h_
            x0 = pos - col_width / 2.0
 
            transform = (
                Affine2D()
                .translate(-xmin, -ymin)
                .scale(scale_x, scale_y)
                .translate(x0, y0)
            )
            patch_path = tp.transformed(transform)
            ax.add_patch(
                PathPatch(patch_path, facecolor=color_map.get(base, "black"),
                          edgecolor="none")
            )
            y0 += letter_h
 
    # -------------------------------------------------- cosmetics
    ax.set_xlim(-0.5, seq_len - 0.5)
    ax.set_ylim(0, total_height)
    ax.axis("off")
 
    if output_path is not None:
        fig.savefig(output_path, format="svg")
 
    return fig, ax, props

def plot_rss_gene_chart(data, output_filename, figure_title, rss_type, figsize=(16, 8)):
    """
    Build a nested-bar chart from RSS-sequence taxonomy statistics and save
    it as SVG.
 
    Each row of `data` must be, in order:
        [rss_seq, n_genes, n_orders, n_species, n_haplotypes, n_contigs,
         n_positive_genes, n_negative_genes]
 
    (n_genes, n_contigs, n_positive_genes, and n_negative_genes are
    accepted for convenience -- e.g. if you're reusing the same data rows
    from elsewhere -- but are not plotted.)
 
    For every rss_seq, a single bar is drawn:
 
      - An outer/background bar for haplotype count.
      - A species-count bar drawn at the same full width directly on top
        of it.
      - An orders-count bar drawn at the same full width on top of that.
 
    Since haplotypes >= species >= orders is the expected containment
    relationship, each larger bar remains visible as a band above the
    next one down. Only the haplotype value is labeled, once per bar.
 
    The chart automatically thins bar borders, shrinks label font size,
    and rotates the x-axis/value labels to vertical as the number of
    rss_seq entries grows, so it stays readable with 50+ entries without
    the figure itself changing size.
 
    Parameters
    ----------
    data : list of lists/tuples
        Rows of [rss_seq, n_genes, n_orders, n_species, n_haplotypes,
        n_contigs, n_positive_genes, n_negative_genes].
    output_filename : str
        Path/filename to save the chart to (should end in ".svg").
    figsize : tuple(float, float), optional
        Figure size in inches, width x height. This is a fixed size that
        does NOT change based on how many rss_seq rows are in `data` --
        pass a larger value here yourself if you want more room for a
        very large dataset. Defaults to (16, 8), a wide/short aspect ratio.
 
    Returns
    -------
    str
        The output_filename the chart was saved to.
    """
    if not data:
        raise ValueError("data must contain at least one row")
 
    rss_seqs = [str(row[0]) for row in data]
    # row[1] = n_genes -> not plotted
    #n_orders = [row[2] for row in data]
    n_species = [row[3] for row in data]
    n_haplotypes = [row[4] for row in data]
    # row[5] = n_contigs -> not plotted
    # row[6] = n_positive_genes -> not plotted
    # row[7] = n_negative_genes -> not plotted
 
    n = len(data)
    x = np.arange(n)
 
    # ---------------- layout ----------------
    # Fraction of each rss_seq slot (width 1.0) filled by the bar; the
    # rest is the gap left between adjacent sequences. Higher = less
    # space between sequences.
    category_fill = 0.85
    bar_width = category_fill
 
    # ---------------- readability scaling ----------------
    # As the number of sequences grows, the bar gets visually narrower
    # even though the figure size stays fixed, so thin out borders/labels
    # and switch to vertical text to avoid clutter/overlap.
    if n <= 15:
        tick_fontsize, value_fontsize, edge_width, rotate_labels = 15, 14, 0.8, False
    elif n <= 30:
        tick_fontsize, value_fontsize, edge_width, rotate_labels = 13, 12, 0.6, False
    elif n <= 60:
        tick_fontsize, value_fontsize, edge_width, rotate_labels = 11, 10, 0.4, True
    else:
        tick_fontsize, value_fontsize, edge_width, rotate_labels = 9, 8, 0.2, True
 
    tick_rotation = 90 if rotate_labels else 45
    tick_ha = "center" if rotate_labels else "right"
 
    # ---------------- colors ----------------
    color_haplotypes = "#EA7580"  # outer/background bar
    color_species = "#172869"     # medium blue
    color_orders = "#F8CD9C"      # dark navy blue
 
    edge_color = "white"
 
    fig, ax = plt.subplots(figsize=figsize)
 
    # Nested bar: haplotypes (back) > species > orders (front), all full width.
    ax.bar(x, n_haplotypes, width=bar_width, color=color_haplotypes,
           edgecolor=edge_color, linewidth=edge_width, zorder=2, label="Haplotypes")
    ax.bar(x, n_species, width=bar_width*0.7, color=color_species,
           edgecolor=edge_color, linewidth=edge_width, zorder=3, label="Species")
 
    # ---------------- value labels ----------------
    tallest = max(n_haplotypes)
    pad = tallest * 0.02
    label_rotation = 90 if rotate_labels else 0
 
    # Haplotype label only, above the top of the bar.
    for xi, h in zip(x, n_haplotypes):
        ax.text(xi, h + pad, f"{h:,}", ha="center", va="bottom",
                 fontsize=value_fontsize, rotation=label_rotation, zorder=6, color="#212121")
 
    # ---------------- axes / labels ----------------
    ax.set_xticks(x)
    ax.set_xticklabels(rss_seqs, rotation=tick_rotation, ha=tick_ha, fontsize=tick_fontsize)
    ax.set_ylabel("Count",fontsize=16)
    ax.set_xlabel(rss_type,fontsize=16)
    ax.set_title(figure_title, fontsize=20)
    ax.set_xlim(-0.5, n - 0.5)
    ax.set_ylim(0, tallest * (1.3 if rotate_labels else 1.15))
    ax.grid(axis="y", linestyle="--", alpha=0.3, zorder=0)
    ax.set_axisbelow(True)
 
    for spine in ("top", "right"):
        ax.spines[spine].set_visible(False)
 
    ax.legend(frameon=False, loc="upper right", fontsize=16)
 
    fig.tight_layout()
    fig.savefig(output_filename, format="svg")
    plt.close(fig)
 
    return output_filename

heptamers=[["IGH",[]],["IGL",[]],["TRA",[]],["TRB",[]],["TRD",[]],["TRG",[]]]
nonamers=[["IGH",[]],["IGL",[]],["TRA",[]],["TRB",[]],["TRD",[]],["TRG",[]]]
for rss in rss_data:
    for h in heptamers:
        if rss[7]==h[0]:
            found=False
            for hep in h[1]:
                if rss[9] == hep[0]:
                    found=True
                    hep[1].append(rss)
                    if str(rss[4])=="+":
                        hep[6].append(rss)
                    elif str(rss[4])=="-":
                        hep[7].append(rss)
            if found==False:
                if str(rss[4])=="+":
                    h[1].append([rss[9],[rss],[],[],[],[],[rss],[]])   
                elif str(rss[4])=="-":
                    h[1].append([rss[9],[rss],[],[],[],[],[],[rss]])               
    
    for n in nonamers:
        if rss[7]==n[0]:
            found=False
            for non in n[1]:
                if rss[10] == non[0]:
                    found=True
                    non[1].append(rss)
                    if str(rss[4])=="+":
                        non[6].append(rss)
                    elif str(rss[4])=="-":
                        non[7].append(rss)
            if found==False:
                if str(rss[4])=="+":
                    n[1].append([rss[10],[rss],[],[],[],[],[rss],[]])   
                elif str(rss[4])=="-":
                    n[1].append([rss[10],[rss],[],[],[],[],[],[rss]])

for locus in heptamers:
    if len(locus[1])>0:
        heptamer_data=[]
        for heptamer in locus[1]:
            for rss in rss_data:
                if heptamer[0] == rss[9] and rss[7] == locus[0]:
                    if rss[0].split("/")[0] not in heptamer[2]:
                        heptamer[2].append(rss[0].split("/")[0])
                    if rss[0].split("/")[1] not in heptamer[3]:
                        heptamer[3].append(rss[0].split("/")[1])
                    if rss[0].split("/")[2] not in heptamer[4]:
                        heptamer[4].append(rss[0].split("/")[2])
                    if rss[2] not in heptamer[5]:
                        heptamer[5].append(rss[2])
        
        write_list=[]
        for heptamer in locus[1]:
            if len(heptamer[4])>heptamer_threshold:
                write_row = [heptamer[0],len(heptamer[1]),len(heptamer[2]),len(heptamer[3]),len(heptamer[4]),len(heptamer[5]),len(heptamer[6]),len(heptamer[7])]
                write_list.append(write_row)
        write_list.sort(key=lambda x: x[4], reverse=True)

        with open(data_dir+"/heptamer_counts.csv","w") as write:
            writer=csv.writer(write)
            writer.writerow(["Heptamer","Occurences","Orders","Species","Haplotypes","Contigs","Positive Genes","Negative Genes"])
            for wr in write_list:
                writer.writerow(wr)
                heptamer_data.append(wr)
            write.close()
        
        n = 0
        for h in heptamer_data:
            n = n+h[1]
        print("Heptamer n = value: ",n)    
        stacked_sequence_logo(heptamer_data, weight_by="hap_count", output_path=data_dir+"/v_heptamer_meme.svg", figsize=(14,4))
        plot_rss_gene_chart(heptamer_data, data_dir+"/v_heptamer_bar_chart.svg", "Number of Unique Haplotypes and Species per V Heptamer", "Heptamers")


for locus in nonamers:
    if len(locus[1])>0:
        nonamer_data=[]
        for nonamer in locus[1]:
            for rss in rss_data:
                if nonamer[0] == rss[10] and rss[7] == locus[0]:
                    if rss[0].split("/")[0] not in nonamer[2]:
                        nonamer[2].append(rss[0].split("/")[0])
                    if rss[0].split("/")[1] not in nonamer[3]:
                        nonamer[3].append(rss[0].split("/")[1])
                    if rss[0].split("/")[2] not in nonamer[4]:
                        nonamer[4].append(rss[0].split("/")[2])
                    if rss[2] not in nonamer[5]:
                        nonamer[5].append(rss[2])
        
        write_list=[]
        for nonamer in locus[1]:
            if len(nonamer[4])>nonamer_threshold:
                write_row = [nonamer[0],len(nonamer[1]),len(nonamer[2]),len(nonamer[3]),len(nonamer[4]),len(nonamer[5]),len(nonamer[6]),len(nonamer[7])]
                write_list.append(write_row)
        write_list.sort(key=lambda x: x[4], reverse=True)

        with open(data_dir+"/nonamer_counts.csv","w") as write:
            writer=csv.writer(write)
            writer.writerow(["Nonamer","Occurences","Orders","Species","Haplotypes","Contigs","Positive Genes","Negative Genes"])
            for wr in write_list:
                writer.writerow(wr)
                nonamer_data.append(wr)
            write.close()
        
        n = 0
        for h in nonamer_data:
            n = n+h[1]
        print("Nonamer n = value: ",n)  
        stacked_sequence_logo(nonamer_data, weight_by="hap_count", output_path=data_dir+"/v_nonamer_meme.svg", figsize=(14,4))
        plot_rss_gene_chart(nonamer_data, data_dir+"/v_nonamer_bar_chart.svg", "Number of Unique Haplotypes and Species per V Nonamer", "Nonamers")