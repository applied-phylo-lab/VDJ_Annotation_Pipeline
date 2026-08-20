import os
import csv
import sys
import statistics
import numpy as np
from scipy.stats import entropy
import matplotlib.pyplot as plt
from matplotlib.textpath import TextPath
from matplotlib.font_manager import FontProperties
from matplotlib.patches import PathPatch
from matplotlib.transforms import Affine2D
from scipy.interpolate import make_interp_spline
from scipy.ndimage import gaussian_filter1d
from scipy.signal import find_peaks
from matplotlib.patches import Patch


data_file=sys.argv[1]
data_folder=os.path.dirname(data_file)
locus=sys.argv[2]
threshes=sys.argv[3]
filter_type=str(threshes.split(":")[1].split("-")[0])
if filter_type=="occur":
    filter_num=1
elif filter_type=="hap":
    filter_num=4
elif filter_type=="spec":
    filter_num=3
elif filter_type=="order":
    filter_num=2
hep_thresh=int(threshes.split(":")[1].split("-")[1].replace("h",""))
non_thresh=int(threshes.split(":")[1].split("-")[2].replace("n",""))


running = sys.argv[4]
to_run=[]
for r in running.split(":"):
    to_run.append(r)

data=[]
with open(data_file,"r") as read:
    reader=csv.reader(read)
    header=next(reader)
    for row in reader:
        data.append(row)
    read.close()

contigs=[]
with open("../input_data/contig_list.csv","r") as read:
    reader=csv.reader(read)
    header=next(reader)
    for row in reader:
        if str(row[2])==locus:
            contigs.append(row)
    read.close()

genomes=[]
with open("../input_data/genome_paths.csv","r") as read:
    reader=csv.reader(read)
    header=next(reader)
    for row in reader:
        genomes.append(row)
    read.close()

def rss_fill(rss, gene, output):
    #yeilds [rss seq, list of all genes, list of order names, list of species names, list of haplotype names, list of contig names, list of postive genes, list of negative genes]
    found=False
    for o in output:
        if o[0]==rss:
            found=True
            if gene not in o[1]:
                o[1].append(gene)
            if gene[0].split("/")[0] not in o[2]:
                o[2].append(gene[0].split("/")[0])
            if gene[0].replace("/"+gene[0].split("/")[-1],"") not in o[3]:
                o[3].append(gene[0].replace("/"+gene[0].split("/")[-1],""))
            if gene[0] not in o[4]:
                o[4].append(gene[0])  
            if gene[2] not in o[5]:
                o[5].append(gene[2])
            
            if str(gene[4])=="+":
                o[6].append(gene)
            elif str(gene[4])=="-":
                o[7].append(gene)
    if found==False:
        if str(gene[4])=="+":
            output.append([rss, [gene], [gene[0].split("/")[0]], [gene[0].replace("/"+gene[0].split("/")[-1],"")], [gene[0]], [gene[2]], [gene], []])
        elif str(gene[4])=="-":
            output.append([rss, [gene], [gene[0].split("/")[0]], [gene[0].replace("/"+gene[0].split("/")[-1],"")], [gene[0]], [gene[2]], [], [gene]])

def rss_interpret(rss_list, rss_size):
    print("\n")
    frequency_matrix=[]
    for rs in range(rss_size):
        frequency_matrix.append([["A",0],["C",0],["T",0],["G",0]])
    rss_count=0
    for r in rss_list:
        r[0]=r[0].upper()
        rss_count=rss_count+int(len(r[1]))
        n=0
        for char in r[0]:
            for base in frequency_matrix[n]:
                if base[0]==char:
                    base[1]+=(int(len(r[1])))
            n+=1
        #print(r[0],len(r[1]))
    if rss_count>0:
        true_frequency_matrix=[]
        for fre in frequency_matrix:
            freq_vector=[]
            for base in fre:
                base.append(round(int(base[1])/rss_count,5))
                freq_vector.append(round(int(base[1])/rss_count,5))
            true_frequency_matrix.append(freq_vector)

        #print(true_frequency_matrix,"\n")
        shannon_entropies = entropy(true_frequency_matrix, base=2, axis=1)
        #for se in shannon_entropies:
        #    print(se)
        mean_shannon = round(sum(shannon_entropies) / len(shannon_entropies),3)
        stdev_shannon = round(statistics.pstdev(shannon_entropies),3)

        print("Mean Shannon: ",mean_shannon,"\nStdev Shannon: ",stdev_shannon)
        return mean_shannon, stdev_shannon
    else:
        return 0, 0

def rss_list_convert(rss_list,rss_filter_thresh,filter_type=4):
    new_rss_list=[]
    for r in rss_list:
        temp_list=[]
        for r1 in r:
            if isinstance(r1,list):
                temp_list.append(len(r1))
            else:
                temp_list.append(r1)
        if temp_list[filter_type]>=rss_filter_thresh:
            new_rss_list.append(temp_list)
    cont=0
    for n in new_rss_list:
        cont=cont+int(n[1])
    print(cont)
    return new_rss_list
 
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

import matplotlib.pyplot as plt
from matplotlib.patches import Patch
import numpy as np

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

def plot_rss_gene_pi_chart(
    data, output_filename, figure_title, rss_type, figsize=(16, 6)
):
    """Build a horizontal row of full pie charts from RSS-sequence statistics

    and save as an SVG.

    Each row of `data` must be, in order:
        [rss_seq, n_genes, n_orders, n_species, n_haplotypes, n_contigs,
         n_positive_genes, n_negative_genes, count_33_36_pm5, count_53_pm5, count_65_67_pm5,
         n_other_genes, count_33_36_pm1, count_53_pm1, count_65_67_pm1]

    Parameters
    ----------
    data : list of lists/tuples
        Sequence taxonomy and region count rows.
    output_filename : str
        Path/filename to save the SVG chart to.
    figure_title : str
        Main title for the figure.
    rss_type : str
        Context label (e.g., '12-RSS' or '23-RSS').
    figsize : tuple(float, float), optional
        Base figure size. Auto-scales width horizontally based on sequence count if using default.

    Returns
    -------
    str
        The output_filename the chart was saved to.
    """
    if not data:
        raise ValueError("data must contain at least one row")

    n = len(data)

    # Automatically expand figure width for large sequence sets if standard default is used
    if figsize == (16, 6) or figsize == (16, 8):
        figsize = (max(16, n * 2.2), 5.5)

    # ---------------- Data Extraction ----------------
    rss_seqs = [str(row[0]) for row in data]
    n_genes = [row[1] for row in data]

    range_33_36_pm5 = [row[8] for row in data]
    range_53_pm5 = [row[9] for row in data]
    range_65_67_pm5 = [row[10] for row in data]
    other_genes = [row[11] for row in data]

    sub_33_36_pm1 = [row[12] for row in data]
    sub_53_pm1 = [row[13] for row in data]
    sub_65_67_pm1 = [row[14] for row in data]

    # ---------------- Color Palette ----------------
    # 33-36 Pair (Blue)
    color_33_pm1 = "#1F4E78"  # Dark Blue (±1)
    color_33_rest = "#5B9BD5"  # Main Blue (Rest of ±5)

    # 53 Pair (Green)
    color_53_pm1 = "#375623"  # Dark Green (±1)
    color_53_rest = "#70AD47"  # Main Green (Rest of ±5)

    # 65-67 Pair (Orange)
    color_65_pm1 = "#833C0C"  # Dark Orange (±1)
    color_65_rest = "#ED7D31"  # Main Orange (Rest of ±5)

    # Other Genes
    color_other = "#8064A2"  # Purple

    slice_colors = [
        color_33_pm1,
        color_33_rest,
        color_53_pm1,
        color_53_rest,
        color_65_pm1,
        color_65_rest,
        color_other,
    ]

    # ---------------- Plot Initialization ----------------
    fig, axes = plt.subplots(1, n, figsize=figsize, squeeze=False)
    axes = axes[0]  # Flatten to 1D subplot array

    for i, ax in enumerate(axes):
        seq = rss_seqs[i]
        total = n_genes[i]

        v33_5 = range_33_36_pm5[i]
        v53_5 = range_53_pm5[i]
        v65_5 = range_65_67_pm5[i]
        v_oth = other_genes[i]

        v33_1 = sub_33_36_pm1[i]
        v53_1 = sub_53_pm1[i]
        v65_1 = sub_65_67_pm1[i]

        # Full pie values: [±1 subcount, remaining ±5 count] for each region
        values = [
            v33_1,
            max(0, v33_5 - v33_1),
            v53_1,
            max(0, v53_5 - v53_1),
            v65_1,
            max(0, v65_5 - v65_1),
            v_oth,
        ]

        if sum(values) > 0:
            # Full Pie Chart (Solid wedges filling to the center)
            ax.pie(
                values,
                radius=1.0,
                colors=slice_colors,
                wedgeprops=dict(edgecolor="white", linewidth=1.2),
                startangle=90,
            )
        else:
            # Placeholder circle if entry has zero total counts
            circle = plt.Circle((0, 0), 1.0, color="#E0E0E0", fill=True)
            ax.add_artist(circle)
            ax.text(
                0,
                0,
                "No Data",
                ha="center",
                va="center",
                fontsize=9,
                color="#666666",
            )

        # Title for each pie chart
        ax.set_title(
            f"{seq}\n(n={total:,})",
            fontsize=10,
            fontweight="bold",
            pad=10,
            color="#212121",
        )

    # ---------------- Global Legend & Aesthetics ----------------
    legend_elements = [
        Patch(facecolor=color_33_rest, edgecolor="w", label="33-36 ±5"),
        Patch(facecolor=color_33_pm1, edgecolor="w", label="33-36 ±1"),
        Patch(facecolor=color_53_rest, edgecolor="w", label="53 ±5"),
        Patch(facecolor=color_53_pm1, edgecolor="w", label="53 ±1"),
        Patch(facecolor=color_65_rest, edgecolor="w", label="65-67 ±5"),
        Patch(facecolor=color_65_pm1, edgecolor="w", label="65-67 ±1"),
        Patch(facecolor=color_other, edgecolor="w", label="Other Genes"),
    ]

    fig.suptitle(f"{figure_title} ({rss_type})", fontsize=18, y=1.05)

    # Place figure-level legend underneath the row of pie charts
    fig.legend(
        handles=legend_elements,
        loc="upper center",
        bbox_to_anchor=(0.5, -0.02),
        ncol=7,
        frameon=False,
        fontsize=11,
    )

    fig.tight_layout()
    fig.savefig(output_filename, format="svg", bbox_inches="tight")
    plt.close(fig)

    return output_filename


if "size" in to_run:
    #d gene size histogram

    def d_gene_smoothed_violin(
        data,
        peaks=None,  # Pass list like ["10-15", "20-25", "30"] or None for auto-detection
        title="D Gene Size Distribution",
        ylabel="D Gene Size (bp)",
        color="#3498db",
        edgecolor="#2c3e50",
        line_color="black",  # Black peak lines and arrows
        max_width=0.4,
        smoothing_sigma=0.8,
        top_n_peaks=3,
        figsize=(7, 9),
        save_path=None,
        show=False,
        ):
        """Smooth violin plot with peak indicator lines, count percentages for single sizes or ranges,
        tight margins, and a centered title relative to the violin.
        """
        # Sort data by gene size
        order = np.argsort([row[0] for row in data])
        sizes = np.array([data[i][0] for i in order], dtype=float)
        counts = np.array([data[i][1] for i in order], dtype=float)
        total_count = np.sum(counts)

        # 1. Smooth counts for peak detection and contour fitting
        counts_smoothed = (
            gaussian_filter1d(counts, sigma=smoothing_sigma)
            if smoothing_sigma > 0
            else counts
        )

        # 2. Dense grid for smooth violin body
        dense_sizes = np.linspace(sizes.min(), sizes.max(), 300)
        spline = make_interp_spline(sizes, counts_smoothed, k=3)
        dense_counts = spline(dense_sizes)
        dense_counts = np.clip(dense_counts, 0, None)

        x_max_violin = max_width / 2
        scaled_counts = (dense_counts / dense_counts.max()) * x_max_violin

        fig, ax = plt.subplots(figsize=figsize)

        # Plot symmetrical violin outline
        ax.fill_betweenx(
            dense_sizes,
            -scaled_counts,
            scaled_counts,
            color=color,
            edgecolor=edgecolor,
            linewidth=1.2,
            alpha=0.8,
        )

        # 3. Process Peak Inputs (or fallback to auto-detection)
        if peaks is None:
            peak_indices, _ = find_peaks(counts_smoothed)
            if len(peak_indices) < top_n_peaks:
                peak_indices = np.argsort(counts)[-top_n_peaks:]
            else:
                peak_indices = peak_indices[
                    np.argsort(counts[peak_indices])[-top_n_peaks:]
                ]
            peak_indices = peak_indices[np.argsort(sizes[peak_indices])]
            peaks_list = [str(int(sizes[idx])) for idx in peak_indices]
        else:
            peaks_list = list(peaks)

        # Tighter x-axis limits to accommodate range labels cleanly
        x_min_bound = -x_max_violin - 0.03
        x_max_bound = x_max_violin + 0.20

        # 4. Draw Black Lines & Percentage Annotations
        for p_item in peaks_list:
            p_str = str(p_item).strip().replace("–", "-").replace("—", "-")

            # Parse range ("num-num") or single base pair ("num")
            if "-" in p_str and not p_str.startswith("-"):
                parts = p_str.split("-")
                low = float(parts[0].strip())
                high = float(parts[1].strip())
            else:
                low = high = float(p_str)

            # Filter data points within target range/size
            if low == high:
                mask = np.isclose(sizes, low)
            else:
                mask = (sizes >= low) & (sizes <= high)

            p_count = np.sum(counts[mask])
            p_pct = (p_count / total_count * 100) if total_count > 0 else 0.0

            # Determine y position for drawing horizontal line and arrow anchor
            if low == high:
                y_pos = low
            else:
                if np.any(mask):
                    # Anchor at max count size within the requested range
                    y_pos = sizes[mask][np.argmax(counts[mask])]
                else:
                    y_pos = (low + high) / 2.0

            # Calculate violin width at y_pos
            y_pos_clamped = np.clip(y_pos, sizes.min(), sizes.max())
            count_at_y = max(0.0, float(spline(y_pos_clamped)))
            p_width = (count_at_y / dense_counts.max()) * x_max_violin

            # Black dashed horizontal line across violin body width
            ax.hlines(
                y=y_pos,
                xmin=-p_width,
                xmax=p_width,
                color=line_color,
                linewidth=1.2,
                linestyle="--",
            )

            # Label formatted with exact input label, total range count, and percentage
            label_text = f"{p_str} bp\nn = {int(round(p_count)):,} ({p_pct:.1f}%)"
            ax.annotate(
                label_text,
                xy=(p_width, y_pos),
                xytext=(x_max_violin + 0.03, y_pos),
                arrowprops=dict(
                    arrowstyle="->",
                    color=line_color,
                    lw=1,
                    shrinkA=0,
                    shrinkB=3,
                ),
                va="center",
                fontsize=9,
                fontweight="bold",
                color="#2c3e50",
            )

        ax.set_ylabel(ylabel, fontsize=12)
        ax.set_xlim(x_min_bound, x_max_bound)
        ax.set_ylim(sizes.min() - 1, sizes.max() + 1)
        ax.set_xticks([])

        # Center title directly over x=0 (middle of the violin plot)
        violin_center_axes = (0 - x_min_bound) / (x_max_bound - x_min_bound)
        ax.set_title(
            title, fontsize=14, fontweight="bold", pad=12, x=violin_center_axes
        )

        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)
        ax.spines["bottom"].set_visible(False)

        fig.tight_layout()
        if save_path:
            fig.savefig(save_path, dpi=150, bbox_inches="tight")
        if show:
            plt.show()

        return fig, ax

    d_gene_sizes=[]
    for d in data:
        found=False
        for s in d_gene_sizes:
            if s[0]==len(d[5]):
                found=True
                s[1]+=1
        if found==False:
            d_gene_sizes.append([len(d[5]),1])
    d_gene_sizes.sort(key=lambda x: x[0])
    n=0
    true_d_gene_sizes=[]
    for d in d_gene_sizes:
        if len(d_gene_sizes) < (n+1):
            if d_gene_sizes[n+1][0]-d_gene_sizes[n][0]!=1:
                for q in range(d_gene_sizes[n+1][0]-d_gene_sizes[n][0]):
                    if d_gene_sizes[n][0]+q == d_gene_sizes[n][0]:
                        true_d_gene_sizes.append(d_gene_sizes[n])
                    else:
                        true_d_gene_sizes.append([d_gene_sizes[n][0]+q, 0])
            else:
                true_d_gene_sizes.append(d_gene_sizes[n])
        else:
            true_d_gene_sizes.append(d_gene_sizes[n])
        n+=1

    for td in true_d_gene_sizes:
        print(td)
    igh_peaks=["33-36","53","65-67"]
    trd_peaks=["13-14","33-36"]
    d_gene_smoothed_violin(true_d_gene_sizes, save_path=data_folder+"/d_gene_lengths_violin.svg", peaks=trd_peaks)


if "rss" in to_run:
    #d rss figures (up and down)    

    # motifs, top hep and non

    up_haps=[]
    down_haps=[]
    up_nons=[]
    down_nons=[]
    for gene in data:
        rss_fill(gene[8], gene, up_haps)
        rss_fill(gene[9], gene, up_nons)
        rss_fill(gene[10], gene, down_haps)
        rss_fill(gene[11], gene, down_nons)

    up_haps.sort(key=lambda x: len(x[filter_num]), reverse=True)
    down_haps.sort(key=lambda x: len(x[filter_num]), reverse=True)
    up_nons.sort(key=lambda x: len(x[filter_num]), reverse=True)
    down_nons.sort(key=lambda x: len(x[filter_num]), reverse=True)         

    print("\nup hap")
    input_up_haps = rss_list_convert(up_haps, hep_thresh)
    print("\ndown hap")
    input_down_haps = rss_list_convert(down_haps, hep_thresh)
    print("\nup non")
    input_up_nons = rss_list_convert(up_nons, non_thresh)
    print("\ndown non")
    input_down_nons = rss_list_convert(down_nons, non_thresh)

    stacked_sequence_logo(input_up_haps, weight_by="hap_count", output_path=data_folder+"/upstream_heptamer_meme.svg", figsize=(14,4))
    stacked_sequence_logo(input_up_nons, weight_by="hap_count", output_path=data_folder+"/upstream_nonamer_meme.svg", figsize=(14,4))
    stacked_sequence_logo(input_down_haps, weight_by="hap_count", output_path=data_folder+"/downstream_heptamer_meme.svg", figsize=(14,4))
    stacked_sequence_logo(input_down_nons, weight_by="hap_count", output_path=data_folder+"/downstream_nonamer_meme.svg", figsize=(14,4))

    def gene_size_acess(rss_list):
        for r in rss_list:
            first_tertile=0
            second_tertile=0
            third_tertile=0
            sixsix=0
            fivethree=0
            threefour=0
            other=0
            for gene in r[1]:
                if 41>=len(gene[5])>=28:
                    first_tertile+=1
                    if 33<=len(gene[5])<=36:
                        threefour+=1
                elif 58>=len(gene[5])>=48:
                    second_tertile+=1
                    if len(gene[5])==53:
                        fivethree+=1
                elif 70>=len(gene[5])>=60:
                    third_tertile+=1
                    if 65<=len(gene[5])<=67:
                        sixsix+=1
                else:
                    other+=1
            r.append(first_tertile)
            r.append(second_tertile)
            r.append(third_tertile)
            r.append(other)
            r.append(threefour)
            r.append(fivethree)
            r.append(sixsix)
    
    gene_size_acess(up_haps)
    gene_size_acess(down_haps)
    gene_size_acess(up_nons)
    gene_size_acess(down_nons)
    bar_up_haps = rss_list_convert(up_haps, hep_thresh, filter_type=filter_num)
    bar_down_haps = rss_list_convert(down_haps, hep_thresh, filter_type=filter_num)
    bar_up_nons = rss_list_convert(up_nons, non_thresh, filter_type=filter_num)
    bar_down_nons = rss_list_convert(down_nons, non_thresh, filter_type=filter_num)

    bar_up_haps.sort(key=lambda x: x[filter_num], reverse=True)
    bar_up_nons.sort(key=lambda x: x[filter_num], reverse=True)
    bar_down_haps.sort(key=lambda x: x[filter_num], reverse=True)
    bar_down_nons.sort(key=lambda x: x[filter_num], reverse=True)

    plot_rss_gene_chart(bar_up_haps, data_folder+"/upstream_heptamer_bar_chart.svg", "Number of Unique Haplotypes and Species per Upstream D Heptamer", "Heptamers")
    plot_rss_gene_chart(bar_down_haps, data_folder+"/downstream_heptamer_bar_chart.svg", "Number of Unique Haplotypes and Species per Downstream D Heptamer", "Heptamers")
    plot_rss_gene_chart(bar_up_nons, data_folder+"/upstream_nonamer_bar_chart.svg", "Number of Unique Haplotypes and Species per Upstream D Nonamer", "Nonamers")
    plot_rss_gene_chart(bar_down_nons, data_folder+"/downstream_nonamer_bar_chart.svg", "Number of Unique Haplotypes and Species per Downstream D Nonamer", "Nonamers")

    '''
    # varience within orders
    # uses mean postional shannon entropy across seqs
    order=[]
    for gene in data:
        found=False
        for ordr in order:
            if ordr[0]==gene[0].split("/")[0]:
                found=True
                ordr[1].append(gene)
        if found==False:
            order.append([gene[0].split("/")[0],[gene]])
    
    order_data=[]
    for c in order:
        if c[0]=="MiscBirds":
            continue
        print(c[0],len(c[1]))
        or_up_haps=[]
        or_down_haps=[]
        or_up_nons=[]
        or_down_nons=[]
        
        haplos=[]
        spieces=[]
        conts=[]
        for cgene in c[1]:
            if cgene[0].split("/")[1] not in spieces:
                spieces.append(cgene[0].split("/")[1])
            if cgene[0].split("/")[2] not in haplos:
                haplos.append(cgene[0].split("/")[2])
            if cgene[2] not in conts:
                conts.append(cgene[2])

            rss_fill(cgene[8], cgene, or_up_haps)
            rss_fill(cgene[9], cgene, or_up_nons)
            rss_fill(cgene[10], cgene, or_down_haps)
            rss_fill(cgene[11], cgene, or_down_nons)    
        
        or_up_haps.sort(key=lambda x: len(x[1]), reverse=True)
        or_down_haps.sort(key=lambda x: len(x[1]), reverse=True)
        or_up_nons.sort(key=lambda x: len(x[1]), reverse=True)
        or_down_nons.sort(key=lambda x: len(x[1]), reverse=True)    

        uh_mean_shannon, uh_stdev_shannon = rss_interpret(or_up_haps, 7)
        dh_mean_shannon, dh_stdev_shannon = rss_interpret(or_down_haps, 7)
        un_mean_shannon, un_stdev_shannon = rss_interpret(or_up_nons, 9)
        dn_mean_shannon, dn_stdev_shannon = rss_interpret(or_down_nons, 9)

        or_input_up_haps = rss_list_convert(or_up_haps)
        or_input_down_haps = rss_list_convert(or_down_haps)
        or_input_up_nons = rss_list_convert(or_up_nons)
        or_input_down_nons = rss_list_convert(or_down_nons)
        
        if os.path.isdir(data_folder+"/orders")==False:
            os.mkdir(data_folder+"/orders")
        if os.path.isdir(data_folder+"/orders/"+c[0])==False:
            os.mkdir(data_folder+"/orders/"+c[0])
        stacked_sequence_logo(or_input_up_haps, weight_by="hap_count", output_path=data_folder+"/orders/"+c[0]+"/upstream_heptamer_meme.svg", figsize=(3.5,1))
        stacked_sequence_logo(or_input_up_nons, weight_by="hap_count", output_path=data_folder+"/orders/"+c[0]+"/upstream_nonamer_meme.svg", figsize=(3.5,1))
        stacked_sequence_logo(or_input_down_haps, weight_by="hap_count", output_path=data_folder+"/orders/"+c[0]+"/downstream_heptamer_meme.svg", figsize=(3.5,1))
        stacked_sequence_logo(or_input_down_nons, weight_by="hap_count", output_path=data_folder+"/orders/"+c[0]+"/downstream_nonamer_meme.svg", figsize=(3.5,1))

        order_data.append([c[0],len(c[1]),len(conts),len(haplos),len(spieces),["UH-"+str(uh_mean_shannon),"DH-"+str(dh_mean_shannon),"UN-"+str(un_mean_shannon),"DN-"+str(dn_mean_shannon)],["UH-"+str(data_folder+"/"+c[0]+"/upstream_heptamer_meme.svg"),"DH-"+str(data_folder+"/"+c[0]+"/downstream_heptamer_meme.svg"),"UN-"+str(data_folder+"/"+c[0]+"/upstream_nonamer_meme.svg"),"DN-"+str(data_folder+"/"+c[0]+"/downstream_nonamer_meme.svg")]])
    for od in order_data:
        print(od)
    sys.exit(1)
    '''

if "tree" in to_run:
    #make a phylo tree to show metrics for haplotypes ( number of genes, *rss varience*, location rel to V, *something with size?*, strand distribution / same postion opposite strand )

    spec_data=[]
    for gen in genomes:
        found=False
        for spec in spec_data:
            if str(gen[0].replace("/"+gen[0].split("/")[-1],"")) == spec[0]:
                found=True
        if found==False:
            spec_data.append([str(gen[0].replace("/"+gen[0].split("/")[-1],"")),[]])

    for gen in genomes:
        for spec in spec_data:
            if str(gen[0].replace("/"+gen[0].split("/")[-1],""))==spec[0]:
                found=False
                for hap in spec[1]:
                    if str(gen[0])==hap[0]:
                        found=True
                if found==False:
                    spec[1].append([str(gen[0]),[]])

    for con in contigs:
        for spec in spec_data:
            if str(con[1].replace("/"+con[1].split("/")[-1],""))==spec[0]:
                for hap in spec[1]:
                    if str(con[1])==hap[0]:
                        found=False
                        for cont in hap[1]:
                            if str(con[0])==cont[0]:
                                found=True
                        if found==False:
                            hap[1].append([str(con[0]),[]])

    for d in data:
        for spec in spec_data:
            if d[0].replace("/"+d[0].split("/")[-1],"")==spec[0]:
                for hap in spec[1]:
                    if d[0]==hap[0]:
                        for cont in hap[1]:
                            if d[2]==cont[0]:
                                cont[1].append(d)

    specs=[]

    single_data=[]
    multi_contig=0
    no_contig=0
    single_contig=0
    total=0
    for spec in spec_data:
        specs.append(spec[0])
        #print("\n"+spec[0])
        for hap in spec[1]:
            total+=1
            #print(hap[0].split("/")[-1])

            if len(hap[1])>1:
                multi_contig+=1
            elif len(hap[1])==0:
                no_contig+=1
            elif len(hap[1])==1:
                single_contig+=1
                data_append=[hap[0],hap[1][0][0],[]]
                num_neg=0
                num_pos=0
                for cont_gene in hap[1][0][1]:
                    data_append[-1].append(cont_gene)
                    if str(cont_gene[4])=="-":
                        num_neg+=1
                    elif str(cont_gene[4])=="+":
                        num_pos+=1
                
                pairs=0
                n=0
                for h1 in hap[1][0][1]:
                    if n+1!=len(hap[1][0][1]) and hap[1][0][1][n][4]!=hap[1][0][1][n+1][4]:
                        if hap[1][0][1][n][4]=="+":
                            current=int(hap[1][0][1][n][3])
                            next_gene=int(hap[1][0][1][n+1][3])-len(hap[1][0][1][n+1][5])
                        elif hap[1][0][1][n][4]=="-":
                            current=int(hap[1][0][1][n][3])-len(hap[1][0][1][n][5])
                            next_gene=int(hap[1][0][1][n+1][3])
                        
                        if current==next_gene:
                            pairs+=1
                    n+=1
                
                data_append.append(num_pos)
                data_append.append(num_neg)
                data_append.append(pairs*2)
                single_data.append(data_append)

    print("\nTotal: ",total,"\nSingle: ",single_contig,"\nMulti: ",multi_contig,"\nNone: ",no_contig)

    '''
    with open("all_species.csv","w",newline="") as write:
        writer=csv.writer(write)
        writer.writerow(["Folder Name","English Name","Scientific Name"])
        for sp in specs:
            spec=sp.split("/")[1].replace("_"," ")
            specr=""
            n1=0
            for sq in spec.split(" "):
                n=0
                sqr=""
                for char in sq:
                    if char.isupper() and n!=0:
                        char=" "+char.lower()
                    sqr=sqr+char
                    n+=1
                if n1!=0:
                    specr=specr+sqr.lower()+" "      
                else:
                    specr=specr+sqr+" "
                n1+=1   
            specr=specr[:-1]         
            writer.writerow([sp.split("/")[1],specr,""])
        write.close()
    '''

    input_data=[]
    for sin in single_data:
        dstream=0
        ustream=0
        vcluster=0
        for gene in sin[2]:
            if gene[13]=="downstream":
                dstream+=1
            elif gene[13]=="upstream":
                ustream+=1
            elif gene[13]=="v_cluster":
                vcluster+=1

        if dstream>ustream and dstream>vcluster:
            rel_location="downstream"
        elif ustream>dstream and ustream>vcluster:
            rel_location="upstream"
        elif vcluster>ustream and vcluster>dstream:
            rel_location="v_cluster"
        else:
            rel_location="no_majority"
        
        con_up_haps=[]
        con_down_haps=[]
        con_up_nons=[]
        con_down_nons=[]

        for cgene in sin[2]:
            rss_fill(cgene[8], cgene, con_up_haps)
            rss_fill(cgene[9], cgene, con_up_nons)
            rss_fill(cgene[10], cgene, con_down_haps)
            rss_fill(cgene[11], cgene, con_down_nons)    

        uh_mean_shannon, uh_stdev_shannon = rss_interpret(con_up_haps, 7)
        dh_mean_shannon, dh_stdev_shannon = rss_interpret(con_down_haps, 7)
        un_mean_shannon, un_stdev_shannon = rss_interpret(con_up_nons, 9)
        dn_mean_shannon, dn_stdev_shannon = rss_interpret(con_down_nons, 9)
                
        input_data.append([sin[0],sin[1],len(sin[2]),sin[3],sin[4],sin[5],float(round(((uh_mean_shannon+dh_mean_shannon)/2),3)),float(round(((un_mean_shannon+dn_mean_shannon)/2),3)),rel_location])

    species_names=[]
    with open("../input_data/species_names.csv","r") as read:
        reader=csv.reader(read)
        header=next(reader)
        for row in reader:
            species_names.append(row)
        read.close()

    for i in input_data:
        for sn in species_names:
            if str(sn[0]) in str(i[0]):
                i.insert(1,sn[2])
                break
        source=i[0]
        i.pop(0)
        i.insert(0,source.split("/")[0])
        i.insert(1,source.split("/")[1])
        i.insert(2,source.split("/")[2])
    
    specs=[]
    for i in input_data:
        found=False
        for s in specs:
            if i[3]==s[0]:
                s[1]+=1
                found=True
        if found==False:
            specs.append([i[3],1])
    
    keep_data=[]
    for i in input_data:
        for s in specs:
            if i[3]==s[0]:
                if s[1]<=2:
                    keep_data.append(i)
                break

    with open(data_folder+"/d_tree_data.csv","w",newline="") as write:
        writer=csv.writer(write)
        writer.writerow(["Order","Species","Haplotype","Sci Name","Contig","Num Genes","Num Positive Genes","Num Negative Genes","Num Double Genes","Heptamer Shannon Entropy","Nonamer Shannon Entropy","Position Relative to V's"])
        for i in keep_data:
            writer.writerow(i)
    