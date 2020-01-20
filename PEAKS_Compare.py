import os
import re
import shutil
import pandas as pd
from venn import venn, pseudovenn
import upsetplot
import matplotlib.pyplot as plt
plt.style.use('seaborn')

all_files = []
labels = []
total_peptides = []
unique_peptides = []
peptide_dict = {}
protein_dict = {}
total_proteins = []
unique_proteins = []


def p_diff(v1, v2):
    num = abs(v1 - v2)
    denom = abs(((v1 + v2) / 2))
    return num/denom * 100


def get_name_files(directory="."):
    for root, dirs, files in os.walk(directory, topdown=True):
        for name in files:
            file_path = os.path.join(root, name)
            if file_path.endswith("-peptides.csv"):
                global all_files
                all_files.append(file_path)

    i = 0
    while i < len(all_files):
        label = input("Please provide a name for the sample "+all_files[i]+"\t")
        global labels
        labels.append(label)
        i += 1


def get_peptides(file_list, label_list):
    for i in range(len(file_list)):
        j = i
        file = file_list[i]
        label = label_list[i]
        df = pd.read_csv(file)
        peptides = df.Peptide.tolist()
        global total_peptides
        total_peptides.append(len(peptides))

        peptides = df[df.Unique == "Y"].Peptide.tolist()
        global unique_peptides
        unique_peptides.append(len(set(peptides)))

        global peptide_dict
        peptide_dict[label] = set(peptides)


def get_proteins(file_list, label_list):
    for i in range(len(file_list)):
        file = file_list[i]
        label = label_list[i]
        df = pd.read_csv(file)
        total_prot = df["Protein Accession"].tolist()
        df = df[df.Unique == "Y"].drop_duplicates(subset="Protein Accession")
        proteins = df["Protein Accession"].tolist()

        global protein_dict
        protein_dict[label] = set(proteins)
        global total_proteins
        total_proteins.append(len(set(total_prot)))
        global unique_proteins
        unique_proteins.append(len(list(set(proteins))))


def plot_peptides(total_peptides, unique_peptides, labels):
    fig, ax = plt.subplots(figsize=(10, 5))

    plot_df = pd.DataFrame({"Total Peptides": total_peptides,
                            "Unique Peptides": unique_peptides}, index=labels)
    plot_df.plot.bar(rot=0, fontsize=13, ax=ax)
    ax.legend(loc='best', bbox_to_anchor=[1, 1], fontsize=13)
    ax.set_title("Total and Unique Peptides", fontsize=20)
    plt.savefig("Total_and_unique_peptides.svg")
    plt.savefig("Total_and_unique_peptides.png", bbox_inches="tight")


def plot_peptide_overlap(peptide_dict, labels):
    fig, ax = plt.subplots(figsize=(10, 10))
    if len(peptide_dict) in range(2, 6):
        venn(peptide_dict, cmap="viridis", ax=ax)
    elif len(peptide_dict) == 6:
        pseudovenn(peptide_dict, cmap="viridis", ax=ax)
    else:
        print("No Peptide Venn Diagram plotted due to invalid number of samples.")
        print("Venn Diagrams require between 2 and 6 samples.")

    ax.set_title("Peptide Overlap", fontsize=20)
    ax.legend(labels=labels, fontsize=15, loc='best', bbox_to_anchor=[1.1, 1])
    fig.savefig("Peptide_overlap_venn.svg")
    fig.savefig("Peptide_overlap_venn.png", bbox_inches="tight")


def plot_peptide_upset(peptide_dict):
    color = '#21918cff'
    plot_df = upsetplot.from_contents(peptide_dict)
    upsetplot.plot(plot_df, sort_by='cardinality', subset_size='auto', facecolor=color)
    # plt.ylim(0, 400)
    plt.title("Distribution of Peptide Overlap")

    plt.savefig("Peptide_upset.svg")
    plt.savefig("Peptide_upset.png")


def plot_proteins(unique_proteins, labels):
    fig, ax = plt.subplots(figsize=(10, 5))
    plot_df = pd.DataFrame({"Total Proteins": total_proteins,
                            "Unique Proteins": unique_proteins}, index=labels)
    plot_df.plot.bar(rot=0, fontsize=13, legend=None, ax=ax)

    y_max, difference = plt.yticks()[0][-1], plt.yticks()[0][-1]-plt.yticks()[0][-2]

    ax.set_ylim(0, y_max+difference)
    ax.set_title("Total Proteins Identified", fontsize=20)

    plt.savefig("Total_proteins.svg")
    plt.savefig("Total_proteins.png")


def plot_protein_overlap(protein_dict, labels):
    fig, ax = plt.subplots(figsize=(10, 10))

    if len(protein_dict) in range(2, 6):
        venn(protein_dict, cmap="viridis", ax=ax)
    elif len(protein_dict) == 6:
        pseudovenn(protein_dict, cmap="viridis", ax=ax)
    else:
        print("No Protein Venn Diagram plotted due to invalid number of samples.")
        print("Venn Diagrams require between 2 and 6 samples.")

    ax.legend(labels=labels, fontsize=15, loc="best", bbox_to_anchor=[1.1, 1])
    ax.set_title("Protein Overlap", fontsize=20)
    plt.savefig("Protein_overlap_venn.svg")
    plt.savefig("Protein_overlap_venn.png", bbox_inches="tight")


def plot_protein_upset(protein_dict):
    color = '#21918cff'
    plot_df = upsetplot.from_contents(protein_dict)
    upsetplot.plot(plot_df, sort_by='cardinality', subset_size='auto', facecolor=color)
    # plt.ylim(0, 60)
    plt.title("Distribution of Protein Overlap")

    plt.savefig("Protein_upset.svg")
    plt.savefig("Protein_upset.png")


def plot_overlay_dist(file_list, labels):
    fig, ax = plt.subplots(figsize=(20, 6))
    colors = [plt.cm.viridis(i/float(len(labels)-1)) for i in range(len(labels))]

    for i in range(len(labels)):
        binned = {}
        observed_mz = []
        df = pd.read_csv(file_list[i])
        mz = df["m/z"].tolist()
        observed_mz.append(mz)

        observed_mz = list(set(mz))

        for mass in observed_mz:
            bin_val = (mass//10)*10
            if bin_val not in binned:
                binned[bin_val] = 1
            else:
                binned[bin_val] += 1

        binned = dict(sorted(binned.items(), key=lambda x: x[0]))
        xs = [x for x in binned.keys()]
        ys = [y for y in binned.values()]
        ax.scatter(xs, ys, color=colors[i])
        ax.plot(xs, ys, linestyle="--", color=colors[i], label=labels[i])
        ax.legend(fontsize=15)
        ax.xaxis.set_tick_params(labelsize=13)
        ax.yaxis.set_tick_params(labelsize=13)
        ax.set_title("Overlaid Distribution of detected m/z", fontsize=20)

    plt.savefig("Line_overlay_dist.svg")
    plt.savefig("Line_overlay_dist.png")


def plot_mult_dist(file_list, labels):

    fig, axs = plt.subplots(len(labels), 1, sharex=True, figsize=(20, 20))
    colors = [plt.cm.viridis(i/float(len(labels)-1)) for i in range(len(labels))]

    y_max = 0

    for i in range(len(labels)):
        binned = {}
        observed_mz = []
        df = pd.read_csv(file_list[i])
        mz = df["m/z"].tolist()
        observed_mz.append(mz)

        observed_mz = list(set(mz))

        for mass in observed_mz:
            bin_val = (mass//10)*10
            if bin_val not in binned:
                binned[bin_val] = 1
            else:
                binned[bin_val] += 1

        binned = dict(sorted(binned.items(), key=lambda x: x[0]))
        if max(binned.values()) > y_max:
            y_max = max(binned.values())//10*10+10

        axs[i].bar(binned.keys(), binned.values(), width=8, color=colors[i])
        axs[i].set_title("Distribution of m/z for sample " + labels[i], fontsize=20)
        axs[i].set_ylim(0, y_max)
        axs[i].set_ylabel("Peptide Count", fontsize=12)
        axs[i].yaxis.set_major_locator(plt.MaxNLocator(5))
        axs[i].yaxis.set_tick_params(labelsize=14)

    ax = plt.gca()
    ax.xaxis.set_tick_params(labelsize=14)
    fig.text(0.5, 0.08, 'm/z', ha='center', size=20)

    plt.savefig("Multi_dist.svg")
    plt.savefig("Multi_dist.png")


def plot_retention(file_list):
    hists = []
    for i in range(len(file_list)):
        file = file_list[i]
        df = pd.read_csv(file)
        rt = df.RT.tolist()

        buckets = {}
        if max(rt) < 20:
            window = 2  # minutes
        window = 5  # minutes
        for item in rt:
            time = (item//window) * window
            if time not in buckets:
                buckets[time] = 1
            else:
                buckets[time] += 1

        buckets = dict(sorted(buckets.items(), key=lambda x: x[0], reverse=False))
        hists.append(buckets)

    grad_info = str(input("Would you like to enter gradient information? [y/n]\t")).lower()
    if grad_info == "y":
        grad_x = []
        grad_y = [None]
        while len(grad_x) != len(grad_y):
            grad_x = str(input("Enter time points for gradient (ex. 0, 5, 10...):\t")).split(",")
            grad_x = [float(x) for x in grad_x]
            grad_y = str(input("Enter %B values for gradient (ex. 0, 5, 10...):\t")).split(",")
            grad_y = [float(y) for y in grad_y]
            if len(grad_x) != len(grad_y):
                print("Error: please make the number of points for time and %B are equal.")
    if len(hists) <= 2:
        n_cols = 2
        n_rows = 1
    else:
        n_cols = 3
        n_rows = len(hists) // n_cols
        if len(hists) > n_cols * n_rows:
            n_rows += 1
    fig, ax = plt.subplots(n_rows, n_cols, figsize=(15, 10), sharey=True, sharex=True)
    colors = [plt.cm.viridis(i/float(25-1), alpha=0.75) for i in range(25)]

    if n_rows == 1:
        i = 0
        for c in range(n_cols):
            x = hists[i].keys()
            y = hists[i].values()
            ax1 = ax[c]
            ax1.bar(x, y, width=window, facecolor=colors[8], edgecolor=colors[8])
            ax1.set_title("Dist of RTs, "+labels[i])
            if grad_info == 'y':
                ax2 = ax1.twinx()
                ax2.plot(grad_x, grad_y, color=colors[6])
                ax2.grid(visible=False)
            i += 1

    else:
        i = 0
        for r in range(n_rows):
            for c in range(n_cols):
                if i > len(hists)-1:
                    break
                x = hists[i].keys()
                y = hists[i].values()
                ax1 = ax[r][c]
                ax1.bar(x, y, width=5,  facecolor=colors[8], edgecolor=colors[8])
                ax1.set_title("Dist of RTs, "+labels[i])
                if grad_info == 'y':
                    ax2 = ax1.twinx()
                    ax2.plot(grad_x, grad_y, color=colors[6])
                    ax2.grid(visible=False)
                i += 1
    plt.savefig("RetentionOverlay.svg")
    plt.savefig("RetentionOverlay.png")


def output_file():
    with open("PEAKS_compare_output.txt", "w") as f:
        f.write("Total peptides identified:" + "\n")
        for sample, num in zip(labels, total_peptides):
            f.write(str(sample)+": "+str(num)+"\n")
        f.write("\n")
        f.write("Unique peptides identified:" + "\n")
        for sample, num in zip(labels, unique_peptides):
            f.write(str(sample)+": "+str(num)+"\n")
        f.write("\n")
        f.write("Total proteins identified:" + "\n")
        for sample, num in zip(labels, total_proteins):
            f.write(str(sample)+": "+str(num)+"\n")
        f.write("\n")
        f.write("Unique proteins identified:" + "\n")
        for sample, num in zip(labels, unique_proteins):
            f.write(str(sample)+": "+str(num)+"\n")

        f.close()


def move_imgs():
    files = os.listdir()
    os.mkdir("PNG")
    os.mkdir("SVG")
    pngs = []
    svgs = []
    for file in files:
        if file.endswith(".png"):
            pngs.append(file)
        if file.endswith(".svg"):
            svgs.append(file)

    for item in pngs:
        shutil.move(item, "PNG")
    for item in svgs:
        shutil.move(item, "SVG")


def main():
    get_name_files()
    get_peptides(all_files, labels)
    # print(peptide_dict)
    get_proteins(all_files, labels)

    plot_peptides(total_peptides, unique_peptides, labels)
    plot_peptide_overlap(peptide_dict, labels)
    # plot_peptide_upset(peptide_dict)

    plot_proteins(unique_proteins, labels)
    plot_protein_overlap(protein_dict, labels)
    # plot_protein_upset(protein_dict)

    plot_mult_dist(all_files, labels)
    plot_overlay_dist(all_files, labels)

    plot_retention(all_files)

    output_file()

    move_imgs()


if __name__ == "__main__":
    main()
