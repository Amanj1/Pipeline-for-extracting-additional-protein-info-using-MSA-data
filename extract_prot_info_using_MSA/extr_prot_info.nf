#!/usr/bin/env nextflow

/*
How to run:
nextflow -C extr_prot_info.config run extr_prot_info.nf -profile amanj
*/

params.input='input_data'


/* input files */
//contig sequences
msa_files = Channel.fromFilePairs("${params.input}/*_ORF1_msa_trimmed.fasta",size:1) 


msa_files.into{
msa_barchart_aa;
msa_lineplot_aa;
msa_barchat_pos_aa;
msa_MMDS_rplot;
msa_cluster
}

process barChart_plot_aa{
  tag {"${sample_id}"}

  publishDir "${params.publish_base_dir}/${sample_id}/plots", mode:'link'

  input:
  set sample_id, msa from msa_barchart_aa
  
  output:
  set sample_id, "${sample_id}_stacked_barChart_aa.png"  into barchart_out
  
  script:
""" 
#!/home/amanj/anaconda3/envs/amanjEnv/bin/python3
import sys
import os
import glob
import re
import string
import csv
import json
import collections
import shutil
import re
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
from matplotlib.colors import ListedColormap
plt_name = "testtest.png"

def split(word): 
    return [char for char in word]

def main():
    lengthOfSeq = 0
    f = open("${msa[0]}", "r")
    seq = f.readlines()
    A = 0
    C = 0
    D = 0
    E = 0
    F = 0
    G = 0
    H = 0
    I = 0
    K = 0
    L = 0
    M = 0
    N = 0
    P = 0
    Q = 0
    R = 0
    S = 0
    T = 0
    V = 0
    W = 0
    Y = 0
    mod_seq = []
    collection = []
    tmp1 = ""
    tmp2 = ""
    start = 0
    for i in range(len(seq)):
        if ">" in seq[i] and start == 0:
            tmp1 = seq[i]
            start = 1
        elif i == len(seq) - 1:
            tmp2 = tmp2 + seq[i]
            collection = [tmp1.replace('\\n',''), tmp2.replace('\\n','')]
            mod_seq.append(collection)
            lengthOfSeq = len(tmp2.replace('\\n',''))         
        elif ">" in seq[i] and start == 1:
            collection = [tmp1.replace('\\n',''), tmp2.replace('\\n','')]
            mod_seq.append(collection)
            tmp1 = seq[i]
            tmp2 = ""
        else:
            tmp2 = tmp2 + seq[i]
            
    aa_seq_list = []
    for mseq in mod_seq:
        aa_list = split(mseq[1])
        for aa in aa_list:
            if aa == 'A':
                A = A + 1
            elif aa == 'C':
                C = C + 1
            elif aa == 'D':
                D = D + 1
            elif aa == 'E':
                E = E + 1
            elif aa == 'F':
                F = F + 1
            elif aa == 'G':
                G = G + 1
            elif aa == 'H':
                H = H + 1
            elif aa == 'I':
                I = I + 1
            elif aa == 'K':
                K = K + 1
            elif aa == 'L':
                L = L + 1
            elif aa == 'M':
                M = M + 1
            elif aa == 'N':
                N = N + 1
            elif aa == 'P':
                P = P + 1
            elif aa == 'Q':
                Q = Q + 1
            elif aa == 'R':
                R = R + 1
            elif aa == 'S':
                S = S + 1
            elif aa == 'T':
                T = T + 1
            elif aa == 'V':
                V = V + 1
            elif aa == 'W':
                W = W + 1
            elif aa == 'Y':
                Y = Y + 1
        aa_seq_list.append([mseq[0],[A,C,D,E,F,G,H,I,K,L,M,N,P,Q,R,S,T,V,W,Y]])
        A = 0
        C = 0
        D = 0
        E = 0
        F = 0
        G = 0
        H = 0
        I = 0
        K = 0
        L = 0
        M = 0
        N = 0
        P = 0
        Q = 0
        R = 0
        S = 0
        T = 0
        V = 0
        W = 0
        Y = 0
    ID = []
    container = []
    A1 = []
    C1 = []
    D1 = []
    E1 = []
    F1 = []
    G1 = []
    H1 = []
    I1 = []
    K1 = []
    L1 = []
    M1 = []
    N1 = []
    P1 = []
    Q1 = []
    R1 = []
    S1 = []
    T1 = []
    V1 = []
    W1 = []
    Y1 = []
    A1 = []
    C1 = []
    D1 = []
    E1 = []
    F1 = []
    G1 = []
    H1 = []
    I1 = []
    K1 = []
    L1 = []
    M1 = []
    N1 = []
    P1 = []
    Q1 = []
    R1 = []
    S1 = []
    T1 = []
    V1 = []
    W1 = []
    Y1 = []
    ID.append('categories')
    A1.append('A')
    C1.append('C')
    D1.append('D')
    E1.append('E')
    F1.append('F')
    G1.append('G')
    H1.append('H')
    I1.append('I')
    K1.append('K')
    L1.append('L')
    M1.append('M')
    N1.append('N')
    P1.append('P')
    Q1.append('Q')
    R1.append('R')
    S1.append('S')
    T1.append('T')
    V1.append('V')
    W1.append('W')
    Y1.append('Y')
    for i in range(len(aa_seq_list)):
        ID.append(aa_seq_list[i][0].replace(">","").replace("${sample_id}_",""))
        index = 0
        A1.append(aa_seq_list[i][1][index])
        index = index + 1
        C1.append(aa_seq_list[i][1][index])
        index = index + 1
        D1.append(aa_seq_list[i][1][index])
        index = index + 1
        E1.append(aa_seq_list[i][1][index])
        index = index + 1
        F1.append(aa_seq_list[i][1][index])
        index = index + 1
        G1.append(aa_seq_list[i][1][index])
        index = index + 1
        H1.append(aa_seq_list[i][1][index])
        index = index + 1
        I1.append(aa_seq_list[i][1][index])
        index = index + 1
        K1.append(aa_seq_list[i][1][index])
        index = index + 1
        L1.append(aa_seq_list[i][1][index])
        index = index + 1
        M1.append(aa_seq_list[i][1][index])
        index = index + 1
        N1.append(aa_seq_list[i][1][index])
        index = index + 1
        P1.append(aa_seq_list[i][1][index])
        index = index + 1
        Q1.append(aa_seq_list[i][1][index])
        index = index + 1
        R1.append(aa_seq_list[i][1][index])
        index = index + 1
        S1.append(aa_seq_list[i][1][index])
        index = index + 1
        T1.append(aa_seq_list[i][1][index])
        index = index + 1
        V1.append(aa_seq_list[i][1][index])
        index = index + 1
        W1.append(aa_seq_list[i][1][index])
        index = index + 1
        Y1.append(aa_seq_list[i][1][index])
    container = [A1,C1,D1,E1,F1,G1,H1,I1,K1,L1,M1,N1,P1,Q1,R1,S1,T1,V1,W1,Y1]
    #print(container)
    #print(ID)
    df = pd.DataFrame(columns=ID, data=container)
    df.set_index('categories')\
    .reindex(df.set_index('categories').sum().sort_values(ascending=False).index, axis=1)\
    .T.plot(kind='bar', stacked=True, colormap=ListedColormap(sns.color_palette("Paired", 20)), figsize=(24,12))
    plt.title("Amino acid composition for ${sample_id} with sequence length of " + str(int(lengthOfSeq)) + " and " + str(len(ID) - 1) + " number of sequences" )
    plt.xlabel("Sequence id")
    plt.ylabel("Number of occurrences")
    plt.xticks(size = 8, rotation=15)
    plt.legend(loc='center left', bbox_to_anchor=(1, 0.5))
    #plt.show()
    plt.savefig("${sample_id}_stacked_barChart_aa.png")
main()
"""
}

process barchart_pos_comp_aa{
  tag {"${sample_id}"}

  publishDir "${params.publish_base_dir}/${sample_id}/plots", mode:'link'

  input:
  set sample_id, msa from msa_barchat_pos_aa
  
  output:
  set sample_id, "${sample_id}_stacked_barchart_aa_for_each_pos_in_msa.png"  into msa_barchat_pos_aa_out
  
  script:
"""
#!/home/amanj/anaconda3/envs/amanjEnv/bin/python3
import sys
import os
import glob
import re
import string
import csv
import json
import collections
import shutil
import re
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
from matplotlib.colors import ListedColormap

plt_name = "${sample_id}_stacked_barchart_aa_for_each_pos_in_msa.png"

f1 = open("${msa[0]}", "r")
msa = f1.readlines()
tmp = ""
all_seq = []
nr_seq = 0
for m in msa:
    if ">" in m:
        nr_seq = nr_seq + 1
        all_seq.append(tmp)
        tmp = ""
    else:
        tmp = tmp + m.strip('\\n')
all_seq.append(tmp)
all_seq.pop(0)
tmp = []
collection = []
graph_data = []
for i in range(len(all_seq[0])):
    for j in range(len(all_seq)):
        tmp.append(all_seq[j][i])
    collection.append(tmp)
    tmp = []
A = 0
C = 0
D = 0
E = 0
F = 0
G = 0
H = 0
I = 0
K = 0
L = 0
M = 0
N = 0
P = 0
Q = 0
R = 0
S = 0
T = 0
V = 0
W = 0
Y = 0
ID = []
container = []
A1 = []
C1 = []
D1 = []
E1 = []
F1 = []
G1 = []
H1 = []
I1 = []
K1 = []
L1 = []
M1 = []
N1 = []
P1 = []
Q1 = []
R1 = []
S1 = []
T1 = []
V1 = []
W1 = []
Y1 = []
aa_seq_list = []
count = 1
for c in collection:
    for aa in c:
        if aa == 'A':
            A = A + 1
        elif aa == 'C':
            C = C + 1
        elif aa == 'D':
            D = D + 1
        elif aa == 'E':
            E = E + 1
        elif aa == 'F':
            F = F + 1
        elif aa == 'G':
            G = G + 1
        elif aa == 'H':
            H = H + 1
        elif aa == 'I':
            I = I + 1
        elif aa == 'K':
            K = K + 1
        elif aa == 'L':
            L = L + 1
        elif aa == 'M':
            M = M + 1
        elif aa == 'N':
            N = N + 1
        elif aa == 'P':
            P = P + 1
        elif aa == 'Q':
            Q = Q + 1
        elif aa == 'R':
            R = R + 1
        elif aa == 'S':
            S = S + 1
        elif aa == 'T':
            T = T + 1
        elif aa == 'V':
            V = V + 1
        elif aa == 'W':
            W = W + 1
        elif aa == 'Y':
            Y = Y + 1
    aa_seq_list.append([count,[A,C,D,E,F,G,H,I,K,L,M,N,P,Q,R,S,T,V,W,Y]])
    count = count + 1
    A = 0
    C = 0
    D = 0
    E = 0
    F = 0
    G = 0
    H = 0
    I = 0
    K = 0
    L = 0
    M = 0
    N = 0
    P = 0
    Q = 0
    R = 0
    S = 0
    T = 0
    V = 0
    W = 0
    Y = 0
A1 = []
C1 = []
D1 = []
E1 = []
F1 = []
G1 = []
H1 = []
I1 = []
K1 = []
L1 = []
M1 = []
N1 = []
P1 = []
Q1 = []
R1 = []
S1 = []
T1 = []
V1 = []
W1 = []
Y1 = []
ID.append('categories')
A1.append('A')
C1.append('C')
D1.append('D')
E1.append('E')
F1.append('F')
G1.append('G')
H1.append('H')
I1.append('I')
K1.append('K')
L1.append('L')
M1.append('M')
N1.append('N')
P1.append('P')
Q1.append('Q')
R1.append('R')
S1.append('S')
T1.append('T')
V1.append('V')
W1.append('W')
Y1.append('Y')
for i in range(len(aa_seq_list)):
    ID.append(aa_seq_list[i][0])
    index = 0
    A1.append(aa_seq_list[i][1][index])
    index = index + 1
    C1.append(aa_seq_list[i][1][index])
    index = index + 1
    D1.append(aa_seq_list[i][1][index])
    index = index + 1
    E1.append(aa_seq_list[i][1][index])
    index = index + 1
    F1.append(aa_seq_list[i][1][index])
    index = index + 1
    G1.append(aa_seq_list[i][1][index])
    index = index + 1
    H1.append(aa_seq_list[i][1][index])
    index = index + 1
    I1.append(aa_seq_list[i][1][index])
    index = index + 1
    K1.append(aa_seq_list[i][1][index])
    index = index + 1
    L1.append(aa_seq_list[i][1][index])
    index = index + 1
    M1.append(aa_seq_list[i][1][index])
    index = index + 1
    N1.append(aa_seq_list[i][1][index])
    index = index + 1
    P1.append(aa_seq_list[i][1][index])
    index = index + 1
    Q1.append(aa_seq_list[i][1][index])
    index = index + 1
    R1.append(aa_seq_list[i][1][index])
    index = index + 1
    S1.append(aa_seq_list[i][1][index])
    index = index + 1
    T1.append(aa_seq_list[i][1][index])
    index = index + 1
    V1.append(aa_seq_list[i][1][index])
    index = index + 1
    W1.append(aa_seq_list[i][1][index])
    index = index + 1
    Y1.append(aa_seq_list[i][1][index])
container = [A1,C1,D1,E1,F1,G1,H1,I1,K1,L1,M1,N1,P1,Q1,R1,S1,T1,V1,W1,Y1]
#print(container)
#print(ID)
df = pd.DataFrame(columns=ID, data=container)
df.set_index('categories')\
.reindex(df.set_index('categories').sum().index, axis=1)\
.T.plot(kind='bar', stacked=True, colormap=ListedColormap(sns.color_palette("Paired", 20)), figsize=(24,12))
plt.title("amino acid composition for ${sample_id}" + " with " + str(len(aa_seq_list))+ " positions and " +str(nr_seq)+ " number of sequences")
plt.xlabel("sequence position")
plt.ylabel("number of occurrences")
plt.xticks(size = 8, rotation=15)
plt.legend(loc='center left', bbox_to_anchor=(1, 0.5))
#plt.show()
plt.savefig(plt_name)
f1.close()
"""
}

process lineplot_aa{
  tag {"${sample_id}"}

  publishDir "${params.publish_base_dir}/${sample_id}/plots", mode:'link'

  input:
  set sample_id, msa from msa_lineplot_aa
  
  output:
  set sample_id, "${sample_id}_lineplot_aa_for_each_pos_in_msa.png"  into msa_lineplot_aa_out
  
  script:
"""
#!/home/amanj/anaconda3/envs/amanjEnv/bin/python3
import matplotlib.pyplot as plt
f1 = open("${msa[0]}", "r")
msa = f1.readlines()
tmp = ""
all_seq = []
nr_seq = 0
for m in msa:
    if ">" in m:
        nr_seq = nr_seq + 1
        all_seq.append(tmp)
        tmp = ""
    else:
        tmp = tmp + m.strip('\\n')
all_seq.append(tmp)
all_seq.pop(0)
tmp = []
collection = []
graph_data = []
for i in range(len(all_seq[0])):
    for j in range(len(all_seq)):
        tmp.append(all_seq[j][i])
    collection.append(tmp)
    tmp = []
for c in collection:
    tmp = []
    tmp2 = [] 
    [tmp.append(x) for x in c if x not in tmp]
    if '-' in tmp:
        for t in tmp:
            if "-" not in t:
                tmp2.append(t) 
        graph_data.append(len(tmp2))
    else:
        graph_data.append(len(tmp))


#print(graph_data)
fig, ax = plt.subplots(figsize=(20, 10))
plt.plot(graph_data)
title = "Amino acid variation/position for ${sample_id} with total sequence length of " + str(len(all_seq[0])) + " and " + str(nr_seq) + " number of sequences"
plt.title(title)
plt.ylabel('Number of amino acid variation')
plt.xlabel('Amino acid position in multiple sequence alignment')
#plt.show()
plt.savefig("${sample_id}_lineplot_aa_for_each_pos_in_msa.png")
f1.close()
"""
}

process mmds_plot{
  tag {"${sample_id}"}

  publishDir "${params.publish_base_dir}/${sample_id}/plots", mode:'link'

  input:
  set sample_id, msa from msa_MMDS_rplot
  
  output:
  set sample_id, "${sample_id}_mmds_plot.pdf"   into MMDS_rplot_out
  
  script:
""" 
#!/usr/bin/Rscript

library("rgl")
library("amap")
library("e1071")
library("scales")
library("cluster")
library("rgl")
library(bios2mds)
#library(tidyverse)
### MMDS analysis

fasta <- import.fasta("${msa[0]}", aa.to.upper = TRUE, gap.to.dash = TRUE)
dif_fasta <- mat.dif(fasta, fasta)
mmds <- mmds(dif_fasta)

png(file="${sample_id}_mmds_plot.png", width=600, height=350)
mmds.2D <- mmds.2D.plot(mmds, legend = FALSE)
dev.off()

file.rename("Rplots.pdf", "${sample_id}_mmds_plot.pdf")

"""
}

process MMseqs2_clustering{

  tag {"${sample_id}"}

  publishDir "${params.publish_base_dir}/${sample_id}/cluster_results", mode:'link'

  input:
  set sample_id, msa from msa_cluster
  
  output:
  set sample_id, "${sample_id}_clusterRes_cluster.tsv", "${sample_id}_clusterRes_all_seqs.fasta", "${sample_id}_clusterRes_rep_seq.fasta"  into cluster_res
  
  script:
"""
mmseqs easy-linclust ${msa[0]} clusterRes tmp 
sleep 2
mv clusterRes_cluster.tsv "${sample_id}_clusterRes_cluster.tsv"
mv clusterRes_all_seqs.fasta "${sample_id}_clusterRes_all_seqs.fasta"
mv clusterRes_rep_seq.fasta "${sample_id}_clusterRes_rep_seq.fasta"
"""
}

cluster_res.into{
clstr_barchart;
clstr_stacked_barchart
}

process clstr_barchart_plot{
  tag {"${sample_id}"}

  publishDir "${params.publish_base_dir}/${sample_id}/plots", mode:'link'

  input:
  set sample_id, clstrTSV, clstrFA, clstrOT from clstr_barchart

  output:
  set sample_id, "${sample_id}_cluster_list.txt", "${sample_id}_barChart_of_cluster_data.png"  into cluster_plot_out

  script:
"""
#!/home/amanj/anaconda3/envs/amanjEnv/bin/python3

import numpy as np
import matplotlib
import matplotlib.pyplot as plt

f1 = open("${clstrTSV}", "r")
clstr = f1.readlines()
f2 = open("${sample_id}_cluster_list.txt", "a")

pr_clstr = ""
collection = []
tmp = []
for c in clstr:
    print(c.split()[0], pr_clstr)
    if c.split()[0] != pr_clstr:
        collection.append(tmp)
        tmp = []
        tmp.append(c)
        pr_clstr = c.split()[0]
    else:
        tmp.append(c)
        pr_clstr = c.split()[0]
collection.append(tmp)
collection.pop(0)

count = 1
clstr_size = []
clstr_name = []
for cn in collection:
    clstr_size.append(len(cn))
    clstr_name.append("C" + str(count))
    f2.write("cluster " + str(count) + " with " + str(len(cn)) + " number of sequences:\\n")
    count = count + 1
    for seq in cn:
        f2.write(seq.split('\\t')[1])
    f2.write('\\n')
f1.close()
f2.close()

x = np.arange(len(clstr_size))  # the label locations
width = 0.35  # the width of the bars

fig, ax = plt.subplots(figsize=(20,10))
rects1 = ax.bar(x, clstr_size, width)

# Add some text for labels, title and custom x-axis tick labels, etc.
ax.set_ylabel('Number of sequences')
ax.set_title('Clustering of amino acid sequnces based on msa data\\nIn total ' + str(len(clstr_name)) +' clusters')
ax.set_xticks(x)
ax.set_xticklabels(clstr_name)
ax.legend().remove()

def autolabel(rects):
    #Attach a text label above each bar in *rects*, displaying its height.
    for rect in rects:
        height = rect.get_height()
        ax.annotate('{}'.format(height),
                    xy=(rect.get_x() + rect.get_width() / 2, height),
                    xytext=(0, 3),  # 3 points vertical offset
                    textcoords="offset points",
                    ha='center', va='bottom')
autolabel(rects1)
fig.tight_layout()
plt.savefig("${sample_id}_barChart_of_cluster_data.png")
"""
}

/*
process clstr_stacked_barChart_plot{
  tag {"${sample_id}"}

  publishDir "${params.publish_base_dir}/${sample_id}/plots", mode:'link'

  input:
  set sample_id, clstrTSV, clstrFA, clstrOT from clstr_stacked_barchart 
  
  output:
  set sample_id, "${sample_id}_stacked_barChart_aa.png"  into stacked_barchart_out
  
  script:
""" 
#!/home/amanj/anaconda3/envs/amanjEnv/bin/python3
import sys
import os
import glob
import re
import string
import csv
import json
import collections
import shutil
import re
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
from matplotlib.colors import ListedColormap
plt_name = "testtest.png"

def split(word): 
    return [char for char in word]

def main():
    lengthOfSeq = 0
    f = open("${msa[0]}", "r")
    seq = f.readlines()
    A = 0
    C = 0
    D = 0
    E = 0
    F = 0
    G = 0
    H = 0
    I = 0
    K = 0
    L = 0
    M = 0
    N = 0
    P = 0
    Q = 0
    R = 0
    S = 0
    T = 0
    V = 0
    W = 0
    Y = 0
    mod_seq = []
    collection = []
    tmp1 = ""
    tmp2 = ""
    start = 0
    for i in range(len(seq)):
        if ">" in seq[i] and start == 0:
            tmp1 = seq[i]
            start = 1
        elif i == len(seq) - 1:
            tmp2 = tmp2 + seq[i]
            collection = [tmp1.replace('\\n',''), tmp2.replace('\\n','')]
            mod_seq.append(collection)
            lengthOfSeq = len(tmp2.replace('\\n',''))         
        elif ">" in seq[i] and start == 1:
            collection = [tmp1.replace('\\n',''), tmp2.replace('\\n','')]
            mod_seq.append(collection)
            tmp1 = seq[i]
            tmp2 = ""
        else:
            tmp2 = tmp2 + seq[i]
            
    aa_seq_list = []
    for mseq in mod_seq:
        aa_list = split(mseq[1])
        for aa in aa_list:
            if aa == 'A':
                A = A + 1
            elif aa == 'C':
                C = C + 1
            elif aa == 'D':
                D = D + 1
            elif aa == 'E':
                E = E + 1
            elif aa == 'F':
                F = F + 1
            elif aa == 'G':
                G = G + 1
            elif aa == 'H':
                H = H + 1
            elif aa == 'I':
                I = I + 1
            elif aa == 'K':
                K = K + 1
            elif aa == 'L':
                L = L + 1
            elif aa == 'M':
                M = M + 1
            elif aa == 'N':
                N = N + 1
            elif aa == 'P':
                P = P + 1
            elif aa == 'Q':
                Q = Q + 1
            elif aa == 'R':
                R = R + 1
            elif aa == 'S':
                S = S + 1
            elif aa == 'T':
                T = T + 1
            elif aa == 'V':
                V = V + 1
            elif aa == 'W':
                W = W + 1
            elif aa == 'Y':
                Y = Y + 1
        aa_seq_list.append([mseq[0],[A,C,D,E,F,G,H,I,K,L,M,N,P,Q,R,S,T,V,W,Y]])
        A = 0
        C = 0
        D = 0
        E = 0
        F = 0
        G = 0
        H = 0
        I = 0
        K = 0
        L = 0
        M = 0
        N = 0
        P = 0
        Q = 0
        R = 0
        S = 0
        T = 0
        V = 0
        W = 0
        Y = 0
    ID = []
    container = []
    A1 = []
    C1 = []
    D1 = []
    E1 = []
    F1 = []
    G1 = []
    H1 = []
    I1 = []
    K1 = []
    L1 = []
    M1 = []
    N1 = []
    P1 = []
    Q1 = []
    R1 = []
    S1 = []
    T1 = []
    V1 = []
    W1 = []
    Y1 = []
    A1 = []
    C1 = []
    D1 = []
    E1 = []
    F1 = []
    G1 = []
    H1 = []
    I1 = []
    K1 = []
    L1 = []
    M1 = []
    N1 = []
    P1 = []
    Q1 = []
    R1 = []
    S1 = []
    T1 = []
    V1 = []
    W1 = []
    Y1 = []
    ID.append('categories')
    A1.append('A')
    C1.append('C')
    D1.append('D')
    E1.append('E')
    F1.append('F')
    G1.append('G')
    H1.append('H')
    I1.append('I')
    K1.append('K')
    L1.append('L')
    M1.append('M')
    N1.append('N')
    P1.append('P')
    Q1.append('Q')
    R1.append('R')
    S1.append('S')
    T1.append('T')
    V1.append('V')
    W1.append('W')
    Y1.append('Y')
    for i in range(len(aa_seq_list)):
        ID.append(aa_seq_list[i][0].replace(">","").replace("${sample_id}_",""))
        index = 0
        A1.append(aa_seq_list[i][1][index])
        index = index + 1
        C1.append(aa_seq_list[i][1][index])
        index = index + 1
        D1.append(aa_seq_list[i][1][index])
        index = index + 1
        E1.append(aa_seq_list[i][1][index])
        index = index + 1
        F1.append(aa_seq_list[i][1][index])
        index = index + 1
        G1.append(aa_seq_list[i][1][index])
        index = index + 1
        H1.append(aa_seq_list[i][1][index])
        index = index + 1
        I1.append(aa_seq_list[i][1][index])
        index = index + 1
        K1.append(aa_seq_list[i][1][index])
        index = index + 1
        L1.append(aa_seq_list[i][1][index])
        index = index + 1
        M1.append(aa_seq_list[i][1][index])
        index = index + 1
        N1.append(aa_seq_list[i][1][index])
        index = index + 1
        P1.append(aa_seq_list[i][1][index])
        index = index + 1
        Q1.append(aa_seq_list[i][1][index])
        index = index + 1
        R1.append(aa_seq_list[i][1][index])
        index = index + 1
        S1.append(aa_seq_list[i][1][index])
        index = index + 1
        T1.append(aa_seq_list[i][1][index])
        index = index + 1
        V1.append(aa_seq_list[i][1][index])
        index = index + 1
        W1.append(aa_seq_list[i][1][index])
        index = index + 1
        Y1.append(aa_seq_list[i][1][index])
    container = [A1,C1,D1,E1,F1,G1,H1,I1,K1,L1,M1,N1,P1,Q1,R1,S1,T1,V1,W1,Y1]
    #print(container)
    #print(ID)
    df = pd.DataFrame(columns=ID, data=container)
    df.set_index('categories')\
    .reindex(df.set_index('categories').sum().sort_values(ascending=False).index, axis=1)\
    .T.plot(kind='bar', stacked=True, colormap=ListedColormap(sns.color_palette("Paired", 20)), figsize=(24,12))
    plt.title("Amino acid composition for ${sample_id} with sequence length of " + str(int(lengthOfSeq)) + " and " + str(len(ID) - 1) + " number of sequences" )
    plt.xlabel("Sequence id")
    plt.ylabel("Number of occurrences")
    plt.xticks(size = 8, rotation=15)
    plt.legend(loc='center left', bbox_to_anchor=(1, 0.5))
    plt.savefig("${sample_id}_stacked_barChart_aa.png")
main()
"""
}
*/











