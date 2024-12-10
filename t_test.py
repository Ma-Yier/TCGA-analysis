import pickle
from search import dataSet
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as mg
import scipy.stats as stats
from scipy.stats import shapiro
from scipy.stats import mannwhitneyu


# Q-Q plot
def QQplot(ca_name, GGT_type, type="tpm", savename="img"):
    GGT = GGT_type
    data = []
    for i in range(0, len(ca_name)):
        if type == "tpm":
            data.append(GGT[ca_name[i]].tpm)
        elif type == "fpkm":
            data.append(GGT[ca_name[i]].fpkm)
        elif type == "fpkmuq":
            data.append(GGT[ca_name[i]].fpkmuq)

    # global params for plot
    fig = plt.figure(figsize=(6, 8))
    plt.matplotlib.rcParams["lines.markersize"] = 0.5
    plt.matplotlib.rcParams["lines.marker"] = "."
    plt.matplotlib.rcParams["lines.linewidth"] = 1.0
    gs = mg.GridSpec(7, 6, hspace=0.5, wspace=0.1)
    axes = {}
    for i in range(0, len(ca_name)):
        axes[i] = fig.add_subplot(gs[i]) 
        stats.probplot(data[i], dist="norm", plot=axes[i])
        axes[i].set_title(f"{ca_name[i]}")
        axes[i].set_xlabel("")
        axes[i].set_ylabel("")
        axes[i].set_xticks([])
        axes[i].set_yticks([])
        print(GGT[ca_name[i]].disease)
    fig.suptitle(f"Q-Q plot for {savename}")
    plt.savefig(savename)
    #plt.show()

        
# Shapiro-Wilk test for normality
def shapiro_test(GGT, ca_name, type="tpm"):
    data = []
    p_values = []
    for i in range(0, len(ca_name)):
        if type == "tpm":
            data.append(GGT[ca_name[i]].tpm)
        elif type == "fpkm":
            data.append(GGT[ca_name[i]].fpkm)
        elif type == "fpkmuq":
            data.append(GGT[ca_name[i]].fpkmuq)
    for i in range(0, len(data)):
        _, p_value = shapiro(data[i])
        p_values.append(p_value)
        if p_value > 0.05:
            print(f"{type}-{ca_name[i]} follows normal distribution.") # print violate list 
    return p_values



# Mann-Whitney U test
def mannwhitneyu_test(GGT, ca_name, type="tpm"):
    datax = []
    ca_list = []
    for i in range(0, len(ca_name)):
        if type == "tpm":
            datax.append(GGT[ca_name[i]].tpm)
        elif type == "fpkm":
            datax.append(GGT[ca_name[i]].fpkm)
        elif type == "fpkmuq":
            datax.append(GGT[ca_name[i]].fpkmuq)
    datay = datax.pop()
    data = []
    nlist = []
    for i in range(0, len(ca_name)-1):
        _, p_value = mannwhitneyu(datax[i], datay, alternative='greater')
        if p_value <= 0.05:
            data.append(datax[i])
            nlist.append(ca_name[i])
    #print(nlist)
    g, gn = sort_mannwhitneyu_test(data, nlist)
    return g, gn


# sort
def sort_mannwhitneyu_test(data_list, name_list):
    if len(data_list) != len(name_list):
        raise TypeError("data_list length should same with name_list.")
    if len(data_list) <= 1:
        return data_list, name_list
    elif len(data_list) == 2:
        _, p_value = mannwhitneyu(data_list[0], data_list[1], alternative='greater')
        if p_value <= 0.05:
            return data_list, name_list
        else:
            return [data_list[1], data_list[0]], [name_list[1], name_list[0]]
    else:
        less = []
        less_name = []
        great = []
        great_name = []
        for i in range(1, len(data_list)):
            _, p_value = mannwhitneyu(data_list[0], data_list[i], alternative='greater')
            if p_value <= 0.05:
                less.append(data_list[i])
                less_name.append(name_list[i])
            else:
                great.append(data_list[i])
                great_name.append(name_list[i])
        small, name_small = sort_mannwhitneyu_test(less, less_name)
        big, name_big = sort_mannwhitneyu_test(great, great_name)
        return big + [data_list[0]] +small, name_big + [name_list[0]] + name_small

if __name__ == "__main__":
    
    # open the results file
    f = open("result.pkl", "rb")
    GGT = pickle.load(f)    # GGT[0] to GGT[4]-> GGT1 GGT2 GGT5 GGT6 GGT7
    ca_name = []
    with open("tcga_abbr.txt", "r", encoding="utf-8") as k:
        while True:
            kline = k.readline()
            if not kline:
                break
            pj_name = kline.strip().split("\t")
            ca_name.append(pj_name[0])

    # draw QQ plot for normality
    GGTfamily = [1, 2, 5, 6, 7]
    for i in range(0, len(GGT)):
        QQplot(ca_name, GGT[i], type="tpm", savename=f"GGT{GGTfamily[i]}-tpm")
        QQplot(ca_name, GGT[i], type="fpkm", savename=f"GGT{GGTfamily[i]}-fpkm")
        QQplot(ca_name, GGT[i], type="fpkmuq", savename=f"GGT{GGTfamily[i]}-fpkmuq")


    # Shapiro-Wilk test
    shapiro_test_pvalues = {}
    for i in range(0, len(GGT)):
        print(f"GGT{GGTfamily[i]}:")
        shapiro_test_pvalues[f"GGT{GGTfamily[i]}-tpm"] = shapiro_test(GGT[i], ca_name, type="tpm")
        shapiro_test_pvalues[f"GGT{GGTfamily[i]}-fpkm"] = shapiro_test(GGT[i], ca_name, type="fpkm")
        shapiro_test_pvalues[f"GGT{GGTfamily[i]}-fpkmuq"] = shapiro_test(GGT[i], ca_name, type="fpkmuq")

    with open("shapiro_result.pkl", "wb") as f:
        pickle.dump(shapiro_test_pvalues, f)

    # GGT5:
    # tpm-HER2 follows normal distribution.
    # GGT7:
    # fpkm-CHOL follows normal distribution.
    # fpkmuq-CHOL follows normal distribution.


    # Mann-Whitney U test
    count_type = ["tpm", "fpkm", "fpkmuq"]
    mannwhitneyu_test_order = {}
    for i in range(0, len(GGT)):
        for j in count_type:
            print(f"GGT{GGTfamily[i]}-{j}:")
            _, glist = mannwhitneyu_test(GGT[i], ca_name, type=j)
            mannwhitneyu_test_order[f"GGT{GGTfamily[i]}_{j}"] = glist
            # print and show the ordered list
            for k in glist:
                print(k)
            print("------------------------")

    '''
    # bar plot
    fig, ax = plt.subplots()
    ax.bar(range(1, len(pvalues_01)+1), pvalues_01, color="blue")
    ax.set_title('Mann-Whitney U test')
    ax.set(xlim=(0, len(pvalues_01)+1), 
                xticks=np.arange(1, len(pvalues_01)+1),
                title='Mann-Whitney U test'
                )
    ca_name.pop()
    plt.xticks(np.arange(1, len(pvalues_01)+1), ca_name, rotation=-70)
    plt.show()
    '''