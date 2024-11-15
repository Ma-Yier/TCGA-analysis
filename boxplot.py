import pickle
import matplotlib.pyplot as plt
import pandas as pd
from search import dataSet
import numpy as np
import math

#plt.ioff()
color = ("lightsalmon", "salmon", "red", "darksalmon", "lightcoral", "indianred",
         "crimson", "firebrick", "darkred", "coral", "orangered", "tomato", "gold", 
        "orange", "darkorange", "lightyellow", "lemonchiffon", "lightgoldenrodyellow", 
         "papayawhip", "moccasin", "peachpuff", "palegoldenrod", "khaki", "darkkhaki", 
         "yellow", "greenyellow", "yellowgreen", "lawngreen", "chartreuse", "lightgreen",
         "forestgreen", "green", "mediumspringgreen", "palegreen" 
         )
#color = ("tomato", "salmon")

with open("result.pkl", "rb") as f:
    result = pickle.load(f)
    ggt_family = ["GGT1", "GGT2", "GGT5", "GGT6", "GGT7"]

    for i in range(0, 5):
        max_value = 0
        # construct new data
        GGTdata = []
        GGTlabel = []
        for item in result[i]:
            if result[i][item].tpm:
                max_value = max(max_value, max(result[i][item].tpm))
                GGTdata.append(result[i][item].tpm)
                GGTlabel.append(item)
        
        # draw boxplot
        fig, ax = plt.subplots()
        #fig.set_size_inches([12.8, 9.6])
        plot = ax.boxplot(
                GGTdata,
                #notch=True, 
                positions=[j for j in range(1,35)], 
                widths=0.5, 
                patch_artist=True,
                #tick_labels=GGTlabel,
                manage_ticks=True,
                showmeans=False, 
                showfliers=True,
                medianprops={"color": "black", "linewidth": 1.5},
                boxprops={"facecolor": "coral", "edgecolor": "black", "linewidth": 1.5},
                whiskerprops={"color": "black", "linewidth": 1.5},
                capprops={"color": "black", "linewidth": 1.5}
                )
        ax.set(xlim=(0, 35), xticks=np.arange(1, 35),
            ylim=(0, math.ceil(max_value)), yticks=np.arange(1, math.ceil(max_value)),
            title=ggt_family[i],
            )
        for patch, color in zip(plot['boxes'], color):
            patch.set_facecolor(color)
        #plt.grid(linestyle="--", alpha=0.3)
        #ax.set_position([0.040, 0.114, 0.940, 0.800])
        #fig.set_dpi(fig.get_dpi()*0.5)
        fig.subplots_adjust(left=0.04, right=0.98, bottom=0.15)
        plt.xticks(np.arange(1, 35), GGTlabel, rotation=-60) #rotate labels
        #plt.show()
        fig.savefig(ggt_family[i])
    
    # draw boxplot
    '''
    plot_notch = ax.boxplot(
        x = GGT1data,         # all data
        notch=True,         # notch or not, default False
        vert=True,          # vertical layout， default True
        widths=0.3,         # box width
        tick_labels=GGT1label,      # box label for each x 
        patch_artist=True,  # fill color for box，default False
        medianprops={       # set middle line attribute
            'linestyle': '--', 'color': 'r', 'linewidth': 1.8
            },
        # showmeans=True,     # show mean point，default False
        # meanline=True,      # show mean line，default False
        # meanprops={         # setter for mean point attr
        #     'marker': 'o', 'markersize': 7.5, 'markeredgewidth': 0.75, 'markerfacecolor': '#b7e1a1', 'markeredgecolor': 'r', 'color': 'k', 'linewidth': 1.5
        #   },
        showflyers=True,    # show flyer，default True
        flyerprops={        # setter flyer attr
            'marker': '^', 'markersize': 6.75, 'markeredgewidth': 0.75, 'markerfacecolor': '#ee5500', 'markeredgecolor': 'k'
        },
        whiskerprops={      # setter vertical line attr
            'linestyle': '--', 'linewidth': 1.2, 'color': '#480656'
        },
        capprops={
            'linestyle': '-', 'linewidth': 1.5, 'color': '#480656'
        },
    )
    '''
    #title_notch = ax.set_title("GGT1")

    # fill color
    #colors = ['pink', 'lightblue', 'lightgreen']
    #for patch, color in zip(plot_notch['boxes'], colors):
    #    patch.set_facecolor(color)

    # set horizontal line and set labels
    #ax.yaxis.grid(True)
    #ax.set_xlabel("disease")
    #ax.set_ylabel("tpm")
    #ax.set(xlim=(0, 35), xticks=np.arange(1, 35),
    #       ylim=(0, math.ceil(max_value)), yticks=np.arange(1, math.ceil(max_value)),
    #       title="GGT1",
    #       )
    #plt.grid(linestyle="--", alpha=0.3)
    #plt.show()
    #fig.savefig('GGT1')