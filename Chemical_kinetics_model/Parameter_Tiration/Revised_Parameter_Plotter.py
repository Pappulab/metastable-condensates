'''
Created on Jul 25, 2023

@author: mina
'''

### This file is used to generate the 2D histograms for the parameter correlations. After
### changing paths, you should be able to run this file as is.

import numpy as np
import matplotlib.pyplot as plt
import Settings_Revised_Eq as Settings

def tick_param_setter(ax, x_direction = 'in', y_direction = 'in'):
    ax.tick_params(axis = 'x', labelsize = Settings.AXES_TICK_SIZE, direction = x_direction,
        length = Settings.TICK_LENGTH, width = Settings.TICK_WIDTH, color = 'black',
        labelcolor = 'black', pad = Settings.TICK_PAD)
    ax.tick_params(axis = 'y', labelsize = Settings.AXES_TICK_SIZE, direction = y_direction,
        length = Settings.TICK_LENGTH, width = Settings.TICK_WIDTH, color = 'black',
        labelcolor = 'black', pad = Settings.TICK_PAD)
    ax.xaxis.labelpad = Settings.X_AXIS_PAD
    ax.yaxis.labelpad = Settings.Y_AXIS_PAD
    for tick in ax.xaxis.get_major_ticks():
        tick.label1.set_fontsize(Settings.AXES_TICK_SIZE)
        tick.label1.set_fontweight('bold')
        tick.label1.set_fontname(Settings.FONT)
    for tick in ax.yaxis.get_major_ticks():
        tick.label1.set_fontsize(Settings.AXES_TICK_SIZE)
        tick.label1.set_fontweight('bold')
        tick.label1.set_fontname(Settings.FONT)

def ax_filler(ax, x, y, z):
    ax.set_xlabel('k_dil_den', fontfamily = Settings.FONT, size = Settings.AXES_LABEL_SIZE, weight = 'bold')
    ax.set_ylabel('k_den_dil', fontfamily = Settings.FONT, size = Settings.AXES_LABEL_SIZE, weight = 'bold')
    #ax.set_ylim(-0.05, 10.05)
    ax.minorticks_off()
    ax.set_xticks(range(0, 15, 2), [0.01, 0.02, 0.05, 0.1, 0.2, 0.5, 1, 2, 5, 10, 20, 50, 100, 200, 500][::2])
    ax.set_yticks(range(0, 15, 2), [0.01, 0.02, 0.05, 0.1, 0.2, 0.5, 1, 2, 5, 10, 20, 50, 100, 200, 500][::2])
    tick_param_setter(ax)
    data_array = np.zeros((15, 15))
    counter = 0
    for i in range(15):
        for j in range(15):
            data_array[i, j] = z[counter]
            counter += 1
    data_array[data_array == 0] = np.nan
    im = ax.imshow(data_array, origin='lower', cmap = 'GnBu')#, vmin=0, vmax=50)
    plt.colorbar(im)

def plotter(x, y, z, title):
    fig = plt.figure(figsize=(12, 9))
    rows = 1
    columns = 1
    gs = fig.add_gridspec(rows, columns)
    ax_matrix = [None] * 1
    ax_matrix[0] = fig.add_subplot(gs[0, 0])
    ax_filler(ax_matrix[0], x, y, z)
    plt.tight_layout()
    plt.title(title)

def show_plot():
    plt.show()

def save_plot(cur_path, title, filetype):
    plt.savefig("".join([cur_path, title, filetype]), bbox_inches='tight')

cur_path = '/project/fava/work/kiersten.ruff/Collaborations/Mittag/2023/ThT_fitting/For_Kiersten/Revised_Eq_Data/'
param_file = '/project/fava/work/kiersten.ruff/Collaborations/Mittag/2023/ThT_fitting/For_Kiersten/Revised_Eq_Data/Parameters.csv'
k_den_dil_list = []
k_dil_den_list = []
t05_list = []
t50_list = []
t95_list = []
s50_list = []
tmax_list = []
smax_list = []
with open(param_file, 'r') as fp:
    for index, line in enumerate(fp):
        if index == 0:
            continue
        else:
            new_line = [float(x) for x in line.split(',')]
            k_den_dil_list.append(new_line[1])
            k_dil_den_list.append(new_line[2])
            t05_list.append(new_line[9])
            t50_list.append(new_line[10])
            t95_list.append(new_line[11])
            s50_list.append(new_line[12])
            tmax_list.append(new_line[13])
            smax_list.append(new_line[14])

title_t05 = 'Time to 5 percent fibrils'
title_t50 = 'Time to 50 percent fibrils'
title_t95 = 'Time to 95 percent fibrils'
title_s50 = 'Slope at 50 percent fibrils'
title_tmax = 'Time to max fibrillization slope'
title_smax = 'Max fibrillization slope'

plt.rcParams['pdf.fonttype'] = 42 # Makes text editiable
plt.rcParams['ps.fonttype'] = 42
plt.rcParams['font.sans-serif'] = "Arial"
plt.rcParams['font.family'] = "sans-serif"

for z_list, title in zip([t05_list, t50_list, t95_list, s50_list, tmax_list, smax_list],
                         [title_t05, title_t50, title_t95, title_s50, title_tmax, title_smax]):

    plotter(k_den_dil_list, k_dil_den_list, z_list, title)
    save_plot(cur_path, title, '.pdf')

