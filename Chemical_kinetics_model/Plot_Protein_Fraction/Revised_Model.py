'''
Created on Jul 12, 2023

@author: mina
'''

import matplotlib.pyplot as plt
import matplotlib as mpl
import Settings_Revised_Eq as Settings

mpl.rcParams['pdf.fonttype'] = 42
mpl.rcParams['axes.linewidth'] = Settings.LINE_WIDTH

### This file does all the heavylifting for running the model. You should not need to change it
### to recapitulate my work, but you will need to change it to modify the model.
### I think most of this file is self-explanatory despite the lack of comments, but feel free to
### ask me about any questions you have. 

class Fibril_Model:
    def __init__(self, k_dense_dilute, k_dilute_dense, k_dilute_fibril_primary,
                 k_dilute_fibril_secondary, k_fibril_elongation, k_fibril_loss,
                 initial_dense_protein, initial_dilute_protein, initial_fibril_num,
                 initial_fibril_protein):
        
        self.step_size = Settings.step_size
        self.num_steps = Settings.num_steps
        self.equil_steps = Settings.equil_steps
        
        self.k_dense_dilute = k_dense_dilute
        self.k_dilute_dense = k_dilute_dense
        self.k_dilute_fibril_primary = k_dilute_fibril_primary
        self.k_dilute_fibril_secondary = k_dilute_fibril_secondary
        self.k_fibril_elongation = k_fibril_elongation
        self.k_fibril_loss = k_fibril_loss
        self.initial_dense_protein = initial_dense_protein
        self.initial_dilute_protein = initial_dilute_protein
        self.initial_fibril_num = initial_fibril_num
        self.initial_fibril_protein = initial_fibril_protein

        self.k_dilute_fibril_primary_initial = k_dilute_fibril_primary
        self.k_dilute_fibril_secondary_initial = k_dilute_fibril_secondary
        self.k_fibril_elongation_initial = k_fibril_elongation
        self.k_fibril_loss_initial = k_fibril_loss
        
        self.dilute_protein = self.initial_dilute_protein
        self.dense_protein = self.initial_dense_protein
        self.fibril_num = self.initial_fibril_num
        self.fibril_protein = self.initial_fibril_protein
        
        self.dilute_protein_list = []
        self.dense_protein_list = []
        self.fibril_protein_list = []
        self.time_list = []
        self.protein_lists = [self.dilute_protein_list, self.dense_protein_list, self.fibril_protein_list]
        
        self.data_dic = {'k_de_di': k_dense_dilute, 'k_di_de': k_dilute_dense,
                         'k1': k_dilute_fibril_primary, 'k2': k_dilute_fibril_secondary,
                         'klong': k_fibril_elongation, 'kloss': k_fibril_loss,
                         'conc': initial_dilute_protein + initial_dense_protein, 't05': 0, 't50': 0,
                         't95': 0, 's50': 0, 'tmax': 0, 'smax': 0}

    def dp_dilute_dt(self):
        return self.k_dense_dilute * self.dense_protein - self.k_dilute_dense * self.dilute_protein \
               - self.k_dilute_fibril_primary * self.dilute_protein \
               - self.k_dilute_fibril_secondary * self.dilute_protein * self.fibril_protein \
               - 2 * self.k_fibril_elongation * self.dilute_protein * self.fibril_num \
               + 2 * self.k_fibril_loss * self.fibril_num

    def dp_dense_dt(self):
        return self.k_dilute_dense * self.dilute_protein - self.k_dense_dilute * self.dense_protein \

    def dp_fibril_num_dt(self):
        return self.k_dilute_fibril_primary * self.dilute_protein \
               + self.k_dilute_fibril_secondary * self.dilute_protein * self.fibril_protein

    def dp_fibril_dt(self):
        return self.k_dilute_fibril_primary * self.dilute_protein \
               + self.k_dilute_fibril_secondary * self.dilute_protein * self.fibril_protein \
               + 2 * self.k_fibril_elongation * self.dilute_protein * self.fibril_num \
               - 2 * self.k_fibril_loss * self.fibril_num

    def data_dic_updater(self, i):
            if 0.03 < self.fibril_protein / (self.fibril_protein + self.dilute_protein + self.dense_protein) < 0.05:
                self.data_dic['t05'] = i * self.step_size
            if 0.48 < self.fibril_protein / (self.fibril_protein + self.dilute_protein + self.dense_protein) < 0.50:
                self.data_dic['t50'] = i * self.step_size
            if 0.93 < self.fibril_protein / (self.fibril_protein + self.dilute_protein + self.dense_protein) < 0.95:
                self.data_dic['t95'] = i * self.step_size

            cur_slope = (self.fibril_protein_list[-1] - self.fibril_protein_list[-2]) / self.step_size
            if 0.48 < self.fibril_protein / (self.fibril_protein + self.dilute_protein + self.dense_protein) < 0.50:
                self.data_dic['s50'] = cur_slope
            if cur_slope > self.data_dic['smax']:
                self.data_dic['tmax'] = i * self.step_size
                self.data_dic['smax'] = cur_slope
    
    def step_forward(self):
        new_dilute_protein = self.step_size * self.dp_dilute_dt()
        new_dense_protein = self.step_size * self.dp_dense_dt()
        new_fibril_num = self.step_size * self.dp_fibril_num_dt()
        new_fibril_protein = self.step_size * self.dp_fibril_dt()
        
        self.dilute_protein = self.dilute_protein + new_dilute_protein
        self.dense_protein = self.dense_protein + new_dense_protein
        self.fibril_num = self.fibril_num + new_fibril_num
        self.fibril_protein = self.fibril_protein + new_fibril_protein
    
    def add_to_lists(self, time_step):
        self.dilute_protein_list.append(self.dilute_protein)
        self.dense_protein_list.append(self.dense_protein)
        self.fibril_protein_list.append(self.fibril_protein)
        self.time_list.append(time_step * self.step_size)

    def equilibrate(self):
        self.k_dilute_fibril_primary = 0
        self.k_dilute_fibril_secondary = 0
        self.k_fibril_elongation = 0
        self.k_fibril_loss = 0

        for i in range(self.equil_steps):
            self.step_forward()
            if self.dilute_protein <= 0:
                break

        self.k_dilute_fibril_primary = self.k_dilute_fibril_primary_initial
        self.k_dilute_fibril_secondary = self.k_dilute_fibril_secondary_initial
        self.k_fibril_elongation = self.k_fibril_elongation_initial
        self.k_fibril_loss = self.k_fibril_loss_initial

        print(self.equil_steps)
                
    def run_program(self):
        self.add_to_lists(0)
        for i in range(self.num_steps):
            self.step_forward()
            self.add_to_lists(i + 1)
            self.data_dic_updater(i)

            #if self.dilute_protein <= 0:
            #    break
    
    def tick_param_setter(self, ax, x_direction = 'in', y_direction = 'in'):
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
    
    def ax_filler(self, ax, title, fibril_check=False, fibril_full_check=False):
        ax.set_xlabel('Time', fontfamily = Settings.FONT, size = Settings.AXES_LABEL_SIZE, weight = 'bold')
        ax.set_ylabel('Amount of protein', fontfamily = Settings.FONT, size = Settings.AXES_LABEL_SIZE, weight = 'bold')
        #ax.set_ylim(-0.05, 10.05)
        ax.minorticks_off()
        self.tick_param_setter(ax)
        if fibril_check:
            fibril_d_list = self.protein_lists[-1]
            fibril_i_list = self.protein_lists[-2]
            new_fibril_list = [d + i for d, i in zip(fibril_d_list, fibril_i_list)]
            ax.errorbar(x = self.time_list, y = new_fibril_list,
                        fmt = '-', color=Settings.color_list[-2],
                        linewidth = Settings.LINE_WIDTH,
                        label='Fibrils')
        elif fibril_full_check:
            fibril_d_list = self.protein_lists[-1]
            fibril_i_list = self.protein_lists[-2]
            new_fibril_list = [d + i for d, i in zip(fibril_d_list, fibril_i_list)]
            ax.errorbar(x = self.time_list, y = new_fibril_list,
                        fmt = '-', color='chocolate',
                        linewidth = Settings.LINE_WIDTH,
                        label='Fibrils')
            for protein_index, protein_list in enumerate(self.protein_lists):
                ax.errorbar(x = self.time_list, y = protein_list,
                            fmt = '-', color=Settings.color_list[protein_index],
                            linewidth = Settings.LINE_WIDTH,
                            label=Settings.protein_label_list[protein_index])
        else:
            for protein_index, protein_list in enumerate(self.protein_lists):
                ax.errorbar(x = self.time_list, y = protein_list,
                            fmt = '-', color=Settings.color_list[protein_index],
                            linewidth = Settings.LINE_WIDTH,
                            label=Settings.protein_label_list[protein_index])

                print(Settings.protein_label_list[protein_index])
                file = open(title+'_'+Settings.protein_label_list[protein_index]+'.txt','w')
                for item in range(0,len(self.time_list)):
                    file.write(str(self.time_list[item])+"\t")
                    file.write(str(protein_list[item])+"\n")
                #hi
        
        lgd = ax.legend(loc = 'upper right', edgecolor = 'black',
                        prop={'size': Settings.LEGEND_SIZE, 'family': Settings.FONT, 'weight': 'bold'})
        lgd.get_frame().set_linewidth(Settings.LINE_WIDTH)

    def plotter(self, title, fibril_check = False, fibril_full_check = False):
        fig = plt.figure(figsize=(12, 9))
        rows = 1
        columns = 1
        gs = fig.add_gridspec(rows, columns)
        self.ax_matrix = [None] * 1
        self.ax_matrix[0] = fig.add_subplot(gs[0, 0])
        self.ax_filler(self.ax_matrix[0], title, fibril_check, fibril_full_check)

        plt.tight_layout()
    
    def annotater(self):
        annotation_dic = {'dil_den': self.k_dilute_dense, 'den_dil': self.k_dense_dilute,
                          'fib_1': self.k_dilute_fibril_primary, 'fib_2': self.k_dilute_fibril_secondary,
                          'fib_long': self.k_fibril_elongation, 'fib_loss': self.k_fibril_loss}
        loc = [0.7, 0.67]
        inc = 0.06
        counter = 0
        for txt, val in annotation_dic.items():
            self.ax_matrix[0].annotate(txt + ' = %s' % float('%.1g' % val),
                                       [loc[0], loc[1] - counter * inc], 
                                       xycoords='axes fraction', size = Settings.LEGEND_SIZE,
                                       fontfamily = Settings.FONT, weight = 'bold')
            counter += 1

    def show_plot(self):
        plt.show()
    
    def save_plot(self, cur_path, title, filetype):
        plt.savefig("".join([cur_path, title, filetype]), bbox_inches='tight')


