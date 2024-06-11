'''
Created on Jul 19, 2023

@author: mina
'''

import Revised_Model
import Settings_Revised_Eq as Settings
import matplotlib.pyplot as plt
from scipy.optimize import minimize

### This file is used to perform the data fits. See more information lower in this file.

class Fibril_Optimization:
    def __init__(self, file_name, free_parameter_checklist, initial_guess,
                 bounds_list, csat, time_cutoff, plot_check, plot_save, plot_title):
        
        self.file_name = file_name
        self.free_parameter_checklist = free_parameter_checklist
        self.csat = csat
        self.plot_check = plot_check
        self.plot_save = plot_save
        self.plot_title = plot_title
        self.free_parameter_list = []
        self.fixed_parameter_list = []
        self.bounds_list = []
        self.data_time_list = []
        
        for index, param in enumerate(initial_guess):
            if self.free_parameter_checklist[index] == 0:
                self.fixed_parameter_list.append(param)
            elif self.free_parameter_checklist[index] == 1:
                self.free_parameter_list.append(param)
                self.bounds_list.append(bounds_list[index])
            else:
                continue

        with open(self.file_name, 'r') as fp:
            for index, line in enumerate(fp):
                new_line = line.split()
                if index == 0:
                    self.conc_list = [float(x) for x in new_line[1:]]
                    self.exp_fibril_list = [[] for x in range(len(self.conc_list))]
                elif index <= time_cutoff:
                    self.data_time_list.append(int(new_line[0]))
                    for fibril_index, fibril_conc in enumerate(new_line[1:]):
                        self.exp_fibril_list[fibril_index].append(int(fibril_conc))
                else:
                    break
    
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

    def residual_calculator(self, exp_fibril_list, comp_fibril_list):
        
        res = 0
        for exp_index, exp_point in enumerate(exp_fibril_list):
            comp_point = comp_fibril_list[0::int(1/Settings.step_size)][exp_index]
            res += (exp_point - comp_point) ** 2
        
        return res
    
    def annotater(self, parameter_list, free_parameter_list):
        annotation_dic = {'den_dil': parameter_list[0], 'dil_den': parameter_list[1], 
                          'fib_1': parameter_list[2], 'fib_2': parameter_list[3],
                          'fib_long': parameter_list[4], 'fib_loss': parameter_list[5],
                          'protein_conc': parameter_list[7],
                          'fluor_scaling': free_parameter_list[-1]}
        loc = [0.6, 0.67]
        inc = 0.06
        counter = 0
        for txt, val in annotation_dic.items():
            self.ax_matrix[0].annotate(txt + ' = %s' % float('%.3g' % val),
                                       [loc[0], loc[1] - counter * inc], 
                                       xycoords='axes fraction', size = Settings.LEGEND_SIZE,
                                       fontfamily = Settings.FONT, weight = 'bold')
            counter += 1

    def ax_filler(self, exp_fibril_list, comp_fibril_list):

        ax = self.ax_matrix[0]
        ax.set_xlabel('Time', fontfamily = Settings.FONT, size = Settings.AXES_LABEL_SIZE, weight = 'bold')
        ax.set_ylabel('Fluorescence intensity', fontfamily = Settings.FONT, size = Settings.AXES_LABEL_SIZE, weight = 'bold')
        #ax.set_ylim(-0.05, 10.05)
        ax.minorticks_off()
        self.tick_param_setter(ax)
        ax.errorbar(x = range(len(exp_fibril_list)), y = exp_fibril_list,
                    fmt = '-', color=Settings.color_list[0],
                    linewidth = Settings.LINE_WIDTH,
                    label='Experimental data')
        ax.errorbar(x = range(len(exp_fibril_list)), y = comp_fibril_list[0:len(exp_fibril_list)],
                    fmt = '-', color=Settings.color_list[1],
                    linewidth = Settings.LINE_WIDTH,
                    label='Model')
        
        lgd = ax.legend(loc = 'upper right', edgecolor = 'black',
                        prop={'size': Settings.LEGEND_SIZE, 'family': Settings.FONT, 'weight': 'bold'})
        lgd.get_frame().set_linewidth(Settings.LINE_WIDTH)

    def save_plot(self, plot_title):
        cur_path = '/project/fava/work/kiersten.ruff/Collaborations/Mittag/2023/ThT_fitting/For_Kiersten/Code_Kiersten_Each_Condition_500_max_text_files/Figures/'
        filetype = '.pdf'
        plt.savefig("".join([cur_path, plot_title, filetype]), bbox_inches='tight')

    def plotter(self, exp_fibril_list, comp_fibril_list):

        fig = plt.figure(figsize=(12, 9))
        rows = 1
        columns = 1
        gs = fig.add_gridspec(rows, columns)
        self.ax_matrix = [None] * 1
        self.ax_matrix[0] = fig.add_subplot(gs[0, 0])
        self.ax_filler(exp_fibril_list, comp_fibril_list)

        #print(len(exp_fibril_list))
        #print(len(comp_fibril_list))

        file = open('exp_fibril_WT_65_Day_3.txt','w') #Change this
        for item in exp_fibril_list:
            file.write(str(item)+"\n")
        file = open('comp_fibril_WT_65_Day_3.txt','w') #Change this
        for item in comp_fibril_list:
            file.write(str(item)+"\n")

        plt.tight_layout()
    
    def fibril_maker(self, free_parameter_list, fixed_parameter_list, csat, plot_check = False):
        
        full_res = 0
        
        dilute_conc_index_list = []
        dilute_conc_val_list = []
        for conc_index, conc_val in enumerate(self.conc_list):
            if conc_val < 64 or conc_val > 66: # Change this
                continue
            else:
                dilute_conc_index_list.append(conc_index)
                dilute_conc_val_list.append(conc_val)
        
        print(dilute_conc_val_list)
        #hi
        for conc_index, conc_val in enumerate(dilute_conc_val_list):
            parameter_list = []
            counter_free = 0
            counter_fixed = 0
            for element in self.free_parameter_checklist[:-1]:
                if element == 0:
                    parameter_list.append(fixed_parameter_list[counter_fixed])
                    counter_fixed += 1
                elif element == 1:
                    parameter_list.append(free_parameter_list[counter_free])
                    counter_free += 1
                elif element == 2:
                    parameter_list.append(dilute_conc_val_list[conc_index])

            revised_model_class = Revised_Model.Fibril_Model(*parameter_list)
            revised_model_class.equilibrate()
            revised_model_class.run_program()
            
            # free_parameter_list[-1] is always a scaling factor to bring us up to fluorescence units

            if self.plot_check:
                self.plotter(self.exp_fibril_list[dilute_conc_index_list[conc_index]], [free_parameter_list[-1] * x for x in revised_model_class.fibril_protein_list][0::int(1/Settings.step_size)])
                self.annotater(parameter_list, free_parameter_list)
                plt.show()

            if self.plot_save:
                if conc_index > 0:
                    plot_title = self.plot_title + '_' + str(conc_index)
                else:
                    plot_title = self.plot_title
                self.plotter(self.exp_fibril_list[dilute_conc_index_list[conc_index]], [free_parameter_list[-1] * x for x in revised_model_class.fibril_protein_list][0::int(1/Settings.step_size)])
                self.annotater(parameter_list, free_parameter_list)
                self.save_plot(plot_title)
            
            full_res += self.residual_calculator(self.exp_fibril_list[dilute_conc_index_list[conc_index]],
                                                 [free_parameter_list[-1] * x for x in revised_model_class.fibril_protein_list])
        
        if self.plot_check or self.plot_save:
            exit()
        print(free_parameter_list)
        print(full_res)
        
        return full_res
        
    def optimizer(self):
    
        self.fit = minimize(self.fibril_maker, self.free_parameter_list,
                       (self.fixed_parameter_list, self.csat, self.plot_check),
                       method = 'nelder-mead', bounds = self.bounds_list,
                       options={'maxfev': 5000})

#===============================================================================
# Order of variables
# k_dense_dilute, k_dilute_dense, k_dilute_fibril_primary,
# k_dilute_fibril_secondary, k_fibril_elongation, k_fibril_loss,
# initial_dense_protein, initial_dilute_protein, initial_fibril_num,
# initial_fibril_protein, fluor_scaling
#===============================================================================

### Data are split up by regime and day. After changing the paths, un-comment a section,
### choose a filename for fitting, and select which parameters should be "free" (1) vs.
### "fixed" (0). The entry corresponding to the dilute phase concentration is set at 2
### because this is directly obtained from the experimental concentration in the dataset.
### The order of variables is listed above. Next, choose a time_cutoff for the experimental
### data (in minutes). Lastly, select whether you want to show a plot of the model and
### experiment using the current parameters (plot_check = 1), save the plot to the disk
### (plot_save = 1), or perform an optimization using the current parameters as the initial
### values (plot_check = 0 and plot_save = 0). You should then be able to run this program.
###
### Important note 1: The way to choose which experimental concentration is used for the fit is
### buried at the top of the fibril_maker function. This is because I kept changing how I wanted
### to select this cutoff. For now, if you want to choose specific concentrations, you can change
### the cutoff values directly in that function. Feel free to update this to work better.
###
### Important note 2: Make sure to use the "For fitting" parameters in Settings_Revised_Eq.py 


#===============================================================================
# #One-phase regime Day 2
# using thtfit_env       
file_name = '/project/fava/work/kiersten.ruff/Collaborations/Mittag/2023/ThT_fitting/For_Kiersten/Revised_Eq_Data/Datasets for fitting/Day 3/WT_data.txt'

# This will change whether it is 1 or 2 phase regime
# The first two values will be 0 if 1 phase, and 1 if 2 phase
free_parameter_checklist = [0, 0, 1, 1, 1, 0, 0, 2, 0, 0, 1]

csat_wt = 71.2
csat_d262n = 121.2
csat_d262v = 115.4
csat = csat_wt

# Plot Title - Change this
plot_title = 'WT_65_Day3_plot_through_501'

# For initial optimization  
#time_cutoff = 501 # use value from plot_exp_ThT - Change this
#plot_check = 0
#plot_save = 0

# Once optimized
time_cutoff = 501 # use value from plot_exp_ThT
plot_check = 0
plot_save = 1
       
#starting guess
if plot_save==0:
    #with open('best_guess_WT_52_Day_3.txt') as myfile: # previous concentrations best guess, Change this
    with open('initial_guess_WT_65_Day_3.txt') as myfile: # previous concentrations best guess didn't work so had to restart from here, Change this
        for line in myfile:
            line = line.strip()
            line = line.split(",")
elif plot_save==1:
    with open('best_guess_WT_65_Day_3.txt') as myfile: # current best fit used to plot, Change this
        for line in myfile:
            line = line.strip()
            line = line.split(",")

initial_guess=[]
for l in range(0,len(line)):
    if l!=7:
        tmp=line[l]
        initial_guess.append(float(tmp))
    else:
        initial_guess.append((line[l]))
#print(initial_guess)
#hi 
  
bounds_list = [(0, None)] * 11
        
cur_opt = Fibril_Optimization(file_name, free_parameter_checklist, initial_guess,
                              bounds_list, csat, time_cutoff, plot_check, plot_save,
                              plot_title)
cur_opt.optimizer()
print(cur_opt.fit)
print(cur_opt.fit.x)
print(cur_opt.fit.fun)
# cur_opt.fit.fun = 6238808.503804573 - update this

best_guess=[]
count=-1
for i in range(0,len(initial_guess)):
    if free_parameter_checklist[i]==1:
        count=count+1
        best_guess.append(cur_opt.fit.x[count])
    else:
        best_guess.append(initial_guess[i])
print(best_guess)


import csv
if plot_save==0:
    with open('best_guess_WT_65_Day_3.txt','w') as myoutfile: # Change this
        for i in range(0,len(best_guess)):
            if i < len(best_guess)-1:
                myoutfile.write(str(best_guess[i])+', ')
            else:
                myoutfile.write(str(best_guess[i])+"\n")
#===============================================================================


