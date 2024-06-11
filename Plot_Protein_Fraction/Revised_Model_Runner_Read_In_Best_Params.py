'''
Created on Jul 19, 2023

@author: mina
'''

import Revised_Model
import pandas as pd

### This file is used for generating fibril plots using the revised model without
### any fitting. This was used to generate the phenomenological plots at the start
### of the slides as well as the parameter correlations in the middle of the slides.
### To generate a plot, un-comment one of the parameter sets below as well as the 
### code block that starts with "fibril_check." "fibril_check" and "fibril_full_check"
### are old flags that tell the program to plot certain curves. They do not work anymore,
### so just set them to 0. If you want to show the plot on the screen, set plot_check = 1.
### If you want to save the plot to the disk, set it to 0. You can then run the program.
### The commented code block at the very bottom of this file was used to generate the
### parameter correlations. You should be able to uncomment it and run it as is (after
### modifying file paths).
###
### Important note: Make sure to use the "For running" parameters in Settings_Revised_Eq.py

#===============================================================================
# # WT two-phase regime
#

currday='1'
currconstruct='WT'
currconc=128

with open('../Code_Kiersten_Each_Condition_500_max_text_files/best_guess_'+currconstruct+'_'+str(currconc)+'_Day_'+currday+'.txt') as myfile: # previous concentrations best guess, Change this
    for line in myfile:
        line = line.strip()
        line = line.split(",")
       
title = currconstruct+'_'+str(currconc)+'_Day_'+currday+'_3phase_protein_amount_equil_0' 
#step_size = .0001
#num_steps = 100000
#equil_steps = 100000
       
k_dense_dilute = float(line[0])
k_dilute_dense = float(line[1])
k_dilute_fibril_primary = float(line[2])
k_dilute_fibril_secondary = float(line[3])
k_fibril_elongation = float(line[4])
k_fibril_loss = 0

initial_conc = currconc
#initial_dilute_protein = 71.2 # used as a test to show that we can't just set this because this is controlled by kin and kout
#initial_dense_protein = currconc-71.2 # used as a test to show that we can't just set this because this is controlled by kin and kout
initial_dilute_protein = initial_conc * k_dense_dilute / (k_dilute_dense + k_dense_dilute)
initial_dense_protein = initial_conc * k_dilute_dense / (k_dilute_dense + k_dense_dilute)

initial_fibril_num = 0
initial_fibril_protein = 0
#===============================================================================


#===============================================================================
fibril_check = False
fibril_full_check = False
plot_check = 0
cur_path = '/project/fava/work/kiersten.ruff/Collaborations/Mittag/2023/ThT_fitting/For_Kiersten/Code_Kiersten_Plot_Concentrations/Figures/'
filetype = '.pdf'
 
revised_model_class = Revised_Model.Fibril_Model(k_dense_dilute, k_dilute_dense,
                                                  k_dilute_fibril_primary, k_dilute_fibril_secondary,
                                                  k_fibril_elongation, k_fibril_loss,
                                                  initial_dense_protein, initial_dilute_protein,
                                                  initial_fibril_num, initial_fibril_protein)
revised_model_class.equilibrate()
revised_model_class.run_program()
revised_model_class.plotter(title, fibril_check, fibril_full_check)
revised_model_class.annotater()
# 
if plot_check:
    revised_model_class.show_plot()
else:
    revised_model_class.save_plot(cur_path, title, filetype)
#===============================================================================

