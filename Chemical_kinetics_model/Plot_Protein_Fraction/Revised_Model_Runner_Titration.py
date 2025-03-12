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
data_dic = {'k_de_di': [], 'k_di_de': [],
            'k1': [], 'k2': [],
            'klong': [], 'kloss': [],
            'conc': [], 't05': [], 't50': [],
            't95': [], 's50': [], 'tmax': [], 'smax': []}

counter = 0
for den_dil in [0.01, 0.02, 0.05, 0.1, 0.2, 0.5, 1, 2, 5, 10, 20, 50, 100, 200, 500]:
    for dil_den in [0.01, 0.02, 0.05, 0.1, 0.2, 0.5, 1, 2, 5, 10, 20, 50, 100, 200, 500]:
# 
        print(counter)
        
        k_dense_dilute = den_dil
        k_dilute_dense = dil_den
        k_dilute_fibril_primary = 0.001
        k_dilute_fibril_secondary = 0.1
        k_fibril_elongation = 5
        k_fibril_loss = 0
        
        initial_conc = 10
        initial_dilute_protein = initial_conc * k_dense_dilute / (k_dilute_dense + k_dense_dilute)
        initial_dense_protein = initial_conc * k_dilute_dense / (k_dilute_dense + k_dense_dilute)
        initial_fibril_num = 0
        initial_fibril_protein = 0
#         
        fibril_check = False
        fibril_full_check = False
        plot_check = 0
        cur_path = '/project/fava/work/kiersten.ruff/Collaborations/Mittag/2023/ThT_fitting/For_Kiersten/Code_Kiersten_Plot_Concentrations/'
        filetype = '.pdf'
#         
        revised_model_class = Revised_Model.Fibril_Model(k_dense_dilute, k_dilute_dense,
                                                         k_dilute_fibril_primary,
                                                         k_dilute_fibril_secondary,
                                                         k_fibril_elongation, k_fibril_loss,
                                                         initial_dense_protein,
                                                         initial_dilute_protein,
                                                         initial_fibril_num, initial_fibril_protein)
        
        revised_model_class.equilibrate()
        revised_model_class.run_program()
# 
        for key in data_dic.keys():
            data_dic[key].append(revised_model_class.data_dic[key])
#         
        counter += 1
# 
cur_df = pd.DataFrame(dict([(k,pd.Series(v)) for k,v in data_dic.items()]))
print(cur_df)
cur_df.to_csv('/project/fava/work/kiersten.ruff/Collaborations/Mittag/2023/ThT_fitting/For_Kiersten/Code_Kiersten_Plot_Concentrations/Titration_Parameters.csv')
#with pd.ExcelWriter('/project/fava/work/kiersten.ruff/Collaborations/Mittag/2023/ThT_fitting/#For_Kiersten/Code_Kiersten_Plot_Concentrations/Titration_Parameters.xlsx') as writer:
#    cur_df.to_excel(writer, sheet_name='Reference')
#===============================================================================

