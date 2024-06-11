'''
Created on Jul 12, 2023

@author: mina
'''

### This file sets some parameters for the model. Make sure to un-comment the correct
### time parameters at the bottom of this file depending on whether you are using the model
### to generate phenomenological plots (use "For running") or fitting the model to
### experimental data (use "For fitting").

TITLE_SIZE = 18
AXES_LABEL_SIZE = 28
AXES_TICK_SIZE = 28
ANNOTATION_SIZE = 12
POINT_LABEL_SIZE = 12
LEGEND_SIZE = 24
MARKER_SIZE = 10
MARKER_EDGE_WIDTH = 0.5
E_LINE_WIDTH = 3
CAP_SIZE = 5
TICK_LENGTH = 10
TICK_WIDTH = 4
LINE_WIDTH = 2
TICK_PAD = 2
X_AXIS_PAD = 2
Y_AXIS_PAD = 2
FONT = 'Arial'

color_list = ['royalblue', 'crimson', 'orange']
protein_label_list = ['Dilute phase', 'Dense phase', 'Fibrils']

#For fitting
step_size = .1
num_steps = 12000
equil_steps = 10000

#===============================================================================
# #For running
# step_size = .001
# num_steps = 50000
# equil_steps = 20000
#===============================================================================
