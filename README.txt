All use this python enviroment
-conda activate thtfit_env

---------------------------------------------------------------------------------------------------------------
Fitting ThT Data
-Folder: ThT_Fitting_Each_Condition_500_max_test_files
-Run: python Revised_Model_Optimizer_{construct}_{concentration}_Day_{day}.py
-Each file need to change the concentration and reference to construct, concentration, and day
-Also make sure using right csat and optimize for a given time cutoff given curve
-Generally uses best guess from previous concentration for initial conditions - best_guess_{construct}_{prevconcentration}_Day_{day}.txt
-If that doesn't lead to a fit a file called initial_guess_{construct}_{currconcentration}_Day_{day}.txt is used which was sometimes just modulated as need to reduce the residuals

Monomer Fraction in Each Phase
-Folder: Plot_Protein_Fraction
-Run: python Revised_Model_Runner_Read_In_Best_Params.py
-Step size and steps listed in Settings_Revised_Eq.py and were step_size = 0.001, num_steps = 500000, equil_steps = 0
-However since using equation to split dilute and dense based on rates equilibriation is not necessary
-Outputs total protein concentration in each phase as a function of the step

Plot titration of kdilden and kdendil
-Folder: Paramter_Titration
-Run: python Revised_Model_Runner_200_mix.py
-Data generated in Data_kdil_kden_titration
-Plotted using plot_kdilden_kdendil_titration.ipynb
