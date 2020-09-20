# ***Pilosocereus aurisetus* ABC scripts and data**
Piloso_simulate_ms_ABC.py - python script to simulate Summary Statistics for the *P. aurisetus* dataset.

msSS.pl - perl script from Naoki Takebayashi to calculate Summary Statistics directly from ms' output.

models.txt - file generated by Piloso_simulate_ms_ABC.py containing the model used in each simulation. 
Have the same order of the file containing the generated Summary Statistics.

SuSt.txt.zip - zipped file containing the Summary Statistics generated by Piloso_simulate_ms_ABC.py.

Emp.txt - Summary statistics calculated from the empirical dataset.

RunABC.R - R script with commands and outputs to load the simulated and empirical data, perform
cross-validation with pseudo-observed data and predict the most likely model with the empirical data.