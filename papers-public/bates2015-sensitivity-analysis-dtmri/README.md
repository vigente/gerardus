This is the code used for the journal paper:

Title: Sensitivity Analysis of Diffusion Tensor MRI in Simulated Rat Myocardium.
Authors: Joanne Bates, Irvin Teh, Peter Kohl, Jurgen E Schneider and Vicente Grau
Journal: Lecture Notes in Computer Science 2015

Presented at the FIMH conference 2015 in Maastrict, the Netherlands.

How to use the code:

Part 1: Image analysis

1. mean_FA_ADC_images.m used to calculate the FA and ADC for the experimental datasets


Part 2: Sensitivity analysis 

1. Use the 3rd party software, 'Screening (Morris, extended by Campolongo et al 2007)(Matlab)' found at https://sites.google.com/site/jrcsimlab/ to generate the parameters for the sensitivity analysis.
2. Use create_smoldyn_input_file.m to generate the input file for the Smoldyn software (http://www.smoldyn.org/index.html) using the parameters generated in 1.
3. Run the Smoldyn simulations.
4. Convert the Smoldyn output files (.txt files) to Matlab format using import_convert_smoldyn_output_to_matlab.m.
5. Use DT_from_smoldyn_data.m to calculate the ADC and FA for each simulation.
6. Use the 3rd party software, 'Screening (Morris, extended by Campolongo et al 2007)(Matlab)' found at https://sites.google.com/site/jrcsimlab/ to complete the sensitivity analysis.


