# Learning-Induced-Odor-Modulation-of-Neuronal-Activity-in-Auditory-Cortex

This repository contains all the data published in the manuscript by Omri Gilday and Adi Mizrahi (Journal of Neuroscience). 

# Visualizing the data: 
To visualize all the data related to the paper in Matlab do the foloowing:
1. Add all the files in this folder to your Matlab path.
2. Type in Matlab the following command:  odor_sound_anlss
3. In the GUI, load the file that contains the data odor_sound_data.mat


# The data file 
The file containing the data is called odor_sound_data.mat 
This file is large (418 MB) and can be downloaded from this google drive: 
TKTK

All files can also be received by an email reques to the corresponding author.

# Brief explanation about the .m files in this repository 
The class odor_sound_data_holder hold all behavioral data for a single mouse/session combination, 
it also holds a odor_sound_neuron object for each neuron containing the responsiveness and statistics of that neuron in that session

all variables are documented within the class files:
odor_sound_data_holder.m
odor_sound_neuron.m

objects of the class Task_GUI_Parameters contain technical data from the task GUI saved during the behavioral session


