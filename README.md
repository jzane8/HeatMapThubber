# HeatMapThubber

## maxheat.ulp

NOTE: File locations in the .ulp file are globally hardcoded because ULPs don't actually run out of their file location.

To run, execute maxheat.ulp from the board menu of the desired PCB. It will then access each part of the PCB, its "maxheat"
attribute (manually inputted from the board menu), and its area, and output them to the output.txt file in files.
It then calls drawHeat.py, which uses the values in output.txt to create a simulation.


## drawHeat.py

This script finds the strength and size of the heat sources outlined in output.txt and outputs their heatmap.
The simulation was dropped when I realized how complicated programming multi-material heat simulation was, since
the model I was using was set-state. We also realized that FUSION had the functionality we were looking for anyway. 
