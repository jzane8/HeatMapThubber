# HeatMapThubber

NOTE: File locations in the .ulp file are globally hardcoded because ULPs don't actually run out of their file location.

To run, execute maxheat.ulp from the board menu of the desired PCB. It will then access each part of the PCB, its "maxheat"
attribute (manually inputted from the board menu), and its area, and output them to the output.txt file in files.
It then calls drawHeat.py, which uses the values in output.txt to create a simulation. 
