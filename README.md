# Analog hysteretic relay-based direct control of output voltage in DC-AC converters
---
## Control Engineering Master's Thesis

**Main scripts**
- *Design.m* executes the controller and auxiliar circuitry design process, all steps are broken into individually executable sections inside.

NOTE: The LPRS plot calculation is quite heavy, it is recommended to use the parallel toolbox if available. Alternatively, the precision of said plot can be reduced by increasing *fmstep* in the "Input parameters" section.

- *Analysis.m* contains the needed code to load the variables required by the simulation model and for the subsequent result analysis.

NOTE: The simulation is also quite resource intensive, once computed it can be saved to workspace (there is a line in the script for that) so the results can be reused.

**Auxiliar scripts**
- *DesignFilter.m* and *DesignHelperPlots.m* contain code for generating various figures used in the Design section of the thesis.