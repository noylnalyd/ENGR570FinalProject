# Sensitivity Analysis for an FDM Model of Hypothermia During Shock

  This repository acts as a framework for a numerical human thermoregulation model of induced hypothermia, with sensitivity analysis of inputs performed via nested Monte-Carlo. The code uses a finite-difference formulation of the Fiala et. al. (1999) model to emulate temperatures over time. Analysis of inputs like ambulance times and saline temeprature provide a distribution of procedure completion time, a key factor of patient survival. The code may be modified to investigate the performance of a cooling device across a wide array of environments.

## Description

  The medical procedure modeled in this repository is called Emergency Preservation and Resuscitation (EPR). EPR aims to help patients of extreme hypovolemic shock resulting from external blood loss. Common instances include multiple gunshot injuries or traumatic limb amputation. The motivation of EPR is to mitigate or nullify damage to the brain caused by deoxygenation, extending the survival window before a blood transfusion in primary care. The device modeled aims to induce hypothermia in the brain's interior until metabolism is curtailed, nullifying the patient's need for oxygen temporarily. To achieve this effect, a goal of achieving 10 degrees celcius in the Tympanic membrane (Tty) within 30 minutes is considered sufficient. This tight timeframe means that external factors like external heat loss and the arrival time of first responders play a critical role in determining success.
  This repository models a patient's temperatures at 310 nodes in the body based on the model of Fiala et. al. The effects of hypovolemic shock and later induced hypothermia are modeled from Lyon. The device modeled is that of Konstas et. al. which uses a near-freezing isotonic saline injection delivered via an ECMO device into the intracarotid bifurcation to locally cool the brain tissue. Improvements to the model's skin and active controls are used from Westin, along with arterial blood formulation.
  Four parameters are encoded as inputs to the model and analyzed for their effects on reaching the desired tympanic temperature (Tty) in sufficient time. These are the mean surroundings temperature (TsrmOutdoors), mean ambulance / first responder surroundings temperature (TsrmIndoors), time until cooling begins (trecovery), and temperature of the introduced ECMO saline/blood mix (TECMO). The sensitivity of the time to Tty=10C to each of these inputs is measured and analyzed within the context of EPR. TsrmOutdoors is beyond the control of first responders, but is critical to determining the efficacy of EPR in different climates.

## Getting Started

## Running Konstas Case example

### Dependencies

* C++ libraries (math.h, cmath, chrono)

### Executing Program
* Navigate to main repository folder
* Open terminal and execute:
```
g++ runKonstas.cpp -o runKonstas.exe
./runKonstas.exe <TsrmIndoors> <TsrmOutdoors> <trecovery> <TECMO>
```
Provide values for temperature in Kelvin and for ```trecovery``` in seconds. Note that the solver may not be stable for extreme values of temperature or trecovery.

## Running Global Sensitivity Analysis Example

### Dependencies

* Python3 and related modules (pandas, scipy, numpy, and matplotlib
* C++ libraries (math.h, cmath, fstream, string, string.h, sstream)

### Executing Program

* Navigate to 'uq' subdirectory
* In terminal:
```
bash parallelRun.sh nProcesses nSamples
```
(nProcesses is the number of parallel processes to run on) 
(nSamples is the size of each Monte Carlo sum; total computational cost is 2 x nSamples x D, where D=4 parameters in this example)

## Authors

Contributors names and contact info

ex. Dylan Lyon [@noylnalyd](https://github.com/noylnalyd)
ex. Thomas Coons[@tcoonsUM](https://github.com/tcoonsUM)


## License

This project is licensed under the MIT License - see the LICENSE.md file for details
