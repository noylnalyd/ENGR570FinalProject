# Sensitivity Analysis for an FDM Model of Hypothermia During Shock

  This repository acts as a framework for a numerical human thermoregulation model of induced hypothermia, with sensitivity analysis of inputs performed via nested Monte-Carlo. The code uses a finite-difference formulation of the Fiala et. al. (1999) model to emulate temperatures over time. Analysis of inputs like ambulance times and saline temeprature provide a distribution of procedure completion time, a key factor of patient survival. The code may be modified to investigate the performance of a cooling device across a wide array of environments.

## Description

  The medical procedure modeled in this repository is called Emergency Preservation and Resuscitation (EPR). EPR aims to help patients of extreme hypovolemic shock resulting from external blood loss. Common instances include multiple gunshot injuries or traumatic limb amputation. The motivation of EPR is to mitigate or nullify damage to the brain caused by deoxygenation, extending the survival window before a blood transfusion in primary care. The device modeled aims to induce hypothermia in the brain's interior until metabolism is curtailed, nullifying the patient's need for oxygen temporarily. To achieve this effect, a goal of achieving 10 degrees celcius in the Tympanic membrane (Tty) within 30 minutes is considered sufficient. This tight timeframe means that external factors like external heat loss and the arrival time of first responders play a critical role in determining success.
  This repository models a patient's temperatures at 310 nodes in the body based on the model of Fiala et. al. The effects of hypovolemic shock and later induced hypothermia are modeled from Lyon. The device modeled is that of Konstas et. al. which uses a near-freezing isotonic saline injection delivered via an ECMO device into the intracarotid bifurcation to locally cool the brain tissue. Improvements to the model's skin and active controls are used from Westin, along with arterial blood formulation.
  Four parameters are encoded as inputs to the model and analyzed for their effects on reaching the desired tympanic temperature (Tty) in sufficient time. These are the mean surroundings temperature (TsrmOutdoors), mean ambulance / first responder surroundings temperature (TsrmIndoors), time until cooling begins (trecovery), and temperature of the introduced ECMO saline/blood mix (TECMO). The sensitivity of the time to Tty=10C to each of these inputs is measured and analyzed within the context of EPR. TsrmOutdoors is beyond the control of first responders, but is critical to determining the efficacy of EPR in different climates.

## Getting Started

### Dependencies

* Describe any prerequisites, libraries, OS version, etc., needed before installing program.
* ex. Windows 10

### Installing

* How/where to download your program
* Any modifications needed to be made to files/folders

### Executing program

* How to run the program
* Step-by-step bullets
```
code blocks for commands
```
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

## Help

Any advise for common problems or issues.
```
command to run if program contains helper info
```

## Authors

Contributors names and contact info

ex. Dominique Pizzie  
ex. [@DomPizzie](https://twitter.com/dompizzie)

## Version History

* 0.2
    * Various bug fixes and optimizations
    * See [commit change]() or See [release history]()
* 0.1
    * Initial Release

## License

This project is licensed under the [NAME HERE] License - see the LICENSE.md file for details

## Acknowledgments

Inspiration, code snippets, etc.
* [fiala]()
* [awesome-readme](https://github.com/matiassingers/awesome-readme)
* [PurpleBooth](https://gist.github.com/PurpleBooth/109311bb0361f32d87a2)
* [dbader](https://github.com/dbader/readme-template)
* [zenorocha](https://gist.github.com/zenorocha/4526327)
* [fvcproductions](https://gist.github.com/fvcproductions/1bfc2d4aecb01a834b46)
