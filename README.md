# Pulse_Duration_Calibration
Python scripts and data to investigate the calibration of pulse duration for spheroidal approximation to cone model CMEs in solar wind models. Uses the HUXt solar wind model.

HUXt code is written in Python 3.9.13 and has a range of dependencies, which are listed in the `requirements.txt` and `environment.yml` files. Because of these dependencies, the simplest way to work with this code is to use `conda` to create a virtual environment for `HUXt`. We recommend using and up-to-date version of [miniconda](https://docs.anaconda.com/free/miniconda/index.html). With the anaconda prompt, in the root directory of  this repository, this can be done as:
```
>>conda env create -f environment.yml
>>conda activate huxt
``` 



Then run code/calibrate.py
