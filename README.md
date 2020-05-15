# Data and code for "Current detection using a Josephson parametric upconverter"

- [Data and code for &quot;Current detection using a Josephson parametric upconverter&quot;](#data-and-code-for-quotcurrent-detection-using-a-josephson-parametric-upconverterquot)
  - [General](#general)
  - [Directory Content](#directory-content)
  - [Data processing and figure generation](#data-processing-and-figure-generation)

This directory contains all the code and data necessary to produce the figures of the paper entitled "Current detection using a Josephson parametric upconverter" written by Felix E. Schmidt, Daniel Bothner, Ines C. Rodrigues, Mario F. Gely, Mark D. Jenkins and Gary A. Steele.

## General

Raw measurement data is structured into subdirectories, each with a `.dat` ASCII file of the measured data, a copy of the measurement `.py` script, a `.dat.gitids.txt` file containing the git hashes of the measurement and analysis libraries, and a `.meta.txt` file for loading the data in [spyview](https://github.com/gsteele13/spyview).
Processed data is usually in the form of `.pkl` files, containing either dictionaries, lists or pandas Dataframes.
Simulation data consists of the simulated data in the `.dat` and the `.meta.txt` file for loading in spyview.
All figures are compiled with _matplotlib_ v3.1.1.

## Directory Content

data_final/
>Processed data, to be loaded for making the final figures.

data_plots/
>Processed data.

data_processed/
>Processed data.

data_raw/
>The raw measurement data.

data_sensitivity/
>Processed data.

Duffing/
>Simulated Duffing responses.

figures/
>Code for generating the individual figures for the main text and the Supplemental Material.

input-output formalism/
>_Mathematica_ notebooks and resulting python source code for calculating the sideband amplitudes.

src/
>This is where all of the source code for data processing is located.
>The zip files `stlab-0.2.zip` and `stlabutils-0.1.zip` are the necessary libraries for measurement and datafile handling (see [stlabutils](https://github.com/steelelab-delft/stlabutils) and [stlab](https://github.com/steelelab-delft/stlab) for latest versions).

*.py, *.ipynb
>Python files or jupyter notebooks to generate the information stored in data_processed/, data_plots/, data_sensitivity/ and data_final/.

README.md
>This file.

python_info.txt
>Generated via `python_info.py`. Contains all modules installed.

tree_info.txt
>Generated via `tree > tree_info.txt`. Tree file with overview of entire directory.

## Data processing and figure generation

The python scripts located in the top-level folder are executed in the below order to generate the final figures from the raw measurement data.
For this, simply execute `bash runall.sh`, which performs the required execution and zips all figures.

1. analysis_classes.py
2. calibration_inputpower.py
3. calibration_gainchain.py
4. Investigating_datafit.py
5. processing_fitpars_vs_Is.py
6. processing_fitpars_vs_Is_F30.py
7. processing_data_F33B_peaks.py
8. Investigating_lossrates.py
9. model_inputoutput_vs_Iset.py
10. processing_data_F33C_4000_peaks.py
11. processing_data_F33C_all_peaks.py
12. model_inputoutput_vs_Delta.py
13. processing_data_F33E_all_peaks.py
14. simulation_ResponseSingleTone_Power.py
15. model_inputoutput_vs_Ppump_Duffing.py
16. CPW_calcs.py
17. JJarrayCPW_analytical_1um_01calcs.ipynb
18. JJarrayCPW_analytical_1um_02sims.ipynb
19. JJarrayCPW_analytical_1um_03proc.ipynb
20. JJarrayCPW_analytical_1um_04sens.ipynb
21. simulation_ResponseSingleTone_Power_all.py
22. model_inputoutput_vs_Ppump_Duffing_vectorized.py
23. Plot_data_Fig1.ipynb
24. Plot_data_Fig2.ipynb
25. Plot_data_Fig4.ipynb
26. figures/fig*/Figure*.ipynb
27. figures/SMfigs/SM*.ipynb
