#!/bin/bash

# executes all jupyter notebooks to perform full data analysis from raw data to final figures

echo
echo "Warning: This will take quite some time!"
read -p "Are you sure you want to continue? [Y/N] " answer
echo
if [ $answer == "N" ]
then
    echo "Aborting"
elif [ $answer != "Y" ]
then
    echo "Invalid answer. Aborting"
else
    echo "Running data analysis"
    echo
    read -p "What should be the timeout in second? " maxtime 
    echo "Waiting for 10s. Last chance to abort..."
    sleep 10
    
    echo "Executing first set of py files"
    python analysis_classes.py
    python calibration_inputpower.py
    python calibration_gainchain.py
    python Investigating_datafit.py
    python processing_fitpars_vs_Is.py
    python processing_fitpars_vs_Is_F30.py
    python processing_data_F33B_peaks.py
    python Investigating_lossrates.py
    python model_inputoutput_vs_Iset.py
    python processing_data_F33C_4000_peaks.py
    python processing_data_F33C_all_peaks.py
    python model_inputoutput_vs_Delta.py
    python processing_data_F33E_all_peaks.py
    python simulation_ResponseSingleTone_Power.py
    python model_inputoutput_vs_Ppump_Duffing.py
    python CPW_calcs.py
    
    echo "Executing JJarrayCPW ipynbs"
    jupyter nbconvert --to notebook --inplace --execute --ExecutePreprocessor.timeout=$maxtime -y JJarrayCPW_analytical_1um_01calcs.ipynb
    jupyter nbconvert --to notebook --inplace --execute --ExecutePreprocessor.timeout=$maxtime -y JJarrayCPW_analytical_1um_02sims.ipynb
    jupyter nbconvert --to notebook --inplace --execute --ExecutePreprocessor.timeout=$maxtime -y JJarrayCPW_analytical_1um_03proc.ipynb
    jupyter nbconvert --to notebook --inplace --execute --ExecutePreprocessor.timeout=$maxtime -y JJarrayCPW_analytical_1um_04sens.ipynb
    
    echo "Executing last two python files"
    python simulation_ResponseSingleTone_Power_all.py
    python model_inputoutput_vs_Ppump_Duffing_vectorized.py
    
    
    echo "Generating data for figures"
    jupyter nbconvert --to notebook --inplace --execute --ExecutePreprocessor.timeout=$maxtime -y Plot_data_Fig1.ipynb
    jupyter nbconvert --to notebook --inplace --execute --ExecutePreprocessor.timeout=$maxtime -y Plot_data_Fig2.ipynb
    jupyter nbconvert --to notebook --inplace --execute --ExecutePreprocessor.timeout=$maxtime -y Plot_data_Fig4.ipynb
    
    echo "Generating figures"
    bash compilefigs.sh
    
    echo "Documentation"
    python python_info.py
    bash tree_info.sh

fi
