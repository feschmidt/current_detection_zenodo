#!/bin/bash

# bash file to be called for executing all figures

#rm -f Figure*.md
#rm -f SMFigure*.md
cd figures/fig1
jupyter nbconvert --to notebook --inplace --execute -y Figure*.ipynb
cd ../fig2
jupyter nbconvert --to notebook --inplace --execute -y Figure*.ipynb
cd ../fig3
jupyter nbconvert --to notebook --inplace --execute -y Figure*.ipynb
cd ../fig4
jupyter nbconvert --to notebook --inplace --execute -y Figure*.ipynb
cd ../fig5
jupyter nbconvert --to notebook --inplace --execute -y Figure*.ipynb
cd ../SMfigs
jupyter nbconvert --to notebook --inplace --execute -y SM*.ipynb
cd ..
#jupytext --to markdown Figure*.ipynb
#jupytext --to markdown SMFigure*.ipynb

#cd figdump/
bash figzip.sh
cd ..
