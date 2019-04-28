This is repository is for the final project for CSE 397 Fall 2018.

Dependence:

    pytorch

    numpy

    sympy

    scipy

    Cython

Installation:

    cd ./pyStructModel/

    python3 setup.py build_ext --inplace

    cp pyStructModel.cpython-37m-darwin.so ../sampleGenerator/

Directory:

+-- cvxRelax

    # Experiments on sum-of-squares techniques

+-- Naive

    # Naively training

+-- sampleGenerator

    # generating training / test samples 

+-- pyStructModel

    # python wrapper for the code for full stuctural model

Try it out:

    cd ./Naive

    python3 training.py

    python3 curvePlot.py
