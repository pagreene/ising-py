# ising-py

Overview:
--------
This is an implementation of the metropolis algorithm in c with a python interface. It has the speed of c but the ease of use of python. It is also designed to run many simulations in parallel on the head node. As of writing this, it does not support running on a cluster, and I do not intend to make those changes myself, however others may.

Install:
-------
To install, simply download the files. If you need to modify the c-code, you can download the c file and run make. Otherwise you can just download the libIsing.so file.

Other requirements include:
- python 2.7 (not tested with any other version)
- numpy
- matplotlib

Only standard libraries are used in the c code.

Example:
-------
An example usage may be found in the file testCode.py
