#
# SP: This test is a bit specific and involves a Python script with ase and spglib dependence
#

The dependendencies have been checked with a clean conda environment
conda create --name test python=3.9
conda activate test

pip install ase
pip install spglib

In case you are on a cluster without remote access, you may have to install it differently.
In which case, you also need to comment out the two lines in ../run-ph.sh and look for case number 12.

