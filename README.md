# TERRA Standalone

To this date (2018-07-12), this is a running version of the TERRA standalone. For a detailed documentation, please see the docs folder.

## SNOWPOLINO Branch
This is a branch of TERRA standalone with implementation of the new MCH snow-cover scheme (SNOWPOLINO). The development is ongoing and therefore this seperate branch shall remain as such for atleast a few months or more. 
Principal changes from TERRA standlone:
### From a developer perspective
- significant changes to tsa_input.f90 to including reading in of snowpolino variables from grib files
- added new subroutines for snowpolino: sfc_snow.f90, sfc_snow_init.f90, sfc_snow_data.f90, sfc_snow_utilities.f90
- added all the necessary 'block' variables for snowpolino variables
### From a user perspective
- added a switch **lsnow** -> True activates snowpolino
- added a switch **lsnow_coldstart** -> True if the initial surface data grib file does NOT have snowpolino variables. Set it to False if otherwise.
- ke_snow controls the maximum number of snow layers for snowpolino (along with ekaterina's earlier multi-layer snow model)

## Compilation
Make sure you have the necessary compilers, grib-api and netCDF support installed on your machine. Then use the makefile in the src directory.

cd src
make


## Running
Use the run script runTSA.sh. Make sure you adjust the path specification variables in the beginning
of the script to your file structure. You can alter all namelist switches in this script.

### Poor man's parallelization
Be aware that currently a poor man's parallelization is implemented for the TSA. The domain is divided into subdomains which are then distributed to individual batch jobs (of course you can change that depending on the architecture of your machine). If you do not whish to do these calculations in individual sub-domains, choose nreg_x=nreg_y=1.


## Post processing
As you get one folder per parallelization tile, you will have to patch the output together. This is done with the merge scripts in the tools folder. Again, make sure you adjust the settings to your machine/ file structure. The script merge_last_date.sh can be used in case you just want to merge one output file.


## TO DO
- Consolidate with v. 5.05. or even unify with ICON-TERRA
- Proper parallelization/ port to GPU with openACC
- automatization of running and postprocessing.
