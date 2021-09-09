
   Terra Stand Alone for the blocked version of TERRA
   ==================================================

Introduction:
-------------

Due to the blocked data structure in the COSMO-Model, that we now use for the 
parameterizations, also the TSA (Terra Stand Alone) had to be adapted to use
the new COSMO modules sfc_terra.f90 (sfc_terra_data.f90, sfc_terra_init.f90 and
sfc_utilities.f90). Also, the GRIB I/O has been extended in order to use the
eccodes library for I/O.

In this context we also started a general re-factoring of TSA. Up to now we
removed the old COSMO (src_soil_multlay.f90, data_soil.f90) and introduced the new
ones (sfc_xxx.f90). We renamed the TSA-modules terra_<xxx>.f90 to tsa_<xxxyyy>.f90,
as terra (TERRA) stands for the soil model itself.

In the following we sketch the new source code structure and explain how to compile
the TSA binary, which is now called tsa_exec


Source Code of TSA:
-------------------

A listing of the directory "src" shows all files used by TSA:

 - data_block_fields.f90
 - data_constants.f90
 - data_fields.f90
 - data_io.f90
 - data_modelconfig.f90
 - data_parallel.f90
 - data_runcontrol.f90
 - fxtr_definition.f90
 - kind_parameters.f90
 - meteo_utilities.f90
 - sfc_terra_data.f90
 - sfc_terra.f90
 - sfc_terra_init.f90
 - sfc_utilities.f90
 - support_datetime.f90
 - tsa_data.f90
 - tsa_gribio.f90
 - tsa_input.f90
 - tsa_interpol.f90
 - tsa_lmparam.f90
 - tsa_main.f90
 - tsa_output.f90
 - tsa_setup.f90
 - tsa_sfc_interface.f90
 - turb_data.f90
 - utilities.f90


Modules from the COSMO-Model:
-----------------------------

The following modules are taken "as is" from the COSMO-Model Version 5.07
(with only one exception):

 - data_block_fields.f90
 - data_constants.f90
 - data_fields.f90
 - data_io.f90
 - data_modelconfig.f90
 - data_parallel.f90
 - data_runcontrol.f90
 - kind_parameters.f90
 - meteo_utilities.f90
 - sfc_terra_data.f90
 - sfc_terra.f90
 - sfc_terra_init.f90
 - sfc_utilities.f90     : this has changed compared to the original COSMO version
 - turb_data.f90
 - utilities.f90


Modification in sfc_utilities.f90:
  sfc_utilities.f90 uses two variables from the COSMO modules sfc_seaice.f90 and
  sfc_flake_data.f90. Because we did not want to implement these modules also in TSA,
  these two variables are defined directly in sfc_utilities.f90.

Going to a new version of COSMO should be rather easy now by just updating the 
above modules.

Other Support Modules (from fieldextra):
----------------------------------------

 - fxtr_definition.f90
 - support_datetime.f90


Modules for TSA:
----------------

All modules for TSA, which are not imported from elsewhere, now have the prefix "tsa_"
for identification. In the following we list these modules and their "tasks".

 - tsa_data.f90
   contains all additional data (fields, switches, etc.) which do not come from
   the COSMO-Model. The COSMO data-modules should be taken without any changes.

 - tsa_gribio.f90   (former module gribio.f90)
   contains subroutines for processing of Grib:
    => initializes, resets, adds fields and clears "griblist"- a data structure which 
       points to all data fields and information on the PDS for GRIB1 (libDWD) or
       shortname, etc., for eccodes GRIB1 and GRIB2.

    => read_grib (read_grib_eccodes) reads a GRIB file according to "griblist"
    => get_gribinfo, get_grib_info, get_grib_info_icon: to read GRIB information


 - tsa_input.f90   (part of former module terra_io.f90)
   contains subroutines for GRIB input operations
    => read_metforc
    => read_statbin
    => read_const_fields
    => read_const_fields_icon (NEW)
    => read_initial_fields
    => read_initial_fields_icon (NEW)
    => read_lmgrib
    => read_icongrib (NEW)
    => read_radolan

 - tsa_interpol.f90   (former module terra_interpol.f90)

 - tsa_lmparam.f90    (former module terra_lmparam.f90)
   compared to the former module terra_lmparam.f90, new routines are included:
    => parturs_newblock: routine parturs_new, but in blocked data structure
    => stats    (from terra_io.f90)
    => model_abort  (from environment.f90: to avoid usage of this COSMO module)

 - tsa_main.f90       (former terra_TSA.f90)
   Main program for TSA.

 
 - tsa_output.f90
   contains subroutines for output operations
    => grbout
    => grbout_eccodes
    => binout
    => ascout

 - tsa_setup.f90   (former module terra_lmenv.f90 with routine read_namelist)
   contains subroutines that handle the setup of TSA
    => read_namelist  (from former module terra_io.f90)
    => allocate_fields:
    => allocate_block_fields:   new for the blocked data structure
    => init_variables:
    => clean_up:



 - tsa_sfc_interface.f90
   this module is new for the blocked data structure handling and calling TERRA


Compilation:
------------

The compilation procedure has been adapted to the COSMO-Model:

 - Makefile now is no more contained in the "src" directory, but one level up 
   It includes two other files 
     => Fopts:    File with compiler call, options, libraries for linking, etc.
                  This file can be adapted to the special compiler / machine used
                  Some example Fopts file are given in the directory LOCAL

     => ObjDependencies:   Compilation dependencies (made by hand, not automatically!)

   The list of Object Files is included in Makefile


Binary:
-------

The name of the binary is "tsa_exec" and can be started in the run-scripts
