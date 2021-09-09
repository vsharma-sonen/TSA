#!/bin/bash -l

#######################################################################
#####                                                             #####
#####   performs merging of sub-regions back to the main domain   #####
#####   activated by "main_merge.sh                               #####
#####   requires fieldextra.                                      #####
#####                                                             #####
#####   make sure to change settings according to your machine    #####
#####   settings and directories here are set to lema at cscs.ch  #####
#####                                                             #####
#####   the regrid data are set for a small domain                #####
#####   for testing. commented data are for a larger              #####
#####   domain over central europe.                               #####
#####                                                             #####
#####   written by Julian Todter of GUF                           #####
#####   editted by Yiftach Ziv of IMS (XYZ): zivy@ims.gov.il      #####
#####                                                             #####
#######################################################################

# modules needed:
#module load grib_api

ulimit -s unlimited

# general settings:
dir=$1
cur_list=$2
maxnum=$3
resultDir=$4

cur_dir=${dir}/dir_${cur_list}
mkdir -p $cur_dir
cd $cur_dir

echo $dir $cur_list $maxnum $resultDir $cur_dir

running_date=`cat ../${cur_list}`

# writing merging script for each sub-region:
for outdate in $running_date ; do
	date=`echo $outdate | cut -c 4-15`   #date=yyyymmddHHMM
        echo $date :
	num=1

        # starting with 1st sub-region:
	cat > merge_${num}.nl << end1
&RunSpecification
  strict_usage          = .false.
  verbosity             = "silent"
  additional_diagnostic = .false.
/

&GlobalResource
!  dictionary            = "/users/tsm/project_escha/fieldextra/v12_6_0/resources/dictionary_cosmo.txt"
   dictionary            = "/users/regenass/Source/Soil_Scripts/dictionary_TSA.txt"
   grib_definition_path  = "/users/tsm/project_escha/fieldextra/v12_6_0/tools/support/grib_api_definitions_2"
/

&GlobalSettings
!  originating_center    = "zurich"
  default_dictionary    = "cosmo"
  default_model_name    = "cosmo"
!  experiment_tag        = 100
/

&ModelSpecification
  model_name           = "cosmo"
  earth_axis_large     = 6371229.
  earth_axis_small     = 6371229.
/

&Process
  in_file = "${dir}/reg_${num}/out${date}"
! values from namelist * 10^6: st_rlon(pos),st_rlat,end_rlon,end_rlat,dlon,dlat,pollon(pos),pollat   (pos)=always positive=360-value
  out_regrid_target = "rotlatlon,353200000,-4400000,4770000,3330000,10000,10000,190000000,43000000"
!  out_regrid_target = "rotlatlon,353200000,-4400000,353690000,-3910000,10000,10000,190000000,43000000"
  out_regrid_method="next_neighbour,square,2"
  out_file = "${resultDir}/out${date}_merge_${num}"
!  "${outdir}/out${date}_merge_${num}"
  out_type = "GRIB1"
/

&Process in_field = "T_SO", in_regrid_method="next_neighbour,square,2", in_regrid_all=.true./
&Process in_field = "W_SO", in_regrid_method="next_neighbour,square,2", in_regrid_all=.true./
&Process in_field = "W_SNOW", in_regrid_method="next_neighbour,square,2", in_regrid_all=.true./
&Process in_field = "W_I", in_regrid_method="next_neighbour,square,2", in_regrid_all=.true./
&Process in_field = "RUNOFF_S", in_regrid_method="next_neighbour,square,2", in_regrid_all=.true./
&Process in_field = "RUNOFF_G", in_regrid_method="next_neighbour,square,2", in_regrid_all=.true./
&Process in_field = "LHFL_S", in_regrid_method="next_neighbour,square,2", in_regrid_all=.true./
&Process in_field = "SHFL_S", in_regrid_method="next_neighbour,square,2", in_regrid_all=.true./
&Process in_field = "RAIN_GSP", in_regrid_method="next_neighbour,square,2", in_regrid_all=.true./
&Process in_field = "SNOW_GSP", in_regrid_method="next_neighbour,square,2", in_regrid_all=.true./
end1

/users/tsm/project_escha/fieldextra/v12_6_0/bin/fieldextra_gnu_opt merge_${num}.nl
#exit 0
rm fieldextra.diagnostic
rm merge_${num}.nl

# for all other sub-regions:
for num in `seq 2 $maxnum`; do
   let prevnum=$num-1

		cat > merge_${num}.nl << end1
&RunSpecification
  strict_usage          = .false.
  verbosity             = "silent"
  diagnostic_length     = 100
  additional_diagnostic = .false.
/

&GlobalResource
!  dictionary            = "/users/tsm/project_escha/fieldextra/v12_6_0/resources/dictionary_cosmo.txt"
  dictionary            = "/users/regenass/Source/Soil_Scripts/dictionary_TSA.txt"
  grib_definition_path  = "/users/tsm/project_escha/fieldextra/v12_6_0/tools/support/grib_api_definitions_2"
/

&GlobalSettings
!  originating_center    = "zurich"
  default_dictionary    = "cosmo"
  default_model_name    = "cosmo"
!  experiment_tag        = 100
/

&ModelSpecification
  model_name           = "cosmo"
  earth_axis_large     = 6371229.
  earth_axis_small     = 6371229.
/

&Process
  in_file = "${dir}/reg_${num}/out${date}"
  out_type = "INCORE"
/
&Process in_field = "T_SO", tag = "T_SO_new" /
&Process in_field = "W_SO", tag = "W_SO_new" /
&Process in_field = "W_SNOW", tag = "W_SNOW_new" /
&Process in_field = "W_I", tag = "W_I_new" /
&Process in_field = "RUNOFF_S", tag = "RUNOFF_S_new" /
&Process in_field = "RUNOFF_G", tag = "RUNOFF_G_new" /
&Process in_field = "LHFL_S", tag = "LHFL_S_new" /
&Process in_field = "SHFL_S", tag = "SHFL_S_new" /
&Process in_field = "RAIN_GSP", tag = "RAIN_GSP_new" /
&Process in_field = "SNOW_GSP", tag = "SNOW_GSP_new" /
&Process in_field = "T_G", tag = "my_mask" /

&Process
  in_file = "${resultDir}/out${date}_merge_${prevnum}" 
  !${dir}/dir_${cur_list}/out${date}_merge_${prevnum}"
  out_file = "${resultDir}/out${date}_merge_${num}"
  !"${dir}/dir_${cur_list}/out${date}_merge_${num}"
  out_type = "GRIB1"
/
&Process in_field = "T_SO", set_trange_type="none", merge_with="T_SO_new", merge_mask="my_mask>230,my_mask<330" /
&Process in_field = "W_SO", set_trange_type="none", merge_with="W_SO_new", merge_mask="my_mask>230,my_mask<330" /
&Process in_field = "W_SNOW", set_trange_type="none", merge_with="W_SNOW_new", merge_mask="my_mask>230,my_mask<330" /
&Process in_field = "W_I", set_trange_type="none", merge_with="W_I_new", merge_mask="my_mask>230,my_mask<330" /
&Process in_field = "RUNOFF_S", set_trange_type="none", merge_with="RUNOFF_S_new", merge_mask="my_mask>230,my_mask<330" /
&Process in_field = "RUNOFF_G", set_trange_type="none", merge_with="RUNOFF_G_new", merge_mask="my_mask>230,my_mask<330" /
&Process in_field = "LHFL_S", set_trange_type="none", merge_with="LHFL_S_new", merge_mask="my_mask>230,my_mask<330" /
&Process in_field = "SHFL_S", set_trange_type="none", merge_with="SHFL_S_new", merge_mask="my_mask>230,my_mask<330" /
&Process in_field = "RAIN_CON", set_trange_type="none", merge_with="RAIN_CON_new", merge_mask="my_mask>230,my_mask<330" /
&Process in_field = "RAIN_GSP", set_trange_type="none", merge_with="RAIN_GSP_new", merge_mask="my_mask>230,my_mask<330" /
&Process in_field = "SNOW_CON", set_trange_type="none", merge_with="SNOW_CON_new", merge_mask="my_mask>230,my_mask<330" /
&Process in_field = "SNOW_GSP", set_trange_type="none", merge_with="SNOW_GSP_new", merge_mask="my_mask>230,my_mask<330" /
end1
                /users/tsm/project_escha/fieldextra/v12_6_0/bin/fieldextra_gnu_opt merge_${num}.nl
		rm fieldextra.diagnostic
		rm merge_${prevnum}.nl
		rm ${resultDir}/out${date}_merge_${prevnum}
	done

#        rm merge_${num}.nl

#	sleep 2
	mv out${date}_merge_${maxnum} ../${resultDir}/
done

#rm -rf ${cur_dir}
#rm ${cur_list}
