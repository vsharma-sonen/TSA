#!/bin/bash

#######################################################################
#####                                                             #####
#####   main script to merge sub-regions                          #####
#####   back to main domain.                                      #####
#####   activates "merge_gribs_domains.sh"                        #####
#####                                                             #####
#####   make sure to change settings according to your machine    #####
#####   settings and directories here are set to lema at cscs.ch  #####
#####                                                             #####
#####   written by Yiftach Ziv of IMS (XYZ): zivy@ims.gov.il      #####
#####                                                             #####
#######################################################################

module load grib_api

maindir=${SCRATCH}/TSA_Package/out     # Change Accordingly
cd $maindir
resultDir=${maindir}/result
mkdir ${resultDir}

numregs=`ls | grep reg | wc -l`

# Looping over time span of the run month by month.
for year in 2011; do         # Change Accordingly
       for mm in 11 12; do     # Change Accordingly
#            if [ year -eq 2015 ] && [ mm -gt 1 ]; then  continue;  fi   # in order to avoid running redundant months
                cd $maindir/reg_1/
		ls -1 out${year}${mm}* > ../list_${year}${mm}
		cd ../ 

	        cat > merge_gribs_par_${year}${mm}.job << endjob	
#!/bin/tcsh
#SBATCH --job-name=M${year}${mm}
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --output=LOG_merge_${year}${mm}.log
##SBATCH --partition=pp-long
##SBATCH --time=120:00:00
#SBATCH --partition=postproc
#SBATCH --time=02:00:00


# set environmental parameters
setenv MALLOC_MMAP_MAX_ 0
setenv MALLOC_TRIM_THRESHOLD_ 536870912

# Run LM in case directory
srun -n 1 ./merge_gribs_domains.sh ${maindir} list_${year}${mm} ${numregs} ${resultDir}
endjob

                # execute monthly merging: merge_gribs_domains.sh
		sleep 5
		sbatch merge_gribs_par_${year}${mm}.job
	done
done

