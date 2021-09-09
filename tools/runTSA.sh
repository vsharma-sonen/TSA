#!/bin/bash

#######################################################################
#####                                                             #####
#####   runscript for TERRA Stand Alone.                          #####
#####   creates namelist, script for each sub-region              #####
#####   and executes them.                                        #####
#####                                                             #####
#####   make sure to change settings according to your machine.   #####
#####   settings and directories here are set to lema at cscs.ch  #####
#####                                                             #####
#####   the namelist variables are set for a small domain         #####
#####   for testing. commented variables are for a larger         #####
#####   domain over central europe. when running a larger domain  #####
#####   change subreg_x/y to the maximum nodes available          #####
#####                                                             #####
#####                                                             #####
#####   written by Yiftach Ziv of IMS (XYZ): zivy@ims.gov.il      #####
#####                                                             #####
#######################################################################


subreg_x=5                           # number of subregions in x axis for external parallelization
subreg_y=3                            # number of subregions in y axis for external parallelization
let reg_num=$subreg_x*$subreg_y      # total number of sub regions for external parallelization

newdir=REG     # directory for each sub region runscript. !!!Change accordingly
               # when re-running, do not use same $newdir. give a new name and delete (rm -rf) all old $newdir directories
rm -rf oldREG*             # Change accordingly if necessary

#dirName=itsik/Yiftach/TSA/exportPack       # Change accordingly
maindir=$SCRATCH/TSA_Package       # Change accordingly
#cd $maindir

#maindir_workspace=/workspace/${dirName}
#mkdir -p $maindir_workspace
#ln -s $maindir_workspace out           # this will be the output folder. make sure you have enough disc space
maindir_workspace=${maindir}/out
rm -rf $maindir_workspace
mkdir $maindir_workspace		# this will be the output folder. make sure you have enough disc space
					# last 3 commented rows will create output folder in workspace (of lema)

cp UTIL/main_merge.sh out/                  # these scripts are for merging the sub-regions at post-processing
cp UTIL/merge_gribs_domains.sh out/         # make sure they are in the right folder, or copy them yourself

for reg in `seq 1 $reg_num`; do           # creates directories, runscripts and executes them for each sub-region
        echo processing region $reg out of $reg_num ...
        mkdir ${newdir}_${reg}
        cd ${newdir}_${reg}
        touch run_${reg}.job
        chmod u+x run_${reg}.job
        outdir=$maindir_workspace/reg_${reg}
        rm -rf $outdir
        mkdir -p $outdir
        ln -s $outdir out
        cp ../YUSPECIF .                  # this is needed because of a bug that requires YUSPECIF file to pre-exist
                                          # make sure it exists in the right folder for copying
        # writing the sub-region run script:
        # Change settings according to your machine
        cat > run_${reg}.job << end1      
#!/bin/tcsh
#SBATCH --job-name=TSA${reg}
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --output=out/LOG_reg_${reg}.log
#SBATCH --partition=long
#SBATCH --time=120:00:00

rm YU* 

# set environmental parameters
setenv MALLOC_MMAP_MAX_ 0
setenv MALLOC_TRIM_THRESHOLD_ 536870912

# Run LM in case directory
aprun -n 1 ${maindir}/SRC/terra
end1


       # Writing Namelist files for each sub-region:

	cat > INPUT_TERRA << end2
&RUN_TERRA
ydate_ini="2014010100"
ydate_end="2014010200"
!ntstep_max=480
dlon=0.02
dlat=0.02
ie=10
je=6
startlon=-5.7
startlat=-4.0
!ie=520
!je=350
!startlon=-5.7
!startlat=-4.0
pollon=-170
pollat=43
polgam=0.0
nsub_x=${subreg_x}
nsub_y=${subreg_y}
which_subreg=${reg}
dt=180
ke_soil=7
czhls_const= 0.01, 0.03, 0.09, 0.27, 0.81, 2.43, 7.29, 21.87
lvegadapt=.true.
lcheck=.false.
lmelt=.true.
lmelt_var=.true.
lcalc=.true.
lconstvegalb=.true.
lstomata=.true.
itype_evsl=2
itype_hydbound=1
lz0local=.false.
itype_hydparam=1
itype_heatcond=2
itype_root=2
itype_hydcond=1
kexpdec=2.0
crsmin=150.0
/END
&EXTPARA
lhomosoil=.false.
constfilename='${maindir}/ext/external_calmo2_filtered'
lgettcl=.false.
/END
&SOILINIT
lhomoinit=.false.
soilinitdir='${maindir}/init/'
soilinitprefix='ini_'
lmulti_in=.true.
/END
&METFORCING
metfiledir='${maindir}/forc/'
metfileprefix='ter'
tincr_max=0
lhourly_data=.true.
ldestaggeruv=.false.
ntype_atminput=1
ntype_radinput=1
lpar=.true.
ntype_raininput=1
ke_model=60
dz=10
dz_u=10
/END
&OUTPUT
outdir='out/'
outprefix='out'
nout_interval=120
ntype_output=3
/END
end2

        # Wait and run the sub-region runscript
        sleep 5
	sbatch run_${reg}.job
	cd ../
done
