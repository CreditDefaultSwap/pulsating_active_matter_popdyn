#!/bin/bash -l

##############################
#open MP
##############################

#SBATCH -J dbug_nonintLy 
#SBATCH -o dbug_nonintLy_%j.out
#SBATCH -N 1
#SBATCH --ntasks-per-node=128
#SBATCH --time=1-00:00:00 
time_max_unit=2870m; time_step_unit=0.2m; #for controlling sleeping function and sigkill below

#wp: batch =normal passive submission; interactive = interactive submission
#SBATCH -p batch
#SBATCH --qos=normal

##############################

#Fed externally
s=$1
save_folder=$2
save_file=$3 
log_file=$4
load_file=$5

output_file="run_end"$s".dat"

#g_log_file="g_log_"$s".txt"
g_log_file=`printf "g_log_%0.2f.txt" $s`

#load intel libraries
#module load toolchain/intel/2020b

exe="cloning_def.exe"
exename=${exe%.exe}
echo "Executing cloning for s=$s"
echo "--Using $exe"

######basic variables set up##############
#----physical params
 export Lx=3.0
 export Ly=3.0
 export rho=1.2

 export mu_r=1
 export mu_phi=1

 export lambda=0.05
 export omega=10
 export amplitude_phi=0. #initial value of the control param

#----simulation
 export SEED=51234
 export TMAX=5.

 export t_save=0 #time to start saving from
 export t_clonelog_save=9999 #if saving, set this value <= t_save

 export dt=0.0005
 export tau=200
 export INITSIM=30
 export dumpPeriod=0.1

#----calibration
export controlled_dynamics=1
export calibration_file="amph_curve_L0303_1.2.txt"
export g_calibration_log=$g_log_file

#----Cloning 
 export NC=1000
 export SVALUE=$s 

 export nc_i=20
 export nc_f=30
 export nc_step=300

#----files
 export FILE=$output_file
 export CLONES_DIRECTORY=$save_folder 
 export SAVE_FILE=$save_file 
 export LOAD_FILE=$load_file
 export save_trajectories=0 #set to 1 to save trajs
 export save_cloneslog=0 #set to 1 to save trajs 

####################

#Executes code and saves to log
echo "--Now executing to log_file: " $log_file
#time ./$exe 
time ./$exe > $log_file & 

#wait
pidexe=`pgrep $exename`
#echo "pid for $exe is $pidexe"

bcl="bc -l $HOME/.bc/funcs.bc" 
#removes unit for sleep command; single character of m=mins, h=hrs
time_step=${time_step_unit%?}
time_max=${time_max_unit%?}
# number of times to call sleep
time_counter_limit=`echo "floor($time_max/$time_step)" |$bcl` #must be integer; bcl = .bc/funcs = expanded bc functionality
#echo $time_counter_limit 

time_counter=0; 
while $( pgrep $exename > /dev/null ); do 
	sleep $time_step_unit
	echo "sleeping $time_step"
	time_counter=$(($time_counter+1))
	
	if [ $time_counter -ge $time_counter_limit ]; then
		#process still running but time is ran out #gracefully exit #sigterm signal=15
		kill -15 $pidexe
	fi 
done

#Clean up .tmp. files in traj
echo "--cleaning up for $s"
#rm $save_folder"/"*".tmp.dat" #gets number of files too many error
if (( $save_trajectories || $save_cloneslog ))
then
	rm $save_folder"/0"[0-4]????".tmp.dat"
	rm $save_folder"/0"[5-8]????".tmp.dat"
fi

echo "--All done" 

