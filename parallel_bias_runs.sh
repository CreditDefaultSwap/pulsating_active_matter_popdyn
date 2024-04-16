#!/bin/bash
#Script file to launch cloning sims for different values of the bias parameter 's'

#---choose bias parameters
s_start=0.00; s_end=-1.00; ds=-1.0; 

load_file_flag=0 #if loading from file

#Folder where trajectories would be saved
traj_folder="controlled"  

#-------------rest of the code
s=0; #needs to be set as a variable for loop below

bcl="bc -l $HOME/.bc/funcs.bc"  #expanded bc set; install separately

s_array=( ); #array
ns_points=`echo "floor(($s_end-($s_start))/$ds)+1" |$bcl` #must be integer; bcl = .bc/funcs = expanded bc functionality

#echo $ns_points

for (( i=0; i<$ns_points; i++))
do
	s=`echo "$s_start + $ds*$i"| bc -l`
	s_array+=( $s ) #appends to array; + should be next to array name, no space
	#echo ${s_array[$i]} 
	
	#folders & save file
	save_folder=$traj_folder"/traj_"$s    
	mkdir -p $save_folder	

	#save_file="cloning_"$s".dat"
	save_file=`printf "cloning_%0.2f.dat" $s`
	#output log
	#log_file="logfile_"$s".txt";
	log_file=`printf "logfile_%0.2f.txt" $s`
	
	echo $save_file $log_file

	echo "Now submitting slurm job for $s bias"
	#Now pass down to the .slurm script where it will request as a separate job for each s
	if [ $load_file_flag -eq 1 ]; then
		load_file=`printf "cloning_%0.2f.dat" $s`
	    save_file=`printf "cloning_pl_%0.2f.dat" $s`
		echo "*Run with loading from file $load_file"
		./run_setup.slurm $s $save_folder $save_file $log_file $load_file  &
	else 
		./run_setup.slurm $s $save_folder $save_file $log_file  &
	fi
done

#wait 


echo "Done at bash script level!"

#das it?yes |3^) 

