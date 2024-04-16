#include <iostream>
#include <fstream>
#include <random>
#include <chrono>
#include <vector> 
#include <algorithm> 
#include <list>

#ifdef _OPENMP
#include <omp.h>
#endif

#include "env.hpp"
#include "particle.hpp"
#include "maths.hpp"
#include "spline.h"
//signal catching library
#include <csignal>

//-----!!!Global pointers: needed for sigcatch procedure--- 
double* lnX_ptr, * lnX_before_ptr, * lnX_after_ptr; 
bool*   lnX_set_ptr; 
bool * pull_push_finished_ptr, * pull_push_copy_ptr; 
bool * clone_loop_finished_ptr, * control_finished_ptr;
double* iter_max_ptr; 
int* pushOffset_ptr, * pullOffset_ptr;
int* newClones_ptr; //array is just a pointer so can use any int pointer
//-------------------

#include "cloningserial.hpp"

//-----!!!Global pointers variables etc!!!
CloningSerial<Deformables>* clones_ptr;
int* N_ptr;
double* sValue_ptr; 

std::chrono::system_clock::time_point startTime_glb = std::chrono::system_clock::now(); 
//-------------------


//wp:It takes a 'system' class type which denotes the simulation system
//wp: It has been templated to accept any system class, in < argument
template<> void CloningSerial<Deformables>::
	loadState(Deformables* dummy) { ///wp must be specialized for each case 
	std::cout << "## Now loading prev. run, confs etc from file " << loadInput.getInputFile() << std::endl;
	// Load cloning configurations from input file. //wp: based on how saveState is performed 
	
  	//Parameters
	double timeTraj, amphiControl;
    //number of particles in system // gets Parameters_def object then access the num of particles
	int numberParticles=(dummy->getParameters()).getNumberParticles(); 
	//std::cout <<"# particles:" << numberParticles << std::endl;

	//Calculate binary size info
	//wp:Header variables solved
	long int headerLength = 7*sizeof(int) + 4*sizeof(double); 
	// size of system related variables, including pos + orientations
	long int headerParticles = sizeof(std::mt19937);
	//std::cout << "Size of mt19937 "<< headerParticles << std::endl;
	long int particleLength = 3*sizeof(double); //wp: one particle 
	long int particlesLength = numberParticles*particleLength; //wp: All particles
	long int footerParticles =  3*sizeof(double); //wp: traj obsservables + time 
	long int cloneLength = headerParticles+particlesLength+footerParticles;
	long int clonesLength = 2*nc*cloneLength; //wp: 2 nc due to zizag arrswitch implmentation
	long int fileLength = headerLength+ clonesLength; //wp: to make sure file matches what's saved 
	
	if ( loadInput.getFileSize() != fileLength ) { 
		std::cerr << "***Filelength: " << fileLength << " filesize: " <<  loadInput.getFileSize() << std::endl;
		std::cerr << "***loadState(): File size does not match expected size. Exiting" << std::endl;  exit(1); } 

	//wp: Now start reading and saving back to variables
	for (int i=0; i<6; i++){ //wp: will read nc,nc_i,etc, from the constructor output
		 //loadInput.read<int>();
		 std::cout<<loadInput.read<int>() << std::endl;
	} 
	loadInput.read<double>();

	//wp: header
	arrswitch = loadInput.read<int>();
	std::cout << " " << arrswitch << std::endl;

 	timeTraj=loadInput.read<double>();
	outputPsiOffset[0] = timeTraj; //wp: trajectory time so far 
	outputPsiOffset[1] = loadInput.read<double>(); //wp: outputPsi/trajectory time

	amphiControl=loadInput.read<double>();//wp: amphi for control force strenght; It can be 0 i.e. no control 
	
	std::cout << "scgf " << outputPsiOffset[1]/numberParticles << std::endl;
	std::cout << "timeTraj, outputPsi " << timeTraj << " " << outputPsiOffset[1] << std::endl;
	std::cout << "amphi control " << amphiControl << std::endl;

	if(save_particles){
		//wp: first create template clones and then fill them later with actual loaded states 
		for (int i=0;i<2*nc;i++) {//system from dummy --> constructor 2
			//systems[i] = new Deformables(dummy, -1, tau, tau, cloneFilename(i)); //ignore seeding as it will copy the previous generator 
			systems[i] = new Deformables(dummy,cloneFilename(i),tau); //ignore seeding as it will copy the previous generator 
			//std::cout << "Will be written to " << cloneFilename(i) << std::endl; 
		}
	}else{ //no particles or clonelog saving --> constructor 2
		//wp: first create template clones and then fill them later with actual loaded states 
		for (int i=0;i<2*nc;i++) {//system from dummy 
			systems[i] = new Deformables(dummy,tau); //ignore seeding as it will copy the previous generator 
			//std::cout << "Will be written to " << cloneFilename(i) << std::endl; 
		}
	}

	//std::cout << " Now updating clones" << std::endl;
	//wp: Now update clones
  	for (int i=0;i<2*nc;i++) { 
		//std::cout << " Updating clone: " << i <<  std::endl;
		// RANDOM GENERATOR
		systems[i]->setGenerator(loadInput.read<std::mt19937>());

		//std::cout << " Before configuration copies: " << i <<  std::endl;
		// COPY POSITIONS AND ORIENTATION
		for (int j=0; j < numberParticles; j++) {
			for (int dim=0; dim < 2; dim++) {
				(systems[i]->getParticle(j))->position()[dim] = loadInput.read<double>(); 
					//printf("%e ",(systems[i]->getParticle(j))->position()[dim]);
			}
				(systems[i]->getParticle(j))->orientation()[0] = loadInput.read<double>(); 
			//printf("%e\n",(systems[i]->getParticle(j))->orientation()[0]);
		} //wp: num of particles loop

		//std::cout << " After configuration copies before sites: " << i <<  std::endl;
		//wp: writes initial positions and orientations to file
		systems[i]->saveInitialState();   //wp: for writeTraj formatting expectations 
		//if( i < nc) systems[i]->saveInitialState(); //wp: crappy way to keep the nc_i > nc clones behind time progression in the file (due to burn in) 

		//wp: reconstructs sites from stored positions
		systems[i]->reconstructSites();
		
		//std::cout << " After sites: " << i <<  std::endl;
		//wp: observables
		systems[i]->setCurrentTraj(loadInput.read<double>());   // 1/N/omega \sim dot(theta) 
		systems[i]->setCurrent2Traj(loadInput.read<double>());  // 1/N/omega \sim dot(theta)^2 
		systems[i]->setTimeTraj(loadInput.read<double>());  // time performed in trajectory so far 
		//printf("%e %e %e\n",systems[i]->getCurrentTraj(),systems[i]->getCurrent2Traj(),systems[i]->getTimeTraj()); 
		
		//std::cout << "Now saving sync strength" << std::endl;
		(systems[i]->getParametersPointer())->setAmplitude_phi(amphiControl);  //synchronization strength; can be 0
	} //end clone loop
	//wp: if wanted, offset every nth=1000 clones with a random call, thereby 'desynchronizing' any possibly generators with same seed 
	outSyncRandomGenerator();
	printf("\t -->Done loading file.\n");
//wp: das it
}

template<> void CloningSerial<Deformables>::
	saveState() { 
	// Save cloning configurations to output file.

	// CLONING ALGORITHM PARAMETERS
	//saveOutput.write<double>(iter*tau*systems[0]->getTimeStep()); //wp: time total up to last iteration
	saveOutput.write<int>(arrswitch); //wp: which track state is supposed to be used (if i or i+nc)

	//wp: average of timeTraj
	double timeTraj_ave=0;
	//int npoints=2
	//std::vector<double> output(npoints, 0.0);
	for (int i=0; i < nc; i++) { //wp: adds over all the clones then averages
		timeTraj_ave+=systems[i]->getTimeTraj(); //wp: trajectory time so far
	//	output[0] += (finalSystem(i))->getCurrentTraj()/((finalSystem(i))->getTimeTraj());
	//	output[1] += (finalSystem(i))->getCurrent2Traj()/((finalSystem(i))->getTimeTraj());
	} 
	timeTraj_ave/=nc;
	//for (unsigned int j=0;j<npoints;j++) saveOutput.write<double>(output[j]/nc);

	// CLONING OUTPUT
	saveOutput.write<double>(timeTraj_ave); //wp: trajectory time so far
	saveOutput.write<double>(outputPsi); //wp: estimate of the scaled cum. func. so far //wp: already divided by time

	// wp: control parameter 
	saveOutput.write<double>((systems[0]->getParametersPointer())->getAmplitude_phi()); //wp: control strength; same for all clones

	//wp:debug 
	//std::cout << arrswitch << std::endl;
	//std::cout << timeTraj_ave << std::endl;
	//std::cout << outputPsi << std::endl;

    //number of particles in system: gets the Parameters_def copy then access; not a pointer
	int numberParticles=(systems[0]->getParameters()).getNumberParticles(); 

	// CLONES
	for (int i=0; i < 2*nc; i++) {
		// RANDOM GENERATOR
		saveOutput.write<std::mt19937>((systems[i]->getRandomGenerator())->getGenerator());
		// POSITIONS AND ORIENTATIONS
		for (int j=0; j < numberParticles; j++) {
			for (int dim=0; dim < 2; dim++) {
				saveOutput.write<double>((systems[i]->getParticle(j))->position()[dim]);
					//printf("%e ",(systems[i]->getParticle(j))->position()[dim]);
			}
			saveOutput.write<double> ((systems[i]->getParticle(j))->orientation()[0]);
			//printf("%e\n",(systems[i]->getParticle(j))->orientation()[0]);
		}
		// DUMPS
		//saveOutput.write<int>(systems[i]->getDump()[0]);
		saveOutput.write<double>(systems[i]->getCurrentTraj()); //wp: trajectory averaged observables
		saveOutput.write<double>(systems[i]->getCurrent2Traj());
		saveOutput.write<double>(systems[i]->getTimeTraj()); //wp: trajectory time so far
		//printf("%e %e %e\n", systems[i]->getCurrentTraj(),systems[i]->getCurrent2Traj(),systems[i]->getTimeTraj());

	} //wp: clones loop
	//std::cout<<"Done saving file." << std::endl;
	saveOutput.close();
}

//wp: This will compile the full trajectory of each clone track for every iteration of the cloning algorithm. The files are written in xxxx.xxxxxx.dat 
template<> void CloningSerial<Deformables>::
	writeTrajFiles(Write& clonesLog) {	// Write trajectory files for cloning loop.
	printf("\t Writing final cloning trajectories  (xxxx.yyyyyy.dat)\n"); 

	//printf("Test 0.43 %i, 0.45 %i, 0.46 %i, 1.46 %i, 1.44 %i, 1.45 %i\n", (int) 0.43, (int) 0.45, (int) 0.46, (int) 1.46, (int) 1.44, (int) 1.45);
	//wp: (int) typecasting acts like a floor function

  int period = 1;
  iter=iter-Niter_save; //wp: to account for when only a partial of the trajectory time has been saved, up to the end
	if (iter <= 0 ) return; //wp: return in principle because nothing was saved
  // BUILD TRAJECTORIES

  clonesLog.flush();
  //wp: in CLONES_DIRECTORY/clones.log
  Read log(clonesLog.getOutputFile());
 
  /*printf("reading clone.log: \n");
  for (int t=0; t < iter; t++) {
	printf("; %i;",t);
    for (int i=0; i < nc; i++) {
      printf("%i",log.read<int>());
	}
  }
  printf("done reading\n");
  log.close(); log.open();*/
  	
  std::vector<std::vector<int>> trajectories (nc);
  std::vector<int> parents (nc);
  for (int i=0; i < nc; i++) { parents[i] = i; }

  //printf("\t iter %i\n", iter);
  //wp: iter = # of time slices
  for (int t=0; t < iter; t++) {
    for (int i=0; i < nc; i++) {
      int newParent = log.read<int>( -((t + 1)*nc - parents[i])*sizeof(int), std::ios_base::end); //wp: reads from end
	  //printf("t, nc_i, newParent %i %i %i \n", t,i,newParent);
	  //wp: trajectories is an array that reconstructs parent geneology starting from present to past
	  //wp: Note that this, by definition, means is that it is not a simple reading of the cloning but a backwards reconstruction
      trajectories[i].insert(trajectories[i].begin(), newParent);
      parents[i] = newParent;
    }
  }
  log.close();

  //wp: debug: this prints out the descendance line, starting from latest parent up to beginning i.e. should match the clone #
  // 	For example, it could print 5 5 4 3 2 1 for nc_i=1, meaning the parent of clone track 1 was 5 at last point, then 5 again, then 4, ..., then 2 and ultimately 1. So the full parent trajectory, moving forwards in time is 1 2 3 4 5 5
  //for (int i=0; i < nc; i++){
  //  printf("traj: nc_i %i: t, parent \n", i ); 
  //	for (int t=0; t< iter; t++){
  //  	printf("\t\t\t%i %i\n", t, trajectories[i][t]);
  //  }
  //}

  // WRITE .dat FILES
  // This dat structure also has methods to access positions and orientation of traj file 
  std::vector<Dat_def*> dat;
  std::vector<std::vector<double>> currentAve;

  //std::vector<std::vector<double>> activeWork;

  for (int i=0; i < 2*nc; i++) {
	//wp:makes structure 'dat' that stores vectors containing all calculated work/order params etc for entire run on this loop 
	//printf("Before flush i %i clone\n", i);
    systems[i]->flushOutputFile();
    dat.push_back(new Dat_def(systems[i]->getOutputFile(), true));
    dat[i]->close();
    //activeWork.push_back(dat[i]->getActiveWork());
    currentAve.push_back(dat[i]->getCurrents());
  }

	//printf("\t Before writing the xxxx.yyyyyy file\n");
  for (int i=0; i < nc; i++) { 
	//wp: name of final .dat file in a 000r.00000i.dat format
	//wp: Done such that the length is constant i.e. string(#,'char') = char written n times
	//wp: to_string(9).length() = 1 char so e.g. 6-1 = 5 '0' characters below 
    Write output(std::experimental::filesystem::path(
      std::experimental::filesystem::path(clonesDirectory) /
      [](int index, int rIndex) 
        { return
          std::string(4 - std::to_string(rIndex).length(), '0')
            + std::to_string(rIndex) + std::string(".")
          + std::string(6 - std::to_string(index).length(), '0')
            + std::to_string(index) + std::string(".dat"); }
        (i, runIndex)
      ).u8string());

	//printf("AFter the 'experimental' declaration\n");
    // header
    output.write<int>(systems[i]->getNumberParticles());
    output.write<double>(systems[i]->getDensity());
    output.write<double>(systems[i]->getSystemSize());
    output.write<double>(systems[i]->getTimeStep());

    // initial frame //wp: i = which clone position
	//printf("This is i  iter %i %i\n",i,iter);
	dat[i]->getFrames(); //wp: dat is of size 2*nc
	//printf("Frames in file %i\n", dat[i]->getFrames() );
	
	//int initFrame=0; //wp
    int initFrame = std::max((dat[i]->getFrames() - 1) - period*((int) (iter + 1)/2),(long int) 0); //wp: should be 0 at least
	//wp: trajectories contains integers of which clones will be copied for next iteration
    int initClone = trajectories[i][0]; //wp: i = clone track # // the t=0 is actually the last written point in trajectory i.e. moves backwards

	//int nc_test=1; //wp:debug
	//if(i ==nc_test){ 
	//printf("Iter: %i\n", iter);
	//std::cout << output.getOutputFile() << std::endl; printf("frames, initClone, initFrame: %i %i %i\n", dat[i]->getFrames(), initClone, initFrame); //} //wp

    dat[initClone]->open();
	//printf("Before positions and orientations of 1st frame\n");
    for (int p=0; p < systems[i]->getNumberParticles(); p++) { // output all particles
      // POSITIONS
      for (int dim=0; dim < 2; dim++) { // output position in each dimension
        output.write<double>(dat[initClone]->getPosition(initFrame, p, dim));
      }
      // ORIENTATIONS
      output.write<double>(dat[initClone]->getOrientation(initFrame, p)); // output orientation
      // VELOCITIES
      //for (int dim=0; dim < 2; dim++) { // output velocity in each dimension
      //  output.write<double>(dat[initClone]->getVelocity(initFrame, p, dim)); // output velocity
      //}
    }
    dat[initClone]->close();

	//printf("Now doing other frames\n");
    // other frames
    for (int frame=1; frame <= iter; frame++) {
	  int frame_switch=frame;
	  //int frame_switch = loadInput.getInputFile() != "" ? frame+1 : frame;  //wp: if reloading, ensure it begins on the right cloning track (2nc or not)

      int currClone = trajectories[i][frame - 1] + nc*(1 - (frame_switch%2)); //wp:offsets by nc depending on cycle //wp: mod(x,2) = either 1 or 0 for any x
      int currFrame = (dat[currClone]->getFrames() - 1) - ((int) (iter - frame)/2)*period; //wp: (int) acts as a floor function
	
	 //printf("i, currClone, currFrame: %i %i %i\n", i,currClone,currFrame);
	 //if(i ==nc_test) 
		//printf("nc*(1-frame%2) %i\n", nc*(1 - (frame_switch%2)));
		//printf("frames-1 %i int 0.5*(iter-frame) %i, currFrame %i\n",(dat[currClone]->getFrames() - 1), ((int) (iter - frame)/2), currFrame);
	 //if(i ==nc_test) 
		//printf("currClone, currFrame: %i %i\n", currClone, currFrame);

      dat[currClone]->open();
      for (int p=0; p < systems[i]->getNumberParticles(); p++) { // output all particles
        // POSITIONS
        for (int dim=0; dim < 2; dim++) { // output position in each dimension
          output.write<double>(
            dat[currClone]->getPosition(currFrame, p, dim));
        }
        // ORIENTATIONS
        output.write<double>(
          dat[currClone]->getOrientation(currFrame, p)); // output orientation
        // VELOCITIES
        //for (int dim=0; dim < 2; dim++) { // output velocity in each dimension
        //  output.write<double>(
        //    dat[currClone]->getVelocity(currFrame, p, dim)); // output velocity
        //}
      } //wp: number of particles loop
      dat[currClone]->close();

		//wp: Current value
      output.write<double>(currentAve[currClone][currFrame/period - 1]);
	  //if(i ==nc_test){ printf("%f\n", currentAve[currClone][currFrame/period - 1]);}

	  //printf("Right before closing of trajWrite\n");
    } //wp: frames loop

    // close file
    output.close();
  } //wp: clone loop

  // delete pointers to input files
  for (int i=0; i < 2*nc; i++) delete dat[i]; 
}

//wp: Signal catching procedure
void signalHandler( int signum ) {
	std::cout << "**********Interrupt signal (" << signum << ") received.\n";
	//wp: since dt variable, use actual, average value over all clones as an approximation 
	double timeTraj_ave=0, outputPsi;
	int iter_last, iter_max, finalOffset, nc, arrswitch;

	//wp:Gather basic info about simulation
	iter_last=clones_ptr->iter;
	nc=clones_ptr->nc;

	if(iter_last < *iter_max_ptr){ 
		std::cout << "\tSaving for broken cloning loop\n" ;
		arrswitch=clones_ptr->arrswitch;
		if(!*pull_push_finished_ptr && !*pull_push_copy_ptr){ //iteration not finished, so used last generation of clones
			//arrswitch=(arrswitch==0) ? 1:0;
			arrswitch=1-arrswitch;
			//wp: in case it breaks right before or after setting of lnX-->use lnX for fully completed iteration
			if (!*lnX_set_ptr){ 	
				lnX_ptr=lnX_after_ptr;
				std::cout << "	-Saving where cloning loop didn't finish\n";
			}else{ 
				lnX_ptr=lnX_before_ptr;
				std::cout << "	-Saving where cloning finished, but didn't copy\n";
			}
		}else if(*pull_push_copy_ptr && !*pull_push_finished_ptr){ //iteration complete, half way copying. Recopy and then go with new list of clones
			std::cout << "	-Saving where cloning finished, was copying and didn't complete\n";
			#ifdef _OPENMP
			#pragma omp parallel for
			#endif
			for (int i=0; i<nc; i++) {
			  //wp: copy from old clone index "pushOffset" to new clone index "pullOffset"
			  (clones_ptr->systems[*pullOffset_ptr + i])->copyState(clones_ptr->systems[ *pushOffset_ptr+ *(newClones_ptr + i) ]); // clone particles
			  (clones_ptr->systems[*pullOffset_ptr + i])->copyDump(clones_ptr->systems[ *pushOffset_ptr + *(newClones_ptr + i) ]); // clone dumps
			} 

			//wp: will forgo control param update for this step but it's a forgivable estimation, especially at long times where it's basically constant	
			arrswitch=1-arrswitch;	
		}else if(*control_finished_ptr && !*clone_loop_finished_ptr){ //cloning done, just be sure to update arrswitch
			std::cout << "	-Saving where cloning and copying finished, except for arrswitch reassignment\n";
			//wp:arrswitch not yet updated
			arrswitch=1-arrswitch;	
		}//Note that otherwise, arrswitch remains same
		clones_ptr->arrswitch=arrswitch;
	} //else keep whatever was defined in last iteration 

	finalOffset=arrswitch*nc;

	//std::cout << nc << " " << arrswitch << " " << finalOffset << std::endl;
	for (int i=0; i<nc; i++) timeTraj_ave+=(clones_ptr->systems[finalOffset+i])->getTimeTraj();
	timeTraj_ave/=nc; timeTraj_ave-=clones_ptr->outputPsiOffset[0]; //wp: only want time spanned after loading, not entire run time
	//std::cout << timeTraj_ave  << std::endl;

	outputPsi = (*lnX_ptr + (clones_ptr->outputPsiOffset[0])*(clones_ptr->outputPsiOffset[1]))/(timeTraj_ave + clones_ptr->outputPsiOffset[0]);
	clones_ptr->outputPsi = outputPsi; 
	//std::cout << outputPsi  << std::endl;

	// SAVE FINAL CONFIGURATIONS
	//same as above: tell system to begin from likely last finished run
	//Note that this assumes the sigterm will have come *at least* sometime after some cloning cycles have been done
	std::cout << "Writing savefile (if specified) \n";

	if ( (clones_ptr->saveOutput).getOutputFile() != "" ) clones_ptr->saveState();

	//wp: Part that would normally be written if the run had finished normally
	clones_ptr->finalOffset=finalOffset; //used in finalSystem below
	int num_output=2; 

	(clones_ptr->outputOP).assign(num_output, 0.);
	for (int i=0; i < nc; i++) {
		clones_ptr->outputOP[0] += ((clones_ptr->finalSystem(i))->getCurrentTraj())/((clones_ptr->finalSystem(i))->getTimeTraj());
		clones_ptr->outputOP[1] += ((clones_ptr->finalSystem(i))->getCurrent2Traj())/((clones_ptr->finalSystem(i))->getTimeTraj());
	}

	//wp: divides by # of clones
	for (unsigned int j=0;j<num_output;j++) { clones_ptr->outputOP[j] /= nc; }

	//std::cout << "clones.outputPsi " << clones.outputPsi << " " << N  << std::endl;
	
	//wp: Note that this clocks from the start of the code, till the end, rather than just the beginning of doCloning 
	std::chrono::system_clock::time_point endTime = std::chrono::system_clock::now();

	// this horrible syntax just computes elapsed time in seconds, at microsecond resolution
	double outputWalltime;	
 	outputWalltime = 1e-6*std::chrono::duration_cast<std::chrono::microseconds>
	  (endTime - startTime_glb).count(); 

	std::cout << std::endl;
	std::cout << "##last iteration "<< iter_last << std::endl
			  << "##max iteration : "  << *iter_max_ptr << std::endl
			  << "##bias parameter s: "    << *sValue_ptr << std::endl
			  << "#SCGF "  << outputPsi/((*N_ptr)) << std::endl
			  << "#<current> "     << clones_ptr->outputOP[0] << std::endl
			  << "#<current^2> "     << clones_ptr->outputOP[1] << std::endl
			  << "#am_phi final "     << ((clones_ptr->finalSystem(0))->getParametersPointer())->getAmplitude_phi() << std::endl
			  << "##time " << outputWalltime << std::endl; 

	std::cout << "End of sigterm writeup\n" ;
	exit(signum);
}


int main() { 
	//wp: register signal SIGTERM handler; SIGTERM is the 'nicest' signal which allows for program to ditch gracefully  
    signal(SIGTERM, signalHandler); //handler cannot accept arguments 

	// cloning parameters
	double tmax = getEnvDouble("TMAX", 1); // dimensionless time to simulate
	int nc = getEnvInt("NC", 10); // number of clones
	double sValue = getEnvDouble("SVALUE", 0); // biasing parameter
	int seed = getEnvInt("SEED", 0); // master random seed
	int nRuns = getEnvInt("NRUNS", 1); // number of different runs
	int cloneMethod = 2; // should keep this set to 2 (!!)
	int initSim = getEnvInt("INITSIM", 1); // number of initial elementary number of iterations to "randomise" the systems //wp: relaxation steps

	// openMP parameters
	#ifdef _OPENMP
	int threads = getEnvInt("THREADS", -1); // number of threads
	printf("# compiled with openMP\n");
	if ( threads > 0 ) {
		printf("# setting threads %d\n",threads);
		omp_set_num_threads(threads);
	}
	printf("# running on %d threads\n",omp_get_max_threads());
	#endif

	// simulation parameters
	int tau = getEnvInt("tau", 100); // elementary number of steps
	double dt_d = getEnvDouble("dt", 0.001); // time step
	float dt = float(dt_d);

	// simulation
	double temperature = getEnvDouble("temperature", 1); // temperature of the bath 
	
	//wp: calibration file for controlled dynamics
	bool controlled_dynamics= getEnvBool("controlled_dynamics", "false");
	std::string calibration_file= getEnvString("calibration_file", "");
	
	// save trajectories or cloneslog
	std::string clonesDirectory = getEnvString("CLONES_DIRECTORY", ""); // if different than "" then clones trajectories are saved to this directory
	bool save_trajectories=getEnvBool("save_trajectories", "false"); //plain text output
	bool save_cloneslog=getEnvBool("save_cloneslog", "false"); //plain text output
	if(!save_trajectories && !save_cloneslog) clonesDirectory=""; //this will have the effect to not save trajectories, regardless of save time
	//if(!save_trajectories) clonesDirectory=""; //this will have the effect to not save trajectories, regardless of save time

	// load initial cloning state
	std::string loadFile = getEnvString("LOAD_FILE", ""); // if != "": input class from where initial configurations are loaded
	// save final cloning state
	std::string saveFile = getEnvString("SAVE_FILE", ""); // if != "": write class to where final configurations at each call of doCloning are saved

	// VARIABLE DEFINITION
	// data 
	bool dump = true;
	int nc_i = getEnvInt("nc_i",0); //wp: initial clone to turn off 
	int nc_f = getEnvInt("nc_f",nc); //wp: up to final clone to turn off 
	int nc_step = getEnvInt("nc_step",1); //wp: then turns off in steps 

	double t_save = getEnvDouble("t_save",0); //wp: time from which to start svaing trajs. turned into # of iters below 
	double t_clonelog_save = getEnvDouble("t_clonelog_save",0); //wp: time from which to start svaing trajs. turned into # of iters below 
	if( t_save < t_clonelog_save && save_trajectories){ 
		t_save=t_clonelog_save; 
		printf("!!! Can't save trajectories before clonelog starts saving. Setting t_save=t_clonelog_save\n");
	}
	int niter_save = (int) round(t_save/dt_d/tau); //wp: # of iters from which to start saving 
	int niter_save_clone = (int) round(t_clonelog_save/dt_d/tau); //wp: # of iters from which to start saving
	//double dumpPeriod = getEnvDouble("dumpPeriod", 0.1); 

	// system
	double rho = getEnvDouble("rho",1.2); //density of system
	//int Lx = getEnvInt("Lx",4); //# of boxes along x
	//int Ly = getEnvInt("Ly",4); //# of boxes along y 
	double Lx = getEnvDouble("Lx",4); //# of boxes along x
	double Ly = getEnvDouble("Ly",4); //# of boxes along y 

	double mu_r= getEnvDouble("mu_r", 1.); //mu = mobility which also control noise as sqrt(mu..)
	double mu_phi= getEnvDouble("mu_phi", 1.);
	double omega= getEnvDouble("omega", 1.);
	double amplitude_r= getEnvDouble("amplitude_r", 1.); //strength of Uij
	double amplitude_phi= getEnvDouble("amplitude_phi", 1.); // strength of synch 

	double radius=getEnvDouble("radius", 1.);
	double lambda=getEnvDouble("lambda", 0.05); 

	//wp: Define Parameters_def 
	Parameters_def parameters(Lx,Ly,rho,omega,mu_r,mu_phi,amplitude_r,amplitude_phi,radius,lambda,dt,temperature);
	Parameters_def *param_d; //pointer of Parameters_def type
	param_d = & parameters;


	int Nx, Ny, N; 
	Nx=parameters.getNx();
	Ny=parameters.getNy();
	N=parameters.getNumberParticles();

	// output to file //wp:physical and sim params using object from 'Write' class
	std::string filename = getEnvString("FILE", ""); // output file name
	Write output(filename); // output class
	output.write<double>(tmax);
	output.write<int>(nc);
	output.write<double>(sValue);
	output.write<int>(seed);
	output.write<int>(nRuns);
	output.write<int>(cloneMethod);
	output.write<int>(initSim);
	output.write<int>(N);
	output.write<int>(tau);
	output.write<double>(dt_d);
	output.close();

// dummy system
	Deformables dummy(&parameters,seed); 

	printf("## CloningSerial Code: numClones %d runs %d s %.2f \n",nc,nRuns,sValue);
	if(controlled_dynamics) printf("\t With controlled dynamics\n");
	printf("\t tmax %.3f initSim (relaxation) # %d, tau %d Delta.t %.4e\n",tmax,initSim,tau,dt);
	printf("\t Lx %.3f,Ly %.3f, N particles %i\n", Lx, Ly, N); 
	//printf("\t Lx %i,Ly %i, N particles %i\n", Lx, Ly, N); 
	printf("\t rho %.2f , omega %.1f, lambda %.3f, am_phi %.2f\n",rho, omega, lambda, amplitude_phi); 

	printf("\t Saving from t_save %.4f, at nc_i %i, to nc_f %i, with nc_steps %i \n", t_save, nc_i, nc_f, nc_step); 
	if(save_cloneslog) printf("\t Saving clones.log file from %0.2f time\n", t_clonelog_save); 

	///------------cloning objects -----------//wp: class object which will contain cloning stuff
	CloningSerial<Deformables> clones(nc, nc_i, nc_f, nc_step, niter_save, niter_save_clone, save_trajectories, 
									tau, sValue, cloneMethod, clonesDirectory, loadFile, saveFile);	
 	
	//wp:!!----declares global pointers 
 	clones_ptr =& clones;
	N_ptr=& N; sValue_ptr=& sValue;


	printf("\t clonesDir: '%s', loadFile: '%s', saveFile: '%s'\n", clonesDirectory.c_str(), loadFile.c_str(), saveFile.c_str()); 
	//wp: controlled dynamics
	Write g_calibration_log(getEnvString("g_calibration_log", "")); //plain text output


	// *. Load calibration file: Note that this curve depends on system size and must be manually specified given this consideration 
		std::vector<std::vector<double>> input_cal_data;
		input_cal_data=read_file(calibration_file);
		// *. Create spline 
		double cal_x0=0, cal_xmax=0, amplitude_min=0, amplitude_max;
		std::vector<double> X = input_cal_data[0]; std::vector<double> Y = input_cal_data[1];
		cal_x0=X[0], cal_xmax=X.back(); //wp: to respect boundaries later on *inverted* function
		amplitude_min=Y[0]; amplitude_max=Y.back();

		tk::spline spline_input(X,Y,tk::spline::cspline_hermite,true); //wp: true = monotonic, Hermite = less global interference from any one point
		// *. Create x2 refined data which will then serve as the template to be inverted 
		std::vector<double> data_i(2); 
		std::vector<std::vector<double>> data_sp, dataInv; 
		double x0=0, dx=0, x_i, y_i;
		int numpts=2*X.size();

		x0=X[0]; dx=(X.back()-x0)/numpts;
		for (int i=0; i<=numpts; i++){ 
			x_i=x0+i*dx; y_i=spline_input(x_i);
			data_i[0]=(x_i); data_i[1]=(y_i);
			//printf("%i, %e, %e \n",i, x_i, y_i);
			//save to 2D data
			data_sp.push_back(data_i);
		}

		transpose(data_sp);

		// *. Now use data to invert and get final spline function 
		X=data_sp[1]; Y=data_sp[0]; 

		tk::spline spline_inv(X,Y,tk::spline::cspline,true); //wp: true = monotonic 
		//wp----done for splining
	

	// set up the clones
	if ( loadFile == "" ) { clones.init(&dummy, seed); } 
	else { clones.loadState(&dummy); } //wp: load from file

	std::cout << "## master seed " << seed << std::endl;
	printf("## Now starting cloning:..........\n");	

  double sFactor = N*tau*dt;

	for (int run = 0; run<nRuns;run++) {
		//wp: Main cloning algorithm procedure Y_Y
		//wp: !!! The actual code is in the header. Here it just calls the obect and its method so the innards are in cloningserial.hpp
		// go! (this includes generating "random" [different] initial conditions for the clones)
		clones.doCloning(tmax, initSim, 
			//wp: This is the main part of the code to consider:
			//wp: --1. Define new class of system
			//wp: --2. Define weight function
			//wp: --3. Introduce relevant control functions
			
			// ITERATION FUNCTION
			[](Deformables* deformable, int Niter, bool save_traj) { iterate_deformable_WCA_control_strat(deformable, Niter, save_traj); },

			// GET WEIGHT FUNCTION
			[&sFactor, &amplitude_max, &amplitude_min, &controlled_dynamics](Deformables* deformable) { 
				//wp:exp(-sWeight) where sWeight= (-U-sA) and U = controlling forces; A = average quantity of interest
				//wp: s = biasing parameter 'field' (reciprocal temperature :^) 
				double sWeight, obvs_norm; 

				sWeight = sFactor*deformable->getBiasingParameter()*deformable->getCurrent2(); 
				//wp: Order parameter weight function with controlled dynamics 
				//wp: s*tau*N*O - 1/s*tau*N int_t \sum_i (S-S_c); here, controlAve is the integral of (S-S_c) therefore no need to remultiply with 'sFactor' 
				if(controlled_dynamics) sWeight = sWeight - sFactor*deformable->getControlAve(); 

				return sWeight;
			},

			// CONTROL FUNCTION #wp: '\' are line continuations
			[&nc, &spline_inv, &cal_x0, &cal_xmax, &amplitude_max, &amplitude_min, &g_calibration_log, &controlled_dynamics](std::vector<Deformables*>& deformable, int pullOffset, int pushOffset, bool  save_calibration) {
				//wp: This is where to alter the amplitude of the control force amph 
				
				if(controlled_dynamics){//wp:!!! Performs observable calibration/adjustment for a given value of bias s			
		
				// order parameter  //wp: depends on what observable is being used for calibration
				double ordparam = 0;
				//wp: does average over all clones. for reduction is a sum over all parallel cores
				#ifdef _OPENMP
				#pragma omp parallel for reduction (+:ordparam)
				#endif
				for (int i=0; i<nc; i++) {
					ordparam += deformable[pushOffset + i]->getCurrent2Traj()/deformable[pushOffset + i]->getTimeTraj(); //average over traj
				}
				ordparam/= nc; 

				//wp: now use this value to invert back and determine which amplitulde of sync to use from the spline cal curve
				double new_amphi=0;
				new_amphi=spline_inv(ordparam);
				
				if (new_amphi <= cal_x0 ){ //wp: enforcing boundaries given calibration curve
					new_amphi=cal_x0;
				}
				else if( new_amphi	>= cal_xmax){
					new_amphi=cal_xmax;
				} 

				//std::cout << "ordparam clone average " << ordparam << " new amphi " << new_amphi << std::endl;
				//write calibration of the parameter to file
				if ( g_calibration_log.getOutputFile() !=""){ //wp: i.e. a name is specified 
					if(save_calibration){
						g_calibration_log.open();
						g_calibration_log.write_txt<double>(new_amphi); //wp: write to plain text
						g_calibration_log.close(); 
					}
				} 

				//wp: propagate to all clones; Note that this is done to the new copy of clones 
				#ifdef _OPENMP
				#pragma omp parallel for
				#endif
				for (int i=0; i<nc; i++)  (deformable[pullOffset + i]->getParametersPointer())->setAmplitude_phi(new_amphi); 
			
				} //end if 
			} //wp: control function

		); //wp: closes clones.doCloning(... 

		printf("\tFinished cloning. Computing final vals and saving.\n");
	
		//std::cout<< "time traj final " << (clones.finalSystem(0))->getTimeTraj() << std::endl;
		//std::cout<< "curr2 final " <<  (clones.finalSystem(0))->getCurrent2Traj() << std::endl;

		int num_output=2;
		clones.outputOP.assign(num_output, 0.);
		for (int i=0; i < nc; i++) {
			clones.outputOP[0] += ((clones.finalSystem(i))->getCurrentTraj())/((clones.finalSystem(i))->getTimeTraj());
			clones.outputOP[1] += ((clones.finalSystem(i))->getCurrent2Traj())/((clones.finalSystem(i))->getTimeTraj());
		}

		//wp: divides by # of clones
		for (unsigned int j=0;j<num_output;j++) { clones.outputOP[j] /= nc; }

		std::cout << std::endl;
		std::cout << "##bias parameter s: "    << sValue << std::endl
		          << "#SCGF "  << clones.outputPsi/N << std::endl
		          << "#<current> "     << clones.outputOP[0] << std::endl
		          << "#<current^2> "     << clones.outputOP[1] << std::endl
		          << "#am_phi final "     << ((clones.finalSystem(0))->getParametersPointer())->getAmplitude_phi() << std::endl
		          << "##time " << clones.outputWalltime << std::endl; 

		// output to file //wp: this is later read by cloning.py from binary file
		output.open();
		output.write<double>(clones.outputPsi);
		output.write<double>(clones.outputOP[0]); //wp: current =1/N/omega \sum_i dot(phi_i)
		output.write<double>(clones.outputOP[1]); //wp: current =1/N/omega^2 \sum_i dot(phi_i)^2
		output.write<double>(clones.outputWalltime);
		output.close(); 
	} //wp: end of 'nRuns' loop 

} //wp: end of main
