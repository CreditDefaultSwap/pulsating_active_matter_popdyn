// Adapted cloning algorithm from RLJ.  //Further Adapted from code by YEK

#ifndef CLONINGSERIAL_HPP
#define CLONINGSERIAL_HPP

#include <random>
#include <vector>
#include <chrono>
#include <experimental/filesystem>

// note this is called CloningSerial even if it also supports openMP
// all loops over clones should be parallelised via openMP if supported
#ifdef _OPENMP
#include <omp.h>
#endif

#include "dat.hpp"
#include "particle.hpp"
#include "iteration.hpp"
#include "readwrite.hpp"

/////////////
// CLASSES //
/////////////


/*  CLONINGSERIAL
 *  -------------
 *  Cloning algorithm.
 */
//wp: This is the 'cloning' class object. Takes a generic system class to do clones of.
template<class SystemClass> class CloningSerial {
  /*  Contains all the parameters relevant to the cloning algorithm and the
   *  clones themselves.
   */

  public:

    // CONSTRUCTORS

    // lightweight constructor
    CloningSerial(
      int nClones, int nWork, double s, int method = 2,
      std::string cDir = "", std::string lFile = "", std::string sFile = "") :
      nc (nClones), cloneMethod (method), tau (nWork), sValue(s),
      clonesDirectory(cDir), loadInput(lFile), saveOutput(sFile) {

      systems.resize(2*nc);

      // save parameters to file
      saveOutput.write<int>(nc);
      saveOutput.write<int>(cloneMethod);
      saveOutput.write<int>(tau);
      saveOutput.write<double>(sValue); 
    }

	// wp: specifying which clones to save
    CloningSerial(
      int nClones, int nClones_save, int nWork, double s, int method = 2,
      std::string cDir = "", std::string lFile = "", std::string sFile = "") :
      nc (nClones), nc_save (nClones_save), cloneMethod (method), tau (nWork), sValue(s),
      clonesDirectory(cDir), loadInput(lFile), saveOutput(sFile) {

      systems.resize(2*nc);

      // save parameters to file
      saveOutput.write<int>(nc);
      saveOutput.write<int>(cloneMethod);
      saveOutput.write<int>(tau);
      saveOutput.write<double>(sValue);
    }

	// wp: specifying which clones to save
    CloningSerial(
      int nClones, int nClones_i,int nClones_f,int nClones_step, int nWork, double s, int method = 2,
      std::string cDir = "", std::string lFile = "", std::string sFile = "") :
      nc (nClones), nc_i(nClones_i), nc_f(nClones_f), nc_step(nClones_step), cloneMethod (method), tau (nWork), sValue(s),
      clonesDirectory(cDir), loadInput(lFile), saveOutput(sFile) {

      systems.resize(2*nc);

      // save parameters to file //wp: may chooose to avoid writing them 
      saveOutput.write<int>(nc);
      saveOutput.write<int>(nc_i);
      saveOutput.write<int>(nc_f);
      saveOutput.write<int>(nc_step);
      saveOutput.write<int>(cloneMethod);
      saveOutput.write<int>(tau);
      saveOutput.write<double>(sValue);
    }

	// wp: specifying which clones to save + from which time to start svaing from
    CloningSerial(
      int nClones, int nClones_i,int nClones_f,int nClones_step, int niter_save, int nWork, double s, int method = 2,
      std::string cDir = "", std::string lFile = "", std::string sFile = "") :
      nc (nClones), nc_i(nClones_i), nc_f(nClones_f), nc_step(nClones_step), cloneMethod (method), Niter_save(niter_save), tau (nWork), sValue(s),
      clonesDirectory(cDir), loadInput(lFile), saveOutput(sFile) {

      systems.resize(2*nc);

      // save parameters to file
      saveOutput.write<int>(nc);
      saveOutput.write<int>(nc_i);
      saveOutput.write<int>(nc_f);
      saveOutput.write<int>(nc_step);
      saveOutput.write<int>(cloneMethod);
      saveOutput.write<int>(tau);
      saveOutput.write<double>(sValue);
    }

	// wp: specifying clones to save(not implemented) + time to start saving from + if to save trajectories (bool)
    CloningSerial(
      int nClones, int nClones_i,int nClones_f,int nClones_step, int niter_save, int niter_save_clone,  
      bool save_traj, int nWork, double s, int method = 2,
      std::string cDir = "", std::string lFile = "", std::string sFile = "") :
      nc (nClones), nc_i(nClones_i), nc_f(nClones_f), nc_step(nClones_step), cloneMethod (method), 
	  Niter_save(niter_save), Niter_save_clone(niter_save_clone), save_particles(save_traj), 
	  tau (nWork), sValue(s),
      clonesDirectory(cDir), loadInput(lFile), saveOutput(sFile) {

      systems.resize(2*nc);

      // save parameters to file
      saveOutput.write<int>(nc);
      saveOutput.write<int>(nc_i);
      saveOutput.write<int>(nc_f);
      saveOutput.write<int>(nc_step);
      saveOutput.write<int>(cloneMethod);
      saveOutput.write<int>(tau);
      saveOutput.write<double>(sValue);

    }

    // DESTRUCTORS 
    // simple destructor, we delete the systems but the vector deals with itself
    ~CloningSerial() { deleteClones(); }

    // METHODS

	//wp: this creates ######.tmp.dat files starting from 000000 for trajectory reconstruction later
    //wp: will only do so within the 'clonesDirectory' if defined
    std::string cloneFilename(int i) {
      return clonesDirectory == "" ? "" : std::experimental::filesystem::path(
        std::experimental::filesystem::path(clonesDirectory) /
        [](int index)
          { return std::string(6 - std::to_string(index).length(), '0')
            + std::to_string(index) + std::string(".tmp.dat"); }
          (i)
        ).u8string();
    }

    void loadState(SystemClass* system); //wp: load cloning configurations from input file
    void saveState(); // save cloning configurations to output file

    void outSyncRandomGenerator() {
      // just for extra safety we can opt to throw our random number generators out of sync here
      // (the risk is that in a very big population, several clones might have the same seed,
      //  the effect here is to increase the effective range of seeds by safetyFactor
      #if 1
      for (int i=0;i<nc;i++) {
        int dum = 0;
        int safetyFactor = 1000;
        //int k = i % nc;
        if ( i==0 )
          std::cout << "#seed safetyFactor " << safetyFactor << std::endl;
        for (int r=0; r < i % safetyFactor; r++ ) {
          dum += (systems[i]->getRandomGenerator())->random01();
          dum += (systems[i+nc]->getRandomGenerator())->random01();
        }
      }
      #endif
    }

    void init(SystemClass* dummy, int masterSeed); // initialise list of clones

    template<typename F, typename G, typename H> void doCloning(
      double tmax, int initSim,
      F iterate, G getSWeight, H control); // this runs the cloning for total time tmax
    void writeTrajFiles(Write& clonesLog);

    void selectClones(int newClones[], double key[], int pullOffset);

    int binsearch(double key[], double val, int keylength); // this is a binary search used in clone selection

    void deleteClones() { for (int i=0; i<2*nc; i++) delete systems[i]; } // delete clones

    SystemClass* finalSystem(int index) { return systems[finalOffset + index]; }

    // ATTRIBUTES

    std::vector<SystemClass*> systems; // these are the clones, the vector has size (2.nc)
    int const nc; // how many clones
    int nc_save=nc, nc_i, nc_f, nc_step; //wp range of clones to save from i to f in steps
	int Niter_save, Niter_save_clone; //wp: #iteration at which to start saving cloning trajectories
    int const cloneMethod; // this determines which clone selection method to use (default 2?)

    int const tau; // number of simulation steps between cloning steps

    double sValue; // biasing parameter

    int iter = 0; // number of iterations performed
    int iter_save; //wp: # of iteration at which to start saving trajs, and obvs 
    bool save_traj=false, save_particles, save_calibration=false, save_clonelog=false, save_clonelog_time=false; //wp: flag to save trajectories, including observables, and separately the clones.log  file
    int runIndex = 0; // index of run

    std::string const clonesDirectory; // if != "": directory where trajectories are saved

    Read loadInput; // if file != "": input class from where initial configurations are loaded
    Write saveOutput; // if file != "": write class to where final configurations at each call of doCloning are saved

    int arrswitch = 0; // set of clones of considered

    double outputPsi; // this is the estimate of Psi
    double outputPsiOffset[2] {0, 0}; // used to consider values of psi from previous simulation via loadState
    std::vector<double> outputOP; // averages over the trajectories of the different clones of (0) active work (1) force part of the active work (2) orientation part of the active work (3) order parameter
    double outputWalltime; // time taken
    int finalOffset; // used to access the final population at the end of the run

    std::mt19937 cloneTwister; // random numbers

};

//wp: methods of 'CloningSerial' class
template<class SystemClass> void CloningSerial<SystemClass>:: 
	//loadState() {}
	loadState(SystemClass* dummy) {}
  // specialisation to be written in individual cloning *.cpp files

template<class SystemClass> void CloningSerial<SystemClass>:: saveState() {}
  // specialisation to be written in individual cloning *.cpp files

template<class SystemClass> void CloningSerial<SystemClass>::
  init(SystemClass* dummy, int masterSeed) {
  // initialise systems array with 2.nc copies of the "dummy" input
  // and gives a new seed to the random number generator on each copy
  // (also makes sure to clean up any old systems that are already in the array)
  // ... but note this does not initialise the state of the actual systems

  // std::cout << "#init CloningSerial" << std::endl;

  // delete any existing clones
  deleteClones();

  cloneTwister.seed(masterSeed);
  // std::uniform_int_distribution<int> sdist(0, 10000000);  // 7 digit random number
  std::uniform_int_distribution<int> sdist(0, 0x7fffffff);  // 31-bit integer, avoids any ambiguity with signed/unsigned

  int processSeeds[2*nc];  // BIG array of seeds(!), generate them here in order to ensure reproducibility
  for (int i=0;i<2*nc;i++) processSeeds[i] = sdist(cloneTwister);

  if (save_particles){
  #ifdef _OPENMP
  #pragma omp parallel for
  #endif
	  for (int i=0;i<2*nc;i++) {
		//wp: using constructor 5 directly: (Deformables* dummy, int seed, std::string filename, bool dump, double period, double total_iterations) 
		systems[i] = new SystemClass(dummy, processSeeds[i], cloneFilename(i),save_particles,(double) tau, tau); // create new system from copy of dummy, with random seed from processSeeds, computing active work and order parameter for every tau iterations 
	  }
  }
  else{//only saving clonelog or nothing
	//printf("\t\t clones.init: no trajs or clone.log save. Before 'systems' creation\n");
  #ifdef _OPENMP
  #pragma omp parallel for
  #endif
	  for (int i=0;i<2*nc;i++) {
		//constructure without filename (Deformables* dummy, int seed, bool dump, double period, double total_iterations)
		systems[i] = new SystemClass(dummy, processSeeds[i], save_particles,(double) tau, tau); // create new system from copy of dummy, with random seed from processSeeds, computing active work and order parameter for every tau iterations 
	  }
	}
	
	//printf("\t\t clones.init: after initializing systems pointers\n"); 
  outSyncRandomGenerator();
  for (int i=0;i<2*nc;i++) systems[i]->saveInitialState(); // important if trajectories are actually saved //wp: must save 2*nc if writing files 

}

template<class SystemClass> template<typename F, typename G, typename H>
  void CloningSerial<SystemClass>::
  doCloning( double tmax, int initSim, F iterate, G getSWeight, H control) {
  // this runs the cloning for total time tmax
  // input functions:
  // void iterate(SystemClass* system, int Niter): iterate system for Niter iterations //wp: here is where the simulation handle would come in
  // double getSWeight(SystemClass* system): returns product of biasing parameter and trajectory weight over last cloning step
  // void control(std::vector<SystemClass*> systems, int pullOffset, int pushOffset): modifications to controlled dynamics at the end of each cloning step

  // wp:!!!!!!!!!! --MAIN CLONING-- !!!!!!!!!!
  // !! this is the main cloning algorithm
  // wp:!!!!!!!!!! --MAIN CLONING-- !!!!!!!!!!

  //wp: if clone trajectory set -> write particles trajectory
	save_clonelog=!(clonesDirectory ==""); //wp: returns boolean from internal test 
	if (!save_clonelog && save_particles){ 
		printf("\t Looks like trajectories are being saved, but not the clone file (no save folder given) \n"); 
		printf( "\t Nothing saved. Check input. Exiting \n");
		exit(1);	
	}
								

  // clones log
  Write clonesLog(clonesDirectory == "" ? "" :
    std::experimental::filesystem::path
      (std::experimental::filesystem::path(clonesDirectory) / "clones.log") 
    .u8string());

  // slightly fancy c++ clock
  std::chrono::system_clock::time_point startTime =
    std::chrono::system_clock::now();

	//printf("clones.doCloning: tmax %e initSim %i\n", tmax, initSim);
	//printf("\tinitial val t_traj %e c2_traj %e \n", systems[0]->getTimeTraj(), systems[0]->getCurrent2Traj()); 
  // random initial condition for the nc systems that form the current population
  if ( loadInput.getInputFile() == "" ) {
  	printf("\t Doing initial relaxation, with iterations %i \n", initSim);
    arrswitch = 0; // is otherwise set in loadState
    #ifdef _OPENMP
    #pragma omp parallel for
    #endif
	//wp: This loops seems to set every clone to 0 bias, then run for tau steps(times relaxation steps)
    for (int i=0;i<nc;i++) {
      systems[i]->setBiasingParameter(0); // setting 0 biasing parameter (unmodified dynamics to sample initial configurations)
	  //wp:Actual simulation call! tau*initSim sets NITER 
      //iterate(systems[i], tau*initSim); // simulate an elementary number of steps
      iterate(systems[i], tau*initSim, save_traj); // simulate an elementary number of steps
      systems[i]->resetDump(); // reset dumps: important between different runs and to only count the relevant quantities within the cloning framework 
    }
  }
  // else load configurations from file //wp: here done in the .cpp file 

  // biasing parameter
  #ifdef _OPENMP
  #pragma omp parallel for
  #endif
  for (int i=0; i<2*nc; i++){ systems[i]->setBiasingParameter(sValue); systems[i]->set_id(i); }  // setting desired biasing parameter

  double lnX = 0.0, lnX_before=0., lnX_after=0.;  // this is used in our final estimate of psi
  bool lnX_set=false, pull_push_copy=false, pull_push_finished=false, clone_loop_finished=false, control_finished=false;

  //========================================
  //wp: Variables for sig catching procedure
  //wp: !!! --Given to global pointers!!!
  lnX_ptr = &lnX; //wp: !!!!this links the pointer of lnX to a **global** pointer used when sigterm handling
  lnX_before_ptr=&lnX_before;
  lnX_after_ptr=&lnX_after;

  lnX_set_ptr=&lnX_set;
  pull_push_copy_ptr=&pull_push_copy;
  pull_push_finished_ptr=&pull_push_finished;
  clone_loop_finished_ptr=&clone_loop_finished;
  control_finished_ptr=&control_finished;
  //========================================

  std::vector<double> sWeight;
  sWeight.resize(nc);

  double iter_max= tmax / (tau*systems[0]->getTimeStep());
  iter_max_ptr=&iter_max; //wp:!!!global pointer

  for (iter = 0; iter < iter_max; iter++) // for each iteration of the algorithm
  {
	lnX_set=false; //wp: keep track if lnX updated	
    pull_push_finished=false; 
	clone_loop_finished=false;
	control_finished=false;

    std::vector<double> upsilon(nc);  // these are the cloning factors

    // pushOffset refers to the current population, pullOffset will be the new population
    // (the systems array is of size 2.nc, only one half is "active" at each time)
    int pushOffset =  arrswitch   *nc;
    int pullOffset = (1-arrswitch)*nc;

	//wp: done as two separate iter counters so that clones.log can be saved in time independently of the rest
	if (iter == Niter_save_clone){ 
			printf("\t -- saving cloneslog at Niter.. %i\n", iter); 
			save_clonelog_time=save_clonelog*true;
	}

	if (iter == Niter_save){ 
			printf("\t -- saving at Niter.. %i\n", iter); 
			save_traj=save_particles*true;
			save_calibration=true;
	} //wp: only need to do it once


    #ifdef _OPENMP
    #pragma omp parallel for
    #endif
    for (int i = 0; i < nc; i++) // for each clone in the current population
    {
		//wp: No initial burn in since copying from last clone (NITER=tau)
        iterate(systems[pushOffset+i], tau, save_traj); // run dynamics

        sWeight[i] = getSWeight(systems[pushOffset+i]);
        upsilon[i] = exp(-sWeight[i]);
    }


	//wp: 'key' = cumulative sum of epsilons over all clones
    // construct key based on the upsilon params
    std::vector<double> key(nc+1);
    key[0] = 0.0;

	//wp: does the sum of upsilons at nc steps
    for (int i = 0; i < nc; i++) {
      key[i + 1] = key[i] + upsilon[i];
    }
    double totups = key[nc]; //Extract total of upsilons

    //Calculate cloning factor and store as log to avoid very large values
	//wp: totups/nc = Y = 1/nc sum_i Y_i, i \in clones
    double X = double(totups) / double(nc); //wp: this is like the mean val
    
	lnX_before=lnX; //wp: keep track if cloning loop completed or not
    lnX = lnX + log(X);
	lnX_after=lnX; lnX_set=true; //wp: lnx has been updated

    // decide who to copy (from which old clone does each new clone "pull" its state)
    int newClones[nc]; //wp: array of indexes denoting which clone to clone
    selectClones(newClones,key.data(),pullOffset);

	//=======!!!global pointers=======
  	pushOffset_ptr=&pushOffset;
  	pullOffset_ptr=&pullOffset;
  	newClones_ptr=newClones; //array 
	//================================


	//wp: clonesLog = log of which clones gets to be cloned //set from outside
	if (save_clonelog_time) for (int i=0;i<nc;i++) clonesLog.write<int>(newClones[i]);

    // this is the actual cloning step
	//wp: will copy in proportion to what index of clone is referred in 'newClones' e.g. newClones can have repeated values of index indicating that clone is getting more replicas
	
    pull_push_copy=true; 
    #ifdef _OPENMP
    #pragma omp parallel for
    #endif
    for (int i=0; i<nc; i++) {
	  //wp: copy from old clone index "pushOffset" to new clone index "pullOffset"
      systems[pullOffset + i]->copyState(systems[ pushOffset + newClones[i] ]); // clone particles
      systems[pullOffset + i]->copyDump(systems[ pushOffset + newClones[i] ]); // clone dumps
    }
    pull_push_copy=false; //wp:finished copying 
    pull_push_finished=true; 

    // CONTROLLED DYNAMICS
    control(systems, pullOffset, pushOffset, save_calibration);
	control_finished=true;

    arrswitch = 1 - arrswitch; //Set the other set of systems to be used in next time loop //wp: apparently, arrswitch just flips between 0 and 1 as to shift over in the 'systems' array. Example:
	clone_loop_finished=true;
//arrswitch=0 -> |old clones | new clones| -> arrswitch=1
//arrswitch=1 -> | new clones | old clones| -> arrswitch=0 rinse,repeat etc
  } //wp: tmax loop

  finalOffset = arrswitch * nc;

  //wp: lnX = sum_\beta log( Y_\beta) i.e. sum over all time chunks up to tmax
  //wp: last iter=tmax;

  //wp: since dt variable, use actual, average value over all clones as an approximation 
  double timeTraj_ave=0; 
  for (int i=0; i<nc; i++) timeTraj_ave+=systems[finalOffset+i]->getTimeTraj(); 
  timeTraj_ave/=nc; timeTraj_ave-=outputPsiOffset[0]; //wp: only want time spanned after loading, not entire run time
  outputPsi = (double(lnX) + outputPsiOffset[0]*outputPsiOffset[1])/(timeTraj_ave + outputPsiOffset[0]);

  // C++ clocks again
  std::chrono::system_clock::time_point endTime =
    std::chrono::system_clock::now();

  // this horrible syntax just computes elapsed time in seconds, at microsecond resolution
  outputWalltime = 1e-6*std::chrono::duration_cast<std::chrono::microseconds>
    (endTime - startTime).count();
 
  // WRITE TRAJECTORY FILES
  if (clonesDirectory != "" && save_traj) writeTrajFiles(clonesLog);

  // SAVE FINAL CONFIGURATIONS
  if ( saveOutput.getOutputFile() != "" ) saveState();

  // RUN INDEX
  runIndex++; // this line must be at the END of the function
} //wp: do clone method

template<class SystemClass> void CloningSerial<SystemClass>::
  writeTrajFiles(Write& clonesLog) {}
  // specialisation to be written in individual cloning *.cpp files
  //wp: done via the Write class in main

//wp:selects clones based on weight performance 'upsilon'
//wp: apparently returns the index of clones to be copied
template<class SystemClass> void CloningSerial<SystemClass>::
  selectClones(int newClones[], double key[], int pullOffset) {

  std::uniform_real_distribution<double> randc(0,1.0);  // random double in [0,1]

  // this offset determines who is pulling (it is either 0 or nc)
  double totups = key[nc];

  if ( cloneMethod == 1 ) {
    // this is the iid method
    #ifdef _OPENMP
    #pragma omp parallel for
    #endif
    for(int i=0; i<nc; i++) {
      // RLJ: we use the system rng of the destination system to decide where to "pull" from
      double rr = (systems[pullOffset+i]->getRandomGenerator())->random01()
        *totups;
      //cout << rr << " " << totups << " " << offset+i << endl;
      newClones[i] = binsearch(key, rr, nc+1);
    }
  }
  else if ( cloneMethod == 2 ) {
    // this is the eq method (should be default) //wp: equal spacing of upsilon; will ensure new copies conform to approximate proportion for each upsilon value of the clone. In other words, bigger values of upsilon, should in principle engender larger clones up to what is proportion of upsilon_i/upsilon_total 
	//wp: random01 = [0,1] uniform random #same for all clones 
    double alpha = (systems[pullOffset]->getRandomGenerator())->random01()
      *totups/nc;

    #ifdef _OPENMP
    #pragma omp parallel for
    #endif
    for(int i=0; i<nc; i++) {
      double rr = alpha + (i*totups/(double)nc);
      newClones[i] = binsearch(key, rr, nc+1);
    }
  }
  // if not method 1 or 2 then just copy everyone equally (should never happen)
  else { for(int i=0; i<nc; i++) newClones[i] = i; }
}

template<class SystemClass> int CloningSerial<SystemClass>::
  binsearch(double *key, double val, int keylength) {
  // this is a binary search used in clone selection

  int l = 0; /*Left hand limit*/
  int r = keylength - 1; /*Right hand limit*/
  int m; /*Midpoint*/
  int k; /*Element containing value*/
  bool elementfound = false; /*True when the element containing the value has been found*/
  while (elementfound == false)
  {
      m = int(floor(float(l + r) / 2)); /*Calculate midpoint*/
      if (val<key[m]) /*If value lower than midpoint, shift the right hand limit*/
      {
          r = m;
      }
      else /*Otherwise shift the left hand limit*/
      {
          l = m;
      }
      if (l == r - 1) /*Value lies in element between limits*/
      {
          elementfound = true;
      }
  }
  k = r - 1; /*Element index is index of upper limit when disregarding the index of the 0 value at the beginning of the key*/
  return k;
}

#endif
