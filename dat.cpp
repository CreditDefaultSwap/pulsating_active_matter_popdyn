#include <string>
#include <iostream>
#include <vector>
#include <algorithm> 
#include "dat.hpp"

/*****************
 * DAT deformable *
 *****************/

// CONSTRUCTORS

Dat_def::Dat_def(std::string filename, bool loadWork) :
  numberParticles(), density(), systemSize(),
    randomSeed(), timeStep(), framesWork(0), dumpParticles(),
    dumpPeriod(),
  input(filename) {

  // HEADER INFORMATION
  input.read<const int>(&numberParticles);
  input.read<const double>(&density);
  input.read<const double>(&systemSize);
  input.read<const int>(&randomSeed);
  input.read<const double>(&timeStep);
  input.read<const bool>(&dumpParticles);
  input.read<const int>(&dumpPeriod); 
  input.read<const int>(&framesWork); //wp # of frames in which observable is computed

	//printf("\t Dat_Def After header information\n");

  // FILE PARTS LENGTHS
  headerLength = input.tellg();
  particleLength = 3*sizeof(double)*dumpParticles;
  frameLength = numberParticles*particleLength;
  //wp:
  workLength = 1*sizeof(double);

	//printf("\t Dat_Def After lengths\n");
  // ESTIMATION OF NUMBER OF COMPUTED WORK AND ORDER PARAMETER SUMS AND FRAMES
 	//wp: numberWork is the # of times work and order parameters were saved on file, which depends on 'framesWork*dumpPeriod' frequency
	//wp: framesWork = # of frames before work/order params are saved on file; for deformable total_time/dumpPeriod e.g. 1/0.1 = at the 10th frame it is saved 
	numberWork =(input.getFileSize() - headerLength - frameLength)/( framesWork*frameLength + workLength);
	//frames = after removing all the work/order params and header, then # of frames is size/particle frame lenght
  	frames = !dumpParticles ? 0 : 
    (input.getFileSize() - headerLength - numberWork*workLength)/frameLength;

  	numberWork =(input.getFileSize() - headerLength - frameLength)/( framesWork*frameLength + workLength);

  // FILE CORRUPTION CHECK
  if ( input.getFileSize() !=
    headerLength + frames*frameLength + numberWork*workLength ) {
    std::cerr << "Invalid file size." << std::endl;
    exit(1);
  }

  //Saved variables 
  if ( loadWork ) {
    double work; //wp: 'work' is dummy for saved observable
	int test;
    for (int i=0; i < numberWork; i++) {
			
	  //wp: Reads *after* whatever # input it is given
      input.read<double>(&work,
        headerLength                     // header
        + frameLength                    // frame with index 0 //wp initial save frame
        + (1 + i)*framesWork*frameLength // all following packs of framesWork frames //wp: all frame packs used while computing work
        + i*workLength);                 // previous values of the active work//wp: actual position where the quantities are dumped
		//printf("\tDat::\t current: %e\n", work);
      	currents.push_back(work);

    }
  }
}

// DESTRUCTORS

Dat_def::~Dat_def() {}

// METHODS

int Dat_def::getNumberParticles() const { return numberParticles; }
double Dat_def::getDensity() const { return density; }
double Dat_def::getSystemSize() const { return systemSize; }
int Dat_def::getRandomSeed() const { return randomSeed; }
double Dat_def::getTimeStep() const { return timeStep; }
int Dat_def::getFramesWork() const { return framesWork; }

long int Dat_def::getNumberWork() const { return numberWork; }
long int Dat_def::getFrames() const { return frames; }

std::vector<double> Dat_def::getCurrents() { return currents; }


double Dat_def::getPosition(
  int const& frame, int const& particle, int const& dimension) {
  // Returns position of a given particle at a given frame.

  return input.read<double>(
    headerLength                                     // header
    + frame*frameLength                              // other frames
    + particle*particleLength                        // other particles
    + (std::max(frame - 1, 0)/framesWork)*workLength // active work sums (taking into account the frame with index 0)
    + dimension*sizeof(double));                     // dimension
}

double Dat_def::getOrientation(int const& frame, int const& particle){
  // Returns position of a given particle at a given frame.

  return input.read<double>(
    headerLength                                     // header
    + frame*frameLength                              // other frames
    + particle*particleLength                        // other particles
    + (std::max(frame - 1, 0)/framesWork)*workLength // active work sums (taking into account the frame with index 0)
    + 2*sizeof(double));                             // positions //wp:i.e. after the x,y position
}

//end of file vile mile tile rile wile pile bile exile..

