#ifndef DAT_HPP
#define DAT_HPP

#include <string>
#include <vector>

#include "readwrite.hpp"

/////////////
// CLASSES //
/////////////

class Dat_def;

/*  DAT_def
 *  ----
 *  Read files as defined by the Deformables class (see particle.hpp).
 */

class Dat_def {

  public:

    // CONSTRUCTORS

    Dat_def(std::string filename, bool loadOrder = true);

    // DESTRUCTORS

    ~Dat_def();

    // METHODS

    int getNumberParticles() const; // returns number of particles
    double getDensity() const; // returns packing fraction
    double getSystemSize() const; // returns system size
    int getRandomSeed() const; // returns random seed
    double getTimeStep() const; // returns time step
    int getFramesWork() const; // returns number of frames on which to sum the active work before dumping

    long int getNumberWork() const; // returns number of computed work sums
    long int getFrames() const; // returns number of frames

    std::vector<double> getCurrents(); // returns vector of average current 

    double getPosition(
      int const& frame, int const& particle, int const& dimension);
      // Returns position of a given particle at a given frame.
    double getOrientation(int const& frame, int const& particle);
      // Returns position of a given particle at a given frame.

    void close() { input.close(); } // close file stream
    void open() { input.open(); } // open file stream



  private:

    // ATTRIBUTES

    int const numberParticles; // number of rotors
    double const density; //wp: number density
    double const systemSize; //wp: number density
	

    double const timeStep; // time step
    int const dumpPeriod; // period of dumping of orientations in number of frames
    bool const dumpParticles; // positions and orientations dumped in file
    int const randomSeed; // random seed

    Read input; // input class

    long int numberWork; // number of computed work sums
    long int frames; // number of frames
    long int headerLength; // length of header in input file
    long int particleLength; // length the data of a single particle takes in a frame
    long int workLength; // length the data of a single work dump takes in a file 

    long int frameLength; // length the data of a single frame takes in a file
    long int orderLength; // length the data of a single order dump takes in a file

    long int numberOrder; // number of computed order parameter sums
	int framesWork; //wp: # of frames for work related outputs

	//wp: observable related
    std::vector<double> currents; // computed active work sums

};



#endif
