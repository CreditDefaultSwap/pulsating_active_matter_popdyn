#ifndef PARTICLE_HPP
#define PARTICLE_HPP

#include <vector>
#include <string>
#include <random>

#include "maths.hpp"
#include "readwrite.hpp"
//wp: for sites
#include <deque>
/////////////
// CLASSES //
/////////////

//wp:for deformable case
class Parameters_def;
class Particle_def; 
class site;
class Deformables;

/*-------------relating deformable particles-------------*/
class Parameters_def {

  public:

    // CONSTRUCTORS
    Parameters_def();
	//Parameters_def(int Lx, int Ly, double rho, double omega, 
	//			double mu_r,  double mu_phi,  
	//			double amplitude_r, double amplitude_phi, 
	//			double radius,  double lambda, 
	//			double dt, double temperature );
	Parameters_def(double Lx, double Ly, double rho, double omega, 
				double mu_r,  double mu_phi,  
				double amplitude_r, double amplitude_phi, 
				double radius,  double lambda, 
				double dt, double temperature );

	//wp:copy from class object
    Parameters_def(Parameters_def* parameters_def);

    // METHODS
    double getLx() const;
    double getLy() const;
    int getNx() const;
    int getNy() const;
    double getRho() const;
    int getNumberParticles() const;
    double getOmega() const;
    double getMu_r() const; // returns mu_r of particles 
    double getMu_phi() const; // etc 
    double getAmplitude_r() const; 
    double getAmplitude_phi(); 
	void setAmplitude_phi(double am_phi_in); 
    double getRadius() const; 
    double getLambda() const; 
    double getDt() const; 
    double getTemperature() const; 

  private:

    // ATTRIBUTES

	double const Lx; double const Ly; //#box size; also # of cells if radius=1
	int const Nx; int const Ny; //# of cells based on radius; mutable = can be set inside
    double const rho; //# density of system
	int const N; //# of particles
	double const mu_r;  double const mu_phi; //mobilities	
	double const amplitude_r; 
	double amplitude_phi; //strength of feedback force //non-const so can be changed
	//double const amplitude_phi; //strength of feedback force 
	double const omega;
	double const radius;  
	double const lambda; 
	double const dt; 
	double const temperature; 

};

class Particle_def {

  public:

    // CONSTRUCTORS

    Particle_def(double d = 1);
    Particle_def(double x, double y, double ang, double d = 1);
    Particle_def(double x, double y, double ang, double d, int i, int j);

    // METHODS 
    double* position(); // returns pointer to position
    double* orientation(); // returns pointer to orientation
    double* velocity(); // returns pointer to velocity

	double* dorientation(); //wp returns pointer to change in torque 
    double diameter() const; // returns pointer to diameter 
    double* force(); // returns pointer to force 
    double* torque(); // returns pointer to aligning torque 
    int* box_index(); // returns pointer to containing box index 

	//controlled
	double* sync_sum();  // returns pointer to quantities 
	double* feedback_sum();
	double* dsync_sum();
	double* ordparam_ij(); 
	void set_ordparam_ij(double new_ordparam_sumij);

  private:

    // ATTRIBUTES

    double r[2]; // position (2D)
    double theta; // orientation
    double v[2]; // velocity (2D)

    double dtheta; //wp net change in orientation

    double const sigma; // diameter

    double f[2]; // force exerted on particle (2D)
    double gamma; // aligning torque

	//Integers i and j in waves particle code
	int index[2]; // containing box index

	//wp: for controlled quantities
	double sync_a_sum; // sum_j sin(theta_j-theta_i) 
	double feedback_b_sum; // sum_j d_theta_i Uij 
	double dsync_da_sum; // sum_j d_theta_i Uij 
	double ordparam_sumij; ///sum_ij cos(theta_j-theta_i)
};

class site{

private:
  std::vector<int> particle_ids;
  
public:
  //wp: No constructors; only methods
  std::vector<int> get_particle_ids(){
    return particle_ids;
  }

  void remove_particle_id(int id){
    for (int k = 0; k < this->particle_ids.size(); ++k){
      if(particle_ids[k]==id){
    	particle_ids.erase(particle_ids.begin() + k);
      }
    }
  }

  //wp:'deque' is a type of array structure which allows to quickly insert or remove items at the beginning or end
  //wp: unlike pointers, 'deque' allow non-contigous memory storage
  //wp: .push_back() = insert; .pop...() = remove
  void copy_particle_ids_to_deque(std::deque<int> &d){
    for(auto id : this->particle_ids){
      d.push_back(id);
    }
  }

  void add_particle_id(int id){
    this->particle_ids.push_back(id);
  }

};

/*  
 * Deformable particles 
 *  ------
 * Parameter class with all relevant parameters for deformable particles.
 * Creates particle class containing all info about deformable particles.
 * Finally defines a Deformables class which is the object used in cloning and simulation (iteration)
 */ 

	/*Structure of the writen binary file
	*  Parameters are stored in a binary file with the following structure:
	*
	*  [HEADER (see Deformable constructor )]
	*  |  |
	*
	*  [INITIAL FRAME (see System::saveInitialState)] (all double)
	*  ||                    FRAME 0                   ||
	*  ||          PARTICLE 1       | ... | PARTICLE N ||
	*  ||   R   | ORIENTATION | ... |     |     ...    ||
	*  || X | Y |    theta    | ... |     |     ...    ||
	*
	*  [BODY (see Deformables::saveNewState)] (all double)
	*  ||                    FRAME 1 + i*period                 || ... || FRAME 1 + (i + framesWork - 1)*period |~
	*  ||           PARTICLE 1        |  	...    | PARTICLE N || ... ||                  ...                  |~
	*  ||   R   | ORIENTATION   | ... |     ...    ||   ...     ||                  	   ...                  |~
	*  || X | Y |    theta      | ... |     ...    ||   ...     ||                         ...                  |~
	*
	*  ~|                  || ...
	*  ~|                  || ...
	*  ~| dynamic quantity || ...
	*  ~|      current     || ...
	*/

class Deformables {
	public:

    // CONSTRUCTORS

    Deformables();
	Deformables(Parameters_def* parameters_def, int seed); 
	Deformables(Parameters_def* params_d, int seed, std::string filename, bool dump, double period); 
	Deformables(Deformables* dummy, int seed, std::string filename, bool dump, double period);
	Deformables(Deformables* dummy, int seed, std::string filename, bool dump, double period, double total_iterations); 
	Deformables(Deformables* dummy, int seed, bool dump, double period, double total_iterations); //wp:cloning without .tmp.dat

    // cloning constructors
    Deformables(Deformables* dummy,int seed,int tau,std::string filename = "") : //wp:uses constructor 4 above: (...,seed,filename,dump,period) 
	Deformables(dummy, seed, filename, filename != "", tau) 
	{ resetDump(); }

    // cloning constructors #2//total time
    Deformables(Deformables* dummy,int seed, int period, int total_iterations, std::string filename = "") : //wp:uses constructor 5 above: (...,seed,filename,dump,period) 
	Deformables(dummy, seed, filename, filename != "", period, (double) total_iterations) 
	{ resetDump(); }

	//loading from file constructor; leaves sites empty
	Deformables(Deformables* dummy, std::string filename, double period);
	Deformables(Deformables* dummy, double period); //not writing trajs -> no *.tmp.dat files
	
    // DESTRUCTORS


    // METHODS

    int getNumberParticles() const; //returns # of deformable particles 

	Parameters_def getParameters(); //parameter class object
	Parameters_def* getParametersPointer(); //parameter class object pointer

	Particle_def* getParticle(int const& index);//gets particle by index 
	std::vector<Particle_def> getParticles_temp();
	std::vector<Particle_def> getParticles(); //particle_def vector
	std::vector<Particle_def*> getParticlesPointers(); //vector of pointers
	std::vector<Particle_def>* getParticlesPointer();//pointer *of vector*
	std::vector<std::vector<site>> getSites(); //gets sites vector 
	std::vector<std::vector<site>>* getSitesPointer();

	void setSites(std::vector<std::vector<site>> site_ext); //gets sites vector 
	void reconstructSites(); //wp: will re-construstct sites given positions
	//wp: cloning related
	double get_id();
	void set_id(int id_in); 

	double getBiasingParameter(); 
	void setBiasingParameter(double s); 

	double getWork(); 
	double getCurrent(); 
	double getCurrent2(); 
	double getControlAve();

	double getCurrentTraj(); 
	double getCurrent2Traj(); 
	void setCurrentTraj(double current_in); 
	void setCurrent2Traj(double current2_in); 

	Random* getRandomGenerator();
	void 	setGenerator(std::mt19937 rndeng); 

	double getTimeStep() const; 
	double getTimeTraj(); 
	void setTimeTraj(double time_traj_in);

	double getDensity(); 
	double getSystemSize(); 

	bool getSaveFlag(); 
	void setSaveFlag(bool saveParticles); 
	//wp: used for writing
	void flushOutputFile(); 
	std::string getOutputFile(); 
	void setOutputFile(std::string filename); 

	double getDumpPeriod(); 
	int* getDump(); 
	void resetDump(); 
	void copyDump(Deformables* deformable); 

	void saveInitialState(); 

	void copyState(std::vector<Particle_def>& newParticles_def); 
	void copyState(Deformables* deformable); 
	void copyParticles_temp(std::vector<Particle_def>& particles_new);

	void saveNewState(double time_spanned); 
	void saveNewStateExtras_control_strat(std::vector<Particle_def>& particles_old, double time_spanned, double time_total, double dt_step, bool save_traj);

  private:

    // ATTRIBUTES
	Parameters_def param_def; //all paramaters of dynamics

    int const numberParticles; // number of particles 
    double nWork=0;
	
	int id_num; //which clone track it is placed at

    std::vector<Particle_def> particles_def; // vector of particles 
    std::vector<Particle_def*> particles_def_pointers; // vector of particles 

    std::vector<Particle_def> particles_temp; // temp vector of particles 

    //wp: has to be left undeclared, then initiated with constructor
	std::vector<std::vector<site>> sites;

	//wp: cloning related
    std::vector<long int> velocitiesDumps; // locations in output file to dump velocities 

    int const randomSeed; // random seed
    Random randomGenerator; // random number generator

    bool dumpParticles; // dump positions and orientations to output file 
    //bool const dumpParticles; // dump positions and orientations to output file 
    int dumpFrame; // number of frames dumped since last reset
    int save_counter=1; // counter to increase time_period after each save 
    int work_save_period; //integer period at which work/order params are saved and/or averaged out 
 
    double const dumpPeriod; // period of dumping of positions and orientations in number of frames
    double const total_time; // total # of iterations done per trajectory while doing the cloning algo 
	
	double time_period; //simulation time at which data is collected
	double time_total=0; //wp: total simulation time done by the 'iteration_...' function
	double time_traj=0; //wp: total simulation up to the nth call of iteration i.e. cumulative over entire trajectory
	
	double current_sum=0, current_ave=0; // dummy for summing observable of interest; observable average after full sim. iterations
	double current2_sum=0, current2_ave=0; //wp: ave = over single in-between cloning simulation  
	double current_traj=0, current2_traj=0; //wp: traj = over *entire* trajectory 
	
	//control
	double control_sum=0, control_ave=0; //wp: ave = for every cloning step calculations 

    Write output; // output class//wp: to write to file in binary apparently

    double biasingParameter; // biasing parameter [cloning algorithm] 
}; //wp: deformables

////////////////
// PROTOTYPES //
////////////////


///////////////
// FUNCTIONS //
/////////////// 
double getOrderParameter(std::vector<Particle_def>& particles_def); 
void get_control_a_b(std::vector<Particle_def> p_def, Parameters_def params_def, double& control_a2, double& control_ab, double& control_da); 
//double get_control_a_b(std::vector<Particle_def> particles_def, Parameters_def params_def);
int sign_copy(double x);

#endif //for ifndef at start of file
