#include <cmath>
#include <math.h>
#include <vector>
#include <string>
#include <algorithm>
//wp: for sites
#include <deque>
//wp: .hpp  = header files
#include "dat.hpp"
#include "particle.hpp"
#include "maths.hpp"
//#include "iteration.hpp"
//wp: to get cpu id
#include <omp.h>

/////////////
// CLASSES //
/////////////

/***************************************************/
/***************Deformable class******************/
/***************************************************/

/******************
 * Parameters_def *
 ******************/
//wp: class_name::constructor():initiation list
//wp: the initiation list is the list of constants that get initiated for that class object; So below, they'd all get zero val
//wp------constructors
Parameters_def::Parameters_def() : 
  Lx(3), Ly(3), Nx((int) ceil(Lx)), Ny((int) ceil(Ly)), 
rho(1), N(1), omega(0), mu_r(0), mu_phi(0), amplitude_r(0),
    amplitude_phi(0), radius(1), lambda(0), dt(0.001), temperature(1)  {}


//constructor 2: Lx and Ly as doubles
Parameters_def::Parameters_def(double Lx, double Ly, double rho, 
				double omega, double mu_r,  double mu_phi,  
				double amplitude_r, double amplitude_phi, 
				double radius,  double lambda, double dt, double temperature) :
  				Lx(Lx), Ly(Ly), rho(rho), omega(omega), 
				N((int) round(Lx*Ly*rho)), 
				mu_r(mu_r), mu_phi(mu_phi), 
				amplitude_r(amplitude_r), amplitude_phi(amplitude_phi), 
				radius(radius), lambda(lambda), 
				Nx((int) ceil(Lx/radius)), Ny((int) ceil(Ly/radius)), //want to round-up because these will denote # of boxes. If Lx>int(Lx) then add the remant as another box hence need an extra index in sites 
				//Nx((int) floor(Lx/radius)), Ny((int) floor(Ly/radius)),  
				dt(dt), temperature(temperature) {}

//copy from pointer of another but same class object
Parameters_def::Parameters_def(Parameters_def* paramd) :
  				Lx(paramd->getLx()), Ly(paramd->getLy()), 
  				Nx(paramd->getNx()), Ny(paramd->getNy()), 
				N(paramd->getNumberParticles()), 
				rho(paramd->getRho()), omega(paramd->getOmega()), 
				mu_r(paramd->getMu_r()), mu_phi(paramd->getMu_phi()), 
				amplitude_r(paramd->getAmplitude_r()), 
				amplitude_phi(paramd->getAmplitude_phi()), 
				radius(paramd->getRadius()), lambda(paramd->getLambda()),
				dt(paramd->getDt()), temperature(paramd->getTemperature())
				 {} 


//wp-------methods
double Parameters_def::getLx() const { return Lx; } 
double Parameters_def::getLy() const { return Ly; } 
int Parameters_def::getNx() const { return Nx; } 
int Parameters_def::getNy() const { return Ny; } 
double Parameters_def::getRho() const { return rho; } 
double Parameters_def::getOmega() const { return omega; } 
int Parameters_def::getNumberParticles() const { return N; } 
double Parameters_def::getMu_r() const { return mu_r; } 
double Parameters_def::getMu_phi() const { return mu_phi; }
double Parameters_def::getAmplitude_r() const { return amplitude_r; }

double Parameters_def::getAmplitude_phi() { return amplitude_phi; }
void Parameters_def::setAmplitude_phi(double am_phi_in) { amplitude_phi=am_phi_in; }

double Parameters_def::getRadius() const { return radius; }
double Parameters_def::getLambda() const { return lambda; }
double Parameters_def::getDt() const { return dt; }
double Parameters_def::getTemperature() const { return temperature; }

/***************
* Particle def * 
 ***************/

// CONSTRUCTORS

Particle_def::Particle_def(double d) :
  r{0, 0}, theta(0), v{0, 0}, sigma(d), f{0, 0}, gamma (0) {}
Particle_def::Particle_def(double x, double y, double ang, double d) :
  r{x, y}, theta(ang), v{0, 0}, sigma(d), f{0, 0}, gamma (0) {}
Particle_def::Particle_def(double x, double y, double ang, double d, int i, int j) :
  r{x, y}, theta(ang), v{0, 0}, sigma(d), f{0, 0}, gamma (0), index{i, j} {}

// METHODS

double* Particle_def::position() { return &r[0]; } // returns pointer to position
double* Particle_def::orientation() { return &theta; } // returns pointer to orientation 
double* Particle_def::dorientation() { return &dtheta; } //wp returns pointer to change in torque

double* Particle_def::velocity() { return &v[0]; } // returns pointer to velocity


double Particle_def::diameter() const { return sigma; } // returns pointer to diameter


double* Particle_def::force() { return &f[0]; }; // returns pointer to force
double* Particle_def::torque() { return &gamma; } // returns pointer to aligning torque

int* Particle_def::box_index() { return &index[0]; } // returns pointer to i,j containing boxes 


double* Particle_def::sync_sum() { return &sync_a_sum; }; // returns pointer to force
double* Particle_def::feedback_sum() { return &feedback_b_sum; }; //
double* Particle_def::dsync_sum() { return &dsync_da_sum; }; //
double* Particle_def::ordparam_ij() { return &ordparam_sumij; }; //
void Particle_def::set_ordparam_ij(double new_ordparam_sumij){ordparam_sumij=new_ordparam_sumij;}

/****sites has built-in methods****/
///*nothing to put here; class already added on the header .hpp*////

/***************
 * Deformables *
 ***************/ 
Deformables::Deformables() :
	param_def(Parameters_def()), randomSeed(0), randomGenerator(), 
	numberParticles(param_def.getNumberParticles()), 
    particles_def(param_def.getNumberParticles()),  particles_temp(param_def.getNumberParticles()),
	output(""), dumpParticles(false), dumpPeriod(-1),dumpFrame(1),total_time(1),
	time_period(0)
	{}

//constructor from parameters object
Deformables::Deformables(Parameters_def* params_d, int seed) :
	param_def(params_d),randomSeed(seed),randomGenerator(seed),
	numberParticles(param_def.getNumberParticles()),
	particles_def(param_def.getNumberParticles()), particles_temp(param_def.getNumberParticles()),
	particles_def_pointers(), output(""), 
	dumpParticles(false), dumpPeriod(1),dumpFrame(1),total_time(1),time_period(param_def.getDt()),
	sites(param_def.getNx(),std::vector<site>(param_def.getNy())) {   
	//sites(param_def.getLx(),std::vector<site>(param_def.getLy())) {   
		int k = 0, l=0; // Iterator on particles
		int i=0, j=0; // Indices of the box
		double X=0, Y=0; // Sampled coordinates
		int N= getNumberParticles(); 
		double radius=param_def.getRadius();
		double Lx,Ly; Lx=param_def.getLx(); Ly=param_def.getLy();

		// Iterations over the deposited particles //wp: random pos and phase
		while(k<N){ 
			// Initialization //random in [0,1]
			X = randomGenerator.random01()*Lx;
			Y = randomGenerator.random01()*Ly;
			//printf("X,Y %e, %e\n", X,Y);

			// Extract label (i,j) of corresponding box
			i = (int) (X/radius);
			j = (int) (Y/radius); //wp: if wanted to change aspect ratio 'b', it would be (Y/(b*radius)) or viceversa on Lx direction 
			//printf("X,Y, i, j %e, %e, %i, %i\n", X,Y,i,j);

			// Store spatial coordinates
			particles_def[k].position()[0] = X;
			particles_def[k].position()[1] = Y;

			// Store lattice coordinates
			particles_def[k].box_index()[0] = i;
			particles_def[k].box_index()[1] = j;

			// Initial size variables
			particles_def[k].orientation()[0] = 2*M_PI*randomGenerator.random01();

			// Add particle in the system
			sites[i][j].add_particle_id(k);
			k++; 
		} //wp: while loop 
  } //wp:end constructor body

//constructor for cloning algo #1
Deformables::Deformables(Parameters_def* params_d, int seed, std::string filename, bool dump, double period) :
	param_def(params_d),randomSeed(seed),randomGenerator(seed),
	numberParticles(param_def.getNumberParticles()),
	particles_def(param_def.getNumberParticles()), particles_temp(param_def.getNumberParticles()),
	particles_def_pointers(), output(filename), 
	dumpParticles(dump), dumpPeriod(period),dumpFrame(1), total_time(period),time_period(period*param_def.getDt()), 
	sites(param_def.getNx(),std::vector<site>(param_def.getNy())) {   
		int k = 0, l=0; // Iterator on particles
		int i, j; // Indices of the box
		double X, Y; // Sampled coordinates
		X=0; Y=0;
		int N= getNumberParticles(); 
		double radius=param_def.getRadius();
		double Lx,Ly; Lx=param_def.getLx(); Ly=param_def.getLy();

		// Iterations over the deposited particles //wp: random pos and phase
		while(k<N){ 
			// Extract label (i,j) of corresponding box
			i = (int) (X/radius);
			j = (int) (Y/radius); //wp: if wanted to change aspect ratio 'b', it would be (Y/(b*radius)) or viceversa on Lx direction 

			// Store spatial coordinates
			particles_def[k].position()[0] = X;
			particles_def[k].position()[1] = Y;

			// Store lattice coordinates
			particles_def[k].box_index()[0] = i;
			particles_def[k].box_index()[1] = j;

			// Initial size variables
			particles_def[k].orientation()[0] = 2*M_PI*randomGenerator.random01();

			// Add particle in the system
			sites[i][j].add_particle_id(k);
			k++;

			X+=1.;
			if(X >= (double) Lx){ X=0.; Y+=1.;}

 			if(k > Lx*Ly-1){ 
			X=particles_def[l].position()[0]+0.2;
			Y=particles_def[l].position()[1];
			l+=1; 
			}
		}

  } //wp:end constructor body

//constructor for cloning algo #2
Deformables::Deformables(Deformables* dummy, int seed, std::string filename, bool dump, double period) :
	param_def(dummy->getParameters()),randomSeed(seed),randomGenerator(seed),
	numberParticles(param_def.getNumberParticles()),
	particles_def(param_def.getNumberParticles()), particles_temp(param_def.getNumberParticles()),
	particles_def_pointers(), output(filename), 
	dumpParticles(dump), dumpPeriod(period),dumpFrame(1),total_time(period), time_period(period*param_def.getDt()), 
	sites(param_def.getNx(),std::vector<site>(param_def.getNy())) {   
		//Body of constructor

		//*. Loop is to start the sites
		int k = 0, l=0; // Iterator on particles
		int i, j; // Indices of the box
		double X, Y; // Sampled coordinates
		X=0; Y=0;
		int N= getNumberParticles(); 
		double radius=param_def.getRadius();
		double Lx,Ly; Lx=param_def.getLx(); Ly=param_def.getLy();


		//printf("\t inside constructor seed: %i\n", seed);
		// Iterations over the deposited particles //wp: random pos and phase
		while(k<N){ 
			// Initialization //random in [0,1]
			X = randomGenerator.random01()*Lx;
			Y = randomGenerator.random01()*Ly;

			// Extract label (i,j) of corresponding box
			i = (int) (X/radius);
			j = (int) (Y/radius); //wp: if wanted to change aspect ratio 'b', it would be (Y/(b*radius)) or viceversa on Lx direction 

			// Store spatial coordinates
			particles_def[k].position()[0] = X;
			particles_def[k].position()[1] = Y;

			// Store lattice coordinates
			particles_def[k].box_index()[0] = i;
			particles_def[k].box_index()[1] = j;

			// Initial size variables
			particles_def[k].orientation()[0] = 2*M_PI*randomGenerator.random01();

			//	printf("X,Y orientation %e, %e, %e\n", X,Y, particles_def[k].orientation()[0]);

			// Add particle in the system
			sites[i][j].add_particle_id(k);
			k++; 
		} 
		
		//*. write header with system parameters to output file
		output.write<int>(getNumberParticles());
		output.write<double>(getDensity());
		output.write<double>(getSystemSize());
		output.write<int>(randomSeed);
		output.write<double>(getTimeStep());
		output.write<bool>(dumpParticles);
		output.write<int>(dumpPeriod);
		output.close();

		//wp: uncomment copyState below if want all particles to start the same
		// copy positions and orientations and update cell list
		//copyState(dummy);
		// copy dumps
		//copyDump(dummy);
} //wp:end constructor body

//constructor for cloning algo, allowing dumping and total time separately #3
Deformables::Deformables(Deformables* dummy, int seed, std::string filename, bool dump, double period, double total_iterations) :
	param_def(dummy->getParameters()),randomSeed(seed),randomGenerator(seed),
	numberParticles(param_def.getNumberParticles()),
	particles_def(param_def.getNumberParticles()), particles_temp(param_def.getNumberParticles()),

	particles_def_pointers(), output(filename), 
	dumpParticles(dump), dumpPeriod(period), dumpFrame(0), total_time(total_iterations), time_period(period*param_def.getDt()), 
	sites(param_def.getNx(),std::vector<site>(param_def.getNy())) {   
		//Body of constructor

		//*. Loop is to start the sites
		int k = 0, l=0; // Iterator on particles
		int i, j; // Indices of the box
		double X, Y; // Sampled coordinates
		X=0; Y=0;
		int N= getNumberParticles(); 
		double radius=param_def.getRadius();
		double Lx,Ly; Lx=param_def.getLx(); Ly=param_def.getLy();
		double dump_observable;

		//printf("\t inside constructor seed: %i\n", seed);
		// Iterations over the deposited particles //wp: random pos and phase
		while(k<N){ 
			// Initialization //random in [0,1]
			X = randomGenerator.random01()*Lx;
			Y = randomGenerator.random01()*Ly;

			// Extract label (i,j) of corresponding box
			i = (int) (X/radius);
			j = (int) (Y/radius); //wp: if wanted to change aspect ratio 'b', it would be (Y/(b*radius)) or viceversa on Lx direction 

			// Store spatial coordinates
			particles_def[k].position()[0] = X;
			particles_def[k].position()[1] = Y;

			// Store lattice coordinates
			particles_def[k].box_index()[0] = i;
			particles_def[k].box_index()[1] = j;

			// Initial size variables
			particles_def[k].orientation()[0] = 2*M_PI*randomGenerator.random01(); 

			// Add particle in the system
			sites[i][j].add_particle_id(k); 

			k++; 
		} 
		
		//*. Deteremines work/order param frequency of normalization and/or saving
		dump_observable = total_time*param_def.getDt()/time_period;

		//*. write header with system parameters to output file
		if(dumpParticles){ //wp: will write .tmp.dat files for when saving trajectories
			//Write output(filename); //wp:redeclares output with filename -> creates .tmp.dat in disk
			//output.setOutputFile(filename);
			output.open();
			output.write<int>(getNumberParticles());
			output.write<double>(getDensity());
			output.write<double>(getSystemSize());
			output.write<int>(randomSeed);
			output.write<double>(getTimeStep());
			output.write<bool>(dumpParticles);
			output.write<int>(dumpPeriod);
			output.write<int>((int) dump_observable);
			output.close();
		}

		//wp: uncomment copyState below if want all particles to start the same
		// copy positions and orientations and update cell list
		//copyState(dummy);
		// copy dumps
		//copyDump(dummy);
} //wp:end constructor body


//constructor for cloning #4: Do *not* instantiate 'Write output' so no .tmp.dat files are written (no trajs will be saved)
Deformables::Deformables(Deformables* dummy, int seed, bool dump, double period, double total_iterations) :
	param_def(dummy->getParameters()),randomSeed(seed),randomGenerator(seed),
	numberParticles(param_def.getNumberParticles()),
	particles_def(param_def.getNumberParticles()), particles_def_pointers(), particles_temp(param_def.getNumberParticles()), 
	dumpParticles(dump), dumpPeriod(period), dumpFrame(0), total_time(total_iterations), time_period(period*param_def.getDt()), 
	sites(param_def.getNx(),std::vector<site>(param_def.getNy())) {   
		//Body of constructor

		//*. Loop is to start the sites
		int k = 0, l=0; // Iterator on particles
		int i, j; // Indices of the box
		double X, Y; // Sampled coordinates
		X=0; Y=0;
		int N= getNumberParticles(); 
		double radius=param_def.getRadius();
		double Lx,Ly; Lx=param_def.getLx(); Ly=param_def.getLy();
		double dump_observable;

		// Iterations over the deposited particles //wp: random pos and phase
		while(k<N){ 
			// Initialization //random in [0,1]
			X = randomGenerator.random01()*Lx;
			Y = randomGenerator.random01()*Ly;

			// Extract label (i,j) of corresponding box
			i = (int) (X/radius);
			j = (int) (Y/radius); //wp: if wanted to change aspect ratio 'b', it would be (Y/(b*radius)) or viceversa on Lx direction 

			// Store spatial coordinates
			particles_def[k].position()[0] = X;
			particles_def[k].position()[1] = Y;

			// Store lattice coordinates
			particles_def[k].box_index()[0] = i;
			particles_def[k].box_index()[1] = j;

			// Initial size variables
			particles_def[k].orientation()[0] = 2*M_PI*randomGenerator.random01();

			//	printf("X,Y orientation %e, %e, %e\n", X,Y, particles_def[k].orientation()[0]);
			//if(k==N-1)	printf("i,j %i %i\n", i,j);

			// Add particle in the system
			sites[i][j].add_particle_id(k); 

			k++; 
		} 
		
		//*. Deteremines work/order param frequency of normalization and/or saving
		dump_observable = total_time*param_def.getDt()/time_period; 

} //wp:end constructor body


//constructor for cloning algo when loading file
Deformables::Deformables(Deformables* dummy, std::string filename, double period) :
	param_def(dummy->getParameters()),randomSeed(-1),randomGenerator(randomSeed),
	numberParticles(param_def.getNumberParticles()),
	particles_def(numberParticles),particles_temp(param_def.getNumberParticles()), 
	particles_def_pointers(), output(filename), 
	dumpParticles(1), dumpPeriod(period), dumpFrame(0), total_time(period), time_period(period*param_def.getDt()), 
	sites(param_def.getNx(),std::vector<site>(param_def.getNy())) {   
		//Body of constructor

		//*. Loop is to start the sites
		int k = 0, l=0; // Iterator on particles
		int i, j; // Indices of the box
		double X, Y; // Sampled coordinates
		X=0; Y=0;
		int N= getNumberParticles(); 
		double radius=param_def.getRadius();
		double Lx,Ly; Lx=param_def.getLx(); Ly=param_def.getLy();
		double dump_observable;

		// Iterations over the deposited particles //wp: random pos and phase
		while(k<N){  //wp: Leaves sites *empty* but with rand. positions and orientations
			// Initialization //random in [0,1]
			X = randomGenerator.random01()*Lx;
			Y = randomGenerator.random01()*Ly;

			// Store spatial coordinates
			particles_def[k].position()[0] = X;
			particles_def[k].position()[1] = Y;

			// Initial size variables
			particles_def[k].orientation()[0] = 2*M_PI*randomGenerator.random01(); 

			k++; 
		} 
		
		//*. Deteremines work/order param frequency of normalization and/or saving
		dump_observable = total_time*param_def.getDt()/time_period;

		//*. write header with system parameters to output file
		output.write<int>(getNumberParticles());
		output.write<double>(getDensity());
		output.write<double>(getSystemSize());
		output.write<int>(randomSeed);
		output.write<double>(getTimeStep());
		output.write<bool>(dumpParticles);
		output.write<int>(dumpPeriod);
		output.write<int>((int) dump_observable);
		output.close();
}

//constructor for cloning algo when loading file, but not saving traj (so no output(filename) construction) (no *.tmp.dat)
Deformables::Deformables(Deformables* dummy, double period) :
	param_def(dummy->getParameters()),randomSeed(-1),randomGenerator(randomSeed),
	numberParticles(param_def.getNumberParticles()),
	particles_def(numberParticles), particles_def_pointers(), particles_temp(param_def.getNumberParticles()), 
	dumpParticles(0), dumpPeriod(period), dumpFrame(0), total_time(period), time_period(period*param_def.getDt()), 
	sites(param_def.getNx(),std::vector<site>(param_def.getNy())) {   
	//sites(param_def.getLx(),std::vector<site>(param_def.getLy())) {   
		//Body of constructor

		//*. Loop is to start the sites
		int k = 0, l=0; // Iterator on particles
		int i, j; // Indices of the box
		double X, Y; // Sampled coordinates
		X=0; Y=0;
		int N= getNumberParticles(); 
		double radius=param_def.getRadius();
		double Lx,Ly; Lx=param_def.getLx(); Ly=param_def.getLy();
		double dump_observable;

		// Iterations over the deposited particles //wp: random pos and phase
		while(k<N){  //wp: Leaves sites *empty* but with rand. positions and orientations
			// Initialization //random in [0,1]
			X = randomGenerator.random01()*Lx;
			Y = randomGenerator.random01()*Ly; 

			// Store spatial coordinates
			particles_def[k].position()[0] = X;
			particles_def[k].position()[1] = Y; 

			// Initial size variables
			particles_def[k].orientation()[0] = 2*M_PI*randomGenerator.random01(); 

			k++; 
		} 
		
		//*. Deteremines work/order param frequency of normalization and/or saving
		dump_observable = total_time*param_def.getDt()/time_period; 
}

//--------------------
//------methods------- 
Particle_def* Deformables::getParticle(int const& index) {return &(particles_def[index]); }
int Deformables::getNumberParticles() const {return numberParticles;};
Parameters_def Deformables::getParameters() {return param_def;}; 
Parameters_def* Deformables::getParametersPointer() {return & param_def;}; 

std::vector<Particle_def> Deformables::getParticles(){return particles_def;};
std::vector<Particle_def> Deformables::getParticles_temp(){return particles_temp;};
//wp: Returns pointer *of vector*
std::vector<Particle_def>* Deformables::getParticlesPointer(){return & particles_def;};

std::vector<Particle_def*> Deformables::getParticlesPointers(){
	//Adds address of particle pointers to the pointer vector
	//done explicitly by counting as the 'auto' for loop method fails
	for (int k=0; k<numberParticles; k++){
		particles_def_pointers.push_back(& particles_def[k]);
	}
	return particles_def_pointers;
}; 

std::vector<std::vector<site>> Deformables::getSites(){return sites;};
//wp: Returns pointer of sites vector 
std::vector<std::vector<site>>* Deformables::getSitesPointer(){return & sites;};

void Deformables::setSites(std::vector<std::vector<site>> site_ext){ sites = site_ext; return;}; //sets sites vector
void Deformables::reconstructSites(){ //wp: will re-construstct sites given positions
		int k=0,i,j; 
		double X,Y, radius=param_def.getRadius();

		while(k<numberParticles){  
			// Based on stored spatial coordinates
			X=particles_def[k].position()[0];
			Y=particles_def[k].position()[1];

			// Extract label (i,j) of corresponding box
			i = (int) (X/radius);
			j = (int) (Y/radius); //wp: if wanted to change aspect ratio 'b', it would be (Y/(b*radius)) or viceversa on Lx direction 

			// Store lattice coordinates
			particles_def[k].box_index()[0] = i;
			particles_def[k].box_index()[1] = j; 

			// Add particle in the system
			sites[i][j].add_particle_id(k); 

			k++; 
		} 
	//wp: should be it
}; 

//wp: ---- cloning related/standardization 
double Deformables::get_id() { return id_num; }
void Deformables::set_id(int id_in) { id_num=id_in; }

double Deformables::getWork() { return nWork; }

//wp: current: 1/N/omega \sum_i d\theta_i/dt , i = particles
double Deformables::getCurrent() { return current_ave; }
double Deformables::getCurrent2() { return current2_ave; }
//wp: controlled variable = action difference with unbiased for within single cloning integration time
double Deformables::getControlAve() { return control_ave; }

double Deformables::getCurrentTraj() { return current_traj; }
double Deformables::getCurrent2Traj() { return current2_traj; }
void Deformables::setCurrentTraj(double current_in) {current_traj=current_in;}
void Deformables::setCurrent2Traj(double current2_in) {current2_traj=current2_in;}

double Deformables::getDumpPeriod() { return dumpPeriod; }

double Deformables::getBiasingParameter() { return biasingParameter; }
void Deformables::setBiasingParameter(double s) { biasingParameter = s; }

Random* Deformables::getRandomGenerator() { return &randomGenerator; }
void Deformables::setGenerator(std::mt19937 rndeng) { randomGenerator.setGenerator(rndeng); }

double Deformables::getDensity() { return param_def.getRho(); }
double Deformables::getSystemSize() { return param_def.getLx(); }

double Deformables::getTimeStep() const { return param_def.getDt(); }
double Deformables::getTimeTraj() { return time_traj; }
void Deformables::setTimeTraj(double time_traj_in) { time_traj=time_traj_in; }

void Deformables::flushOutputFile() { output.flush(); }
std::string Deformables::getOutputFile() {return output.getOutputFile();}
void Deformables::setOutputFile(std::string filename ) {output.setOutputFile(filename);} //wp: modifies file name

int* Deformables::getDump() { return &dumpFrame; }

bool Deformables::getSaveFlag() { return dumpParticles; } 
void Deformables::setSaveFlag(bool saveParticles) { dumpParticles=saveParticles; } 

void Deformables::resetDump() {
  // Reset time-extensive quantities over trajectory.  
  dumpFrame = 0;
  current_traj =0;
  current2_traj =0;
  time_traj =0;

  particles_temp[0].set_ordparam_ij(0); 
}

//wp: copies total traj quantities, including time_traj from other systems 
void Deformables::copyDump(Deformables* deformable) {
	dumpFrame = deformable->getDump()[0];
	current_traj = deformable->getCurrentTraj();
	current2_traj = deformable->getCurrent2Traj();
	time_traj = deformable->getTimeTraj(); 
}

void Deformables::saveInitialState() {
  // Saves initial state of particles to output file.

  // output
  if ( dumpParticles ) {
    output.open();

    for (int i=0; i < numberParticles; i++) { // output all particles
      // POSITIONS
      for (int dim=0; dim < 2; dim++) { // output position in each dimension
        output.write<double>(particles_def[i].position()[dim]);
      }
      // ORIENTATIONS
      output.write<double>(particles_def[i].orientation()[0]); // output orientation
    }

    output.close();
  }

  // reset dump
  resetDump();
}

void Deformables::copyState(std::vector<Particle_def>& newParticles_def) {
// Copy positions and orientations.

  for (int i=0; i < numberParticles; i++) {
    for (int dim=0; dim < 2; dim++) {
      // POSITIONS
      particles_def[i].position()[dim] = newParticles_def[i].position()[dim];
    }
    // ORIENTATIONS
    particles_def[i].orientation()[0] = newParticles_def[i].orientation()[0];
  } 
  
}

void Deformables::copyState(Deformables* deformable) {
  // Copy positions and orientations.

  //wp: copies *all* properties of particles, including control params etc
  particles_def=deformable->getParticles(); 

  // UPDATING CELL LIST //wp: equiv to using cell lists
  sites=deformable->getSites();

 //wp: for phase speeds
  particles_temp=deformable->getParticles_temp(); 
}

//wp:used in iteration to save to file
//wp: This is adapted to the deformable particles. It only *saves* to file 
// All particles are updated in place by default
void Deformables::saveNewState(double time_spanned){
	// Saves new state of particles //wp, das it 

	////////////
	// SAVING //
	////////////
	dumpFrame++; //wp--just count the # of invokes foraverages later I guess? 

	//wp: saves file every dumpPeriod*dt = time_spanned
	//printf("this is time_spanne and dumpPeriod*dumpFrame %e %e\n", time_spanned, time_period*dumpFrame);
	if ((dumpParticles && time_spanned > time_period )) {
		output.open(); 
		  //printf("Now printing particles to\n");
		  //std::cout<< output.getOutputFile() << std::endl;
	      for(int i=0; i<numberParticles; i++){
			 //printf("%e %e %e \n",particles_def[i].position()[0],particles_def[i].position()[1],particles_def[i].orientation()[0]);
        	  output.write<double>(particles_def[i].position()[0]);
        	  output.write<double>(particles_def[i].position()[1]);
			  // ORIENTATION
		      output.write<double>(particles_def[i].orientation()[0]); 
	      }
	  } 

	  // CLOSE FILE
	  if ( output.is_open() ) { output.close(); } 

}

void Deformables::copyParticles_temp(std::vector<Particle_def>& newParticles_def){
	//particles_temp=particles_new;	
	// Copy positions and orientations.  
	  for (int i=0; i < numberParticles; i++) {
		for (int dim=0; dim < 2; dim++) {
		  // POSITIONS
		  particles_temp[i].position()[dim] = newParticles_def[i].position()[dim]; }
		// ORIENTATIONS
		particles_temp[i].orientation()[0] = newParticles_def[i].orientation()[0];
	  }

}

//wp: Computes orderparam directly from particles.dsync()//requires dsync to not be multiplied by amphi in code
void Deformables::saveNewStateExtras_control_strat(std::vector<Particle_def>& particles_old, double time_spanned, double time_total, double dt_step, bool save_traj){
    // Saves new state of particles 
    double mu_phi=param_def.getMu_phi(), am_phi=param_def.getAmplitude_phi(), mu_temperature_4=0.25/param_def.getTemperature()/mu_phi, mu_phi_2=0.5*mu_phi;
    double control_a2=0, control_ab=0, control_da=0;
    double current2_sum_old=0, control_a2_old=0, control_ab_old=0, control_da_old=0, dt_seg=0; 
	double speed_sum=0; 

    ////////////
    // SAVING //
    ////////////

    if (dumpFrame == 0){ current2_ave=0; current_ave=0; control_sum=0; copyParticles_temp(particles_old); } //wp: resets at the beginning of the call 

    if (dt_step < 2e-5 && time_spanned < time_total ){ //time step too short i.e. an unphysicality in integration. Wait for next point
        //printf("dt step, dt_min %e %e\n", dt_step, dt_min);
        //!!particle configurations must remain untouched! and they are :^)
        for(int i=0; i<numberParticles; i++){
            particles_def[i].feedback_sum()[0]=particles_old[i].feedback_sum()[0]; //old will be kept as the new old in next iteration
            particles_def[i].sync_sum()[0]=particles_old[i].sync_sum()[0]; 
            particles_def[i].dsync_sum()[0]=particles_old[i].dsync_sum()[0]; 
			
			//store old phase configurations
			//particles_def_temp[i].orientation()[0]=particles_old[i].orientation()[0]
        }

        particles_def[1].set_ordparam_ij(particles_old[1].ordparam_ij()[0]+dt_step); //time lapse since last 'good' point 

        if (dumpFrame > 0) dumpFrame++;
        return;
    } 

    // ORDER PARAMETER computation //wp: uses current2_sum variables as dummy for convinience but actually order param =  1/N^2 \sum_i \sum_j cos(theta_i-theta_j), j=/ i 
    for(int i=0; i<numberParticles; i++){

        //old
        control_a2_old += particles_old[i].sync_sum()[0]*particles_old[i].sync_sum()[0]; //wp: pointer[0] = access value of that pointer
        control_ab_old += particles_old[i].feedback_sum()[0]*particles_old[i].sync_sum()[0];
        control_da_old += particles_old[i].dsync_sum()[0];

        //control quantities, summing over i (j sum already implicitly done at the dynamics step)
        control_a2 += particles_def[i].sync_sum()[0]*particles_def[i].sync_sum()[0]; //wp: pointer[0] = access value of that pointer
        control_ab += particles_def[i].feedback_sum()[0]*particles_def[i].sync_sum()[0]; 
        control_da += particles_def[i].dsync_sum()[0]; 
      } 

    //wp: This approach is basically a trapezoidal integration of observables for accuracy. However does require careful storing of t-1 and t+1 quantities.  
    //	  particles(t) store control info about particles(t-1) (due to how it is implemented during the dynamics integration step). 
    //    //    Hence calculating quants for particles(t) correspond to particles(t+1); Technically it belongs to the clone, but adding to the first particle
    //    for convinience
    //wp: saves order parameter which corresponds to the proper time step
	current2_sum_old=-control_da;
    dt_seg=particles_old[1].ordparam_ij()[0]; //time step corresponding to the t-1 segment
    particles_def[0].set_ordparam_ij(current2_sum_old);
    particles_def[1].set_ordparam_ij(dt_step); //stores time step since variable and want to correspond to right interval
	

    if(dumpFrame==0){dt_seg=0;} //ignore t=-1 contribution due to the way this calc is behind t-1 in time 

    //wp: given variable time step, must multiply for every dt_step //using stratonovich convention (var(t-1)+v(t))/2
    current2_ave +=0.5*(sqrt(numberParticles+particles_old[0].ordparam_ij()[0])+sqrt(numberParticles+particles_def[0].ordparam_ij()[0]))*dt_seg; 

    //control var: int_t \sum_i 1/4T(a^2+2*a*b) + 1/2d/d\theta_i a
    control_da_old*=am_phi; control_da*=am_phi; //since not done at the force level
    control_sum += 0.5*( ((control_a2_old+2.*control_ab_old)*mu_temperature_4+mu_phi_2*control_da_old)
                        + (control_a2+2.*control_ab)*mu_temperature_4        +mu_phi_2*control_da)*dt_seg;

	//reset sums for next dynamic step
    current2_sum=0; 

    //wp: saves file every save_counter*dumpPeriod*dt = time_spanned
    if(save_traj){ //wp: to control from which cloning step time it is saved
        if (( dumpParticles && time_spanned > time_total )) { //wp: save only once per iteration call, regardless of # of iterations
            output.open();
            for(int i=0; i<numberParticles; i++){
                output.write<double>(particles_def[i].position()[0]);
                output.write<double>(particles_def[i].position()[1]);
                // ORIENTATION
                output.write<double>(particles_def[i].orientation()[0]);
            }
            //close file
            output.close();
        }
    }

    dumpFrame++; //wp--just counts the # of calls
    //Done full simulation, compute average + save if wanted
    if ( time_spanned >= time_total ){
        //---Last point requires computation of particle_def quantities; save current particle_def quants to old which is what it originally corresponds to
        //---*does not* apply to phase speed calc because avg assumed to carry over a single long traj, rather than sectioned as here
		particles_old[0].set_ordparam_ij(particles_def[0].ordparam_ij()[0]); 
		control_a2_old=control_a2; control_ab_old=control_ab; control_da_old=control_da;
		//--then get the actual control quantities for this last step
		get_control_a_b(particles_def, param_def, control_a2, control_ab, control_da );  //pass by reference; so modified in place

		if( log(abs(control_ab/control_ab_old)) < 1.098612289 ){ //if new repulsion term less than a factor of exp(1.0986)~3  than the old
			particles_def[0].set_ordparam_ij(-control_da); 

			current2_ave +=0.5*(sqrt(numberParticles+particles_old[0].ordparam_ij()[0])+sqrt(numberParticles+particles_def[0].ordparam_ij()[0]))*dt_step; //wp: corrected for full expression in yiwei's work

			control_da*=am_phi;
			control_sum += 0.5*( ((control_a2_old+2.*control_ab_old)*mu_temperature_4+mu_phi_2*control_da_old)
							+ (control_a2+2.*control_ab)*mu_temperature_4        +mu_phi_2*control_da)*dt_step;
        }
        else{//approximate as only last step 
			control_sum += ((control_a2_old+2.*control_ab_old)*mu_temperature_4+mu_phi_2*control_da_old)*dt_step;
            current2_ave+=sqrt(numberParticles+particles_old[0].ordparam_ij()[0])*dt_step;
        } 
        //----------
        
		//=====phase speed calc 
		for (int i=0; i<numberParticles; i++) speed_sum +=(particles_def[i].orientation()[0]-particles_temp[i].orientation()[0]); 

		//speed: (v_old+v_new)/2; since done every time_total=constant time steps
		current_ave=0.5*(particles_temp[0].ordparam_ij()[0]+speed_sum)/numberParticles; 

		//current_ave=(speed_sum)/numberParticles;
		particles_temp[0].set_ordparam_ij(speed_sum); //to be passed on for later calcs

        //wp:obvs for the entire trajectory of cloning
        current_traj += current_ave;  //wp: to be divided by the total time at the end 
        current2_traj += current2_ave/numberParticles;
        time_traj += time_spanned; //wp: keep track of total actual time simulated since dt are variable 

		current_ave = current_ave/time_spanned;
        current2_ave = current2_ave/numberParticles/time_spanned;
        control_ave = control_sum/numberParticles/time_spanned; //wp: in principle divided by N and T, but it will later be multiplied by sNT so cancels out
        if(save_traj){
            if (dumpParticles){ output.open(); output.write<double>(current2_ave); output.close(); }
        }

        dumpFrame=0;
    }

}


///////////////
// FUNCTIONS //
/////////////// 

//computes a=sync b=feedback da=related to ordparam //doesn't pass arg by arg;
void get_control_a_b(std::vector<Particle_def> p_def, Parameters_def params_def, double& control_a2, double& control_ab, double& control_da){ 
  	double rad, rad6, rx, ry, r, r2, r6; // Interparticle distance
 	double amp_rep, amp_phi, amp_align, dsync; // Interaction strengths
	//double radius,  lambda,  amplitude_r,  amplitude_phi,  Ly,  Lx;

	int N = p_def.size();
    double Lx=params_def.getLx();
    double Ly=params_def.getLy();
    //particle props
    double radius = params_def.getRadius();
    double lambda = params_def.getLambda();
    double amplitude_r = params_def.getAmplitude_r();
    double amplitude_phi = params_def.getAmplitude_phi();

	for( int i=0; i < N; i++){
	    p_def[i].feedback_sum()[0]=0;
        p_def[i].sync_sum()[0]=0;
        p_def[i].dsync_sum()[0]=0;
	}

	for (int i=0; i < N; i++) { // loop over particles
		Particle_def* p1=&p_def[i];
 		for (int j=i+1; j < N; j++) { // loop over particles 
			Particle_def* p2=&p_def[j];

			// Interparticle distance
			rx= p2->position()[0] - p1->position()[0];
			ry= p2->position()[1] - p1->position()[1];
			rad = radius/(1+lambda)*(1 + lambda/2*(sin(p1->orientation()[0]) + sin(p2->orientation()[0]))); //wp: looks like a 2 was factored out from within so that radius= \sigma_0

			// Pay attention to boundary conditions
			if(fabs(rx)>(Lx-rad)) rx -= sign_copy(rx)*Lx;
			if(fabs(ry)>(Ly-rad)) ry -= sign_copy(ry)*Ly;

			r2 = rx*rx + ry*ry;
			r = sqrt(r2);

			//wp: add full sync (i.e. does not depend on distance)
			amp_align = amplitude_phi*sin(p2->orientation()[0]-p1->orientation()[0]);
			dsync = -cos(p2->orientation()[0]-p1->orientation()[0]);  //amphi must be added *outside*
			//printf("\tcontrol: amp_align %e, dsync %e\n", amp_align, dsync);

			//control quantities 1
			p1->sync_sum()[0] += amp_align; //wp: adds to a = \sum_j sin(\th_j - \th_i)
			p2->sync_sum()[0] -= amp_align; //wp: adds to a = \sum_j sin(\th_j - \th_i) // i-> j flips sign on sin

			p1->dsync_sum()[0] += dsync ; //wp: adds to da/dtheta = -\sum_j cos(theta_j-theta_i)
			p2->dsync_sum()[0] += dsync; //wp: adds to da/dtheta = -\sum_j cos(theta_j-theta_i)		

			// Add contribution to the force
			if(r<rad){
				// Interaction strengths 
 				r6        = r2*r2*r2;
		        rad6      = rad*rad*rad*rad*rad*rad;

				// Interaction strengths
				amp_rep   = 12*amplitude_r*rad6/r6/r*(rad6/r6 - 1);
				amp_phi   = amp_rep*radius*r*lambda/(2*rad*(1+lambda));

				//control quantities 1
				p1->feedback_sum()[0] += -amp_phi*cos(p1->orientation()[0]); //wp: adds to b = \sum_j dUij/d\th_i

				//control quantities 2
				p2->feedback_sum()[0] += -amp_phi*cos(p2->orientation()[0]); //wp: adds to b = \sum_j dUij/d\th_i  // even to change in i -> j

				//printf("\tcontrol:r<rad:amps: rad6 %e r6 %e r %e \n", rad6,r6,r); 
				//printf("\tcontrol:r<rad:amps: amp_rep %e amp_phi %e \n", amp_rep, amp_phi);
				//printf("\tcontrol:r<rad:force: a_i %e , b_i %e\n", p1->sync_sum()[0], p1->feedback_sum()[0]);
			} 
		} //j loop
	} //i loop 

	//calculate final control quantities sum_i
	//double control_a2=0, control_ab=0, control_da=0;
	control_a2=0, control_ab=0, control_da=0; //since passed by ref, ensure that it is cleared before adding up
	for (int i=0; i<N; i++){
        //control quantities, summing over i (j sum already implicitly done at the dynamics step)
        control_a2 += p_def[i].sync_sum()[0]*p_def[i].sync_sum()[0]; //wp: pointer[0] = access value of that pointer
        control_ab += p_def[i].feedback_sum()[0]*p_def[i].sync_sum()[0];
        control_da += p_def[i].dsync_sum()[0];
        //printf("control_a_b: a_i: %e, b_i: %e da: %e \n", p_def[i].sync_sum()[0], p_def[i].feedback_sum()[0], p_def[i].dsync_sum()[0]);
        //printf("control_a_b: i %i, a_i: %e, b_i: %e da: %e \n",i, p_def[i].sync_sum()[0], p_def[i].feedback_sum()[0], p_def[i].dsync_sum()[0]);
	} 
    //printf("control_a_b: a_i: %e, b_i: %e da: %e \n", control_a2, control_ab, control_da); 
	return;
}

double getOrderParameter(std::vector<Particle_def>& particles_def) {
  //wp: Returns order parameter = 1/N | \sum_j e^j\theta | averaged all the particles.
  double order_partial=0, order_param=0; 

  int N = particles_def.size();
  for (int i=0; i < N; i++) { // loop over particles
  	for (int j=i+1; j < N; j++) { // loop over particles
      order_partial += cos(particles_def[i].orientation()[0]-particles_def[j].orientation()[0]);
    }
  }
  
  order_param=1./N*sqrt(N+2*order_partial);

  return order_param;
}

// Compute the sign of an argument
int sign_copy(double x){ return (x > 0) ? 1 : (-1); }
