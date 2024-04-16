#include <cmath>
#include <math.h>
#include <vector>

#include "iteration.hpp"
#include "particle.hpp"


/////////////////////////////////
// wp: Deformable Particles    //
/////////////////////////////////

// support functions 
// Compute the sign of an argument
int sign(double x){ return (x > 0) ? 1 : (-1); }
int sign_int(int x){ return (x > 0) ? 1 : (-1); }

// Choose next and previous elements with periodic boundary condition
int NEXT(int k, double L){ //wp: if it's at box's edge, return 0; else k+1
  return (k == (int) L-1) ? 0 : (k+1);
}

int PREV(int k, double L){ //wp: if at origin, return image; else k-1 
  return (k == 0) ? ((int) L-1) : (k-1);
}

int IDX(int k, double L, int n){
  return (k == (int) L+n) ? 0 : (k+n);
}

// Simulates system up to a certian time according to the dynamics of deformable particles with WCA potential
// Note that rather than doing fixed iterations it will do fixed time. The actual # of iterations varies due to adjustable time stepping 
void iterate_deformable_WCA_control_strat(Deformables* deformable, int Niter, bool save_traj) {
  	double PI=3.14159265358979323846;

  //Parameters
    //Get object address containing simulation parameters, constants etc
	Parameters_def params_def = deformable->getParameters(); 	
    //number of particles in system
	int N=params_def.getNumberParticles(); 
	//int Lx=params_def.getLx();
	//int Ly=params_def.getLy();
	double Lx=params_def.getLx();
	double Ly=params_def.getLy();
	//# of grids
	int Nx=params_def.getNx();
	int Ny=params_def.getNy(); 
    //particle props
	float radius = params_def.getRadius();
	float lambda = params_def.getLambda();
	float amplitude_r = params_def.getAmplitude_r();
	float amplitude_phi = params_def.getAmplitude_phi();
	float mu_phi = params_def.getMu_phi();
	float mu_r = params_def.getMu_r();
	float omega = params_def.getOmega();

    //wp: simulation related 
    //wp: adapting from original code: Niter = total time for these dynamics
    double time_spanned=0, time_spanned_new=0;
    double dt_d = params_def.getDt();
    double totaltime = Niter*dt_d;
    //float dt= (float) dt_d;
    double dt=dt_d;
    double temperature = params_def.getTemperature(), sqrt_mur_temp=sqrt(2*mu_r*temperature),sqrt_muphi_temp=sqrt(2*mu_phi*temperature);

	//wp: deformable particles and sites needed for simulation
	std::vector<Particle_def> particles_old;  //wp: copy of vectors 
	std::vector<Particle_def>* particles = deformable->getParticlesPointer();  //wp: pointer of vector -> will get modified 
	std::vector<std::vector<site>>* sites = deformable->getSitesPointer();  //wp: pointers of sites 
	Random* rand_gen=deformable->getRandomGenerator(); //wp: does the equivalent as the explicit declaration below but only is seeded once

	//wp:random generators 
	std::vector<double> Noise_r(2*N,0);
	std::vector<double> Noise_phi(N,0); 

	//wp: -------------------MAIN DYNAMICS LOOP------------------------ 
	// Run dynamics until the time iterator reaches the final time //wp: new random # on each call
	while(time_spanned<totaltime){ 
	  // Sample realizations of the noise 
	  for (int k=0;k<N;k++){//wp: rand num for each dimension
			Noise_r[2*k]   =sqrt_mur_temp*(rand_gen->gauss()); 
			Noise_r[2*k+1] =sqrt_mur_temp*(rand_gen->gauss()); 
			Noise_phi[k]   =sqrt_muphi_temp*(rand_gen->gauss()); 
	  }                                   

      //wp: copy old
      particles_old = deformable->getParticles(); 

	  // Move the particles according to the dynamics; // wp: spanned time will add to itself based on var. time step
	  time_spanned_new=MoveParticles_particleloop(*particles,*sites,Noise_r,Noise_phi,N,Nx,Ny,Lx,Ly,radius,lambda,amplitude_r,amplitude_phi,mu_phi,mu_r,omega,PI,dt,time_spanned); 

	 //save to file 
	 deformable->saveNewStateExtras_control_strat(particles_old,time_spanned_new,totaltime,time_spanned_new-time_spanned,save_traj);
	 time_spanned=time_spanned_new;
	 //to try to ensure totaltime is done exactly, rather than left over due to variable time step
	 if( time_spanned+dt >= totaltime){
	    dt=totaltime-time_spanned+1e-14; //to ensure an extra of float precision	
	 } 
	 
	}//end while 
}//end iterator 

//wp: main dynamic function for deformable particles 
double MoveParticles_particleloop(std::vector<Particle_def>& particles, std::vector<std::vector<site>>& sites, std::vector<double> &Noise_r, std::vector<double> &Noise_phi, 
int N, int Nx,int  Ny, double Lx, double Ly, float radius, float lambda, float amplitude_r, float amplitude_phi, float mu_phi, float mu_r, float omega, double PI, double dt, double Time){

	int k1, k2, k, i, j, i_new, j_new;
	double dqx, dqy, dphi; // Increment of coordinates 
	double test, max_amp=0; // Adaptative time stepping
	double dt_g, sdt_g; // Time step

	//particle objects 
	Particle_def* p1; //dummies for calc below 
	Particle_def* p2; //dummies for calc below 

	//resets forces and torques //wp and control
	for(k=0;k<N;k++){
		particles[k].force()[0]=0;
		particles[k].force()[1]=0;
		particles[k].torque()[0]=0; 

		particles[k].feedback_sum()[0]=0;
		particles[k].sync_sum()[0]=0;
		particles[k].dsync_sum()[0]=0;
	}   


	std::deque<int> temp_ids;

  //Double loop approach; can also use a cell list but for small particle size performance is actually better
  for(i=0;i<N;i++){ 
	p1 = & particles[i]; // Structure of particle 1
    for(j=i+1;j<N;j++){ // create a copy of particle ids to utilize newton's third law (correctly!) on same site                    
		p2 = & particles[j]; // Get its structure 
		add_force_fullsync(p1, p2, radius, lambda, amplitude_r, amplitude_phi, Ly, Lx); 
    } 
  } 

	// Compute maximum of interaction amplitude
	for(k=0;k<N;k++){
		test = mu_r*fabs(particles[k].force()[0])/radius;
		if(test>max_amp) max_amp = test;
		test = mu_r*fabs(particles[k].force()[1])/radius;
		if(test>max_amp) max_amp = test;
		test = (omega + mu_phi*particles[k].torque()[0])/(2*PI);
		if(test>max_amp) max_amp = test;
	} 

    // Adaptative time stepping
	if(max_amp * dt<1e-1){ 
		dt_g = dt; } 
    else{ 
		dt_g = 1e-1/max_amp; }

	sdt_g    = sqrt(dt_g);
	Time += dt_g;	


  //wp: Here is the actual time integration step.  // Update coordinates of all particles
  for(k=0;k<N;k++){ // Iterations over the N particles
    
    // Extract position and label of the box containing the particle
    i = (int) (particles[k].position()[0]/radius);
    j = (int) (particles[k].position()[1]/radius);

	// Dynamics 
    //wp: passing by reference so all structures can be accessed by '.' operations 
	dqx  = dt_g*mu_r*particles[k].force()[0] + sdt_g*Noise_r[2*k];
	dqy  = dt_g*mu_r*particles[k].force()[1] + sdt_g*Noise_r[2*k+1];
	dphi = dt_g*(omega + mu_phi*particles[k].torque()[0]) + sdt_g*Noise_phi[k];

    // New coordinates after one iteration
    particles[k].position()[0]  += dqx;
    particles[k].position()[1]  += dqy;
    particles[k].orientation()[0] += dphi;
	//wp: optional: save change in orientation
    //particles[k].dorientation()[0] = omega + mu_phi*particles[k].torque()[0] + Noise_phi[k];

    // Pay attention to the boundary conditions
    if(particles[k].position()[0]>Lx || particles[k].position()[0]<0)  particles[k].position()[0] -= Lx*floor(particles[k].position()[0]/Lx);
		if(particles[k].position()[1]>Ly || particles[k].position()[1]<0)  particles[k].position()[1] -= Ly*floor(particles[k].position()[1]/Ly);

   /* // New label of the box containing the particle
    i_new = (int) (particles[k].position()[0]/radius);
    j_new = (int) (particles[k].position()[1]/radius); 

    //printf("\t\t\t before update sites\n");
    // Update de Sites
    if(i!=i_new || j!=j_new){
     sites[i_new][j_new].add_particle_id(k);
     sites[i][j].remove_particle_id(k);
    }*/
  } //wp:k = N particle loop

  return Time;
} 


/*----------Computes deformable force between pair of particles----*/ 
// Compute the force between particles and add it to the force array + control quantities in the modified observable
// Done this way to optimize pair-based calls
void add_force_fullsync(Particle_def* p1, Particle_def* p2, double radius, double lambda, double amplitude_r, double amplitude_phi, double Ly, double Lx){

  double rad, rad6, rx, ry, r, r2, r6; // Interparticle distance
  double amp_rep, amp_phi, amp_align, dsync; // Interaction strengths
  
  // Interparticle distance
	rx= p2->position()[0] - p1->position()[0];
	ry= p2->position()[1] - p1->position()[1]; 
	rad = radius/(1+lambda)*(1 + lambda/2*(sin(p1->orientation()[0]) 
			+ sin(p2->orientation()[0]))); //wp: looks like a 2 was factored out from within so that radius= \sigma_0 

	// Pay attention to boundary conditions
	if(fabs(rx)>(Lx-rad)) rx -= sign(rx)*Lx; 
	if(fabs(ry)>(Ly-rad)) ry -= sign(ry)*Ly;

	r2 = rx*rx + ry*ry;
	r = sqrt(r2); 

	//wp: add full sync (i.e. does not depend on distance)
	amp_align = amplitude_phi*sin(p2->orientation()[0]-p1->orientation()[0]);		
	dsync = -cos(p2->orientation()[0]-p1->orientation()[0]); //amphi must be applied onSaveNew when saving averages

	p1->torque()[0] += amp_align; 
	p2->torque()[0] -= amp_align;

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
		amp_rep   = 12*amplitude_r*rad6/r6/r*(rad6/r6 - 1);
		amp_phi   = amp_rep*radius*r*lambda/(2*rad*(1+lambda));

		// Interaction on p1
		p1->force()[0] 		 -= amp_rep*rx/r;
		p1->force()[1] 		 -= amp_rep*ry/r;
		p1->torque()[0] += -amp_phi*cos(p1->orientation()[0]); //wp: align + feedback ( derivative gives 'am_phi*cos()')

		//control quantities 1
		p1->feedback_sum()[0] += -amp_phi*cos(p1->orientation()[0]); //wp: adds to b = \sum_j dUij/d\th_i 
		
		// Interaction on p2
		p2->force()[0] 		 += amp_rep*rx/r;
		p2->force()[1] 		 += amp_rep*ry/r;
		p2->torque()[0] += -amp_phi*cos(p2->orientation()[0]);
		
		//control quantities 2
		p2->feedback_sum()[0] += -amp_phi*cos(p2->orientation()[0]); //wp: adds to b = \sum_j dUij/d\th_i  // even to change in i -> j

	}
}

//end of file
