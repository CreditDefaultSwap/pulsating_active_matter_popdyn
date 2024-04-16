#ifndef ITERATION_HPP
#define ITERATION_HPP

#include "particle.hpp"

//----Deformable particles 
int sign(double x); 
int sign_int(int x); 
//site related
int NEXT(int k, double L); //wp: if it's at box's edge, return 0; else k+1 
int PREV(int k, double L); //wp: if at origin, return image; else k-1 
int IDX(int k, double L, int n);


//=======Iterates=========
void iterate_deformable_WCA_control_strat(Deformables* deformable, int Niter, bool save_traj); 


//=======MoveParticles===========
double MoveParticles_particleloop(std::vector<Particle_def>& particles, std::vector<std::vector<site>>& sites, std::vector<double> &Noise_r, std::vector<double> &Noise_phi, int N, int Nx, int Ny, double Lx, double Ly, float radius, float lambda, float amplitude_r, float amplitude_phi, float mu_phi, float mu_r, float omega, double PI, double dt, double Time); 

//=========forces======== 
void add_force_fullsync(Particle_def* p1, Particle_def* p2, double radius, double lambda, double amplitude_r, double amplitude_phi, double Ly, double Lx);



#endif
