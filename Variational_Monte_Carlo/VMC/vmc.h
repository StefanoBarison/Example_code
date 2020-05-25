#ifndef VMC
#define VMC

#include "random.h"
#include <vector> 

//For random number generation
int seed[4];
Random rnd;

//Function to calculate quantum mechanical quantities
double squared_wave(double x,double mu, double sigma);
double E_loc(double x, double mu, double sigma);
double Potential(double x);
double wavefunction(double x , double mu, double sigma);


//Quantities for monte carlo methods

double xold,xnew; // Position on real axes
double best_mu, best_sigma; //Optimization parameters
double delta; //delta for the uniform transition probability
int nblk,nsteps,nsample,nopt; //blocks and steps for each calculation
double l_rate; //learning rate for gradient descent


//Function for monte carlo methods

double Integrate(int nsample,double delta,double mu ,double sigma);
//void Optimize(double nopt,double l_rate,double guess_mu, double guess_sigma);
void Simulated_annealing(double nopt, double temp,double guess_mu, double guess_sigma);
void Input();
std::vector<double> Metropolis_sampling(double nsample,double delta, double x0,double mu, double sigma);

//Functions for data blocking
double error(std::vector<double> & av,std::vector<double> & av2, int n);
std::vector<double> Blocking(int nblocks,std::vector<double> & v);
void DataBlocking(std::vector <double> value, std::vector <double> &sum_prog,std::vector <double> &err_prog);


















#endif