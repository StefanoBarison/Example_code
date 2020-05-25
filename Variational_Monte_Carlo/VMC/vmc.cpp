#include "vmc.h"
#include "random.h"
#include <iostream>
#include <fstream>
#include <cmath>
#include <vector>



using namespace std;

int main(){

	// As first thing, we initialize the quantities useful for our calculations
	Input();

	//Then we want to find a delta for our transition such that the acceptance rate is almost 0.5
	// So let's try a metropolis calculation using the input delta

	//double mu=0.764856,sigma=0.594263;

	//Best values
	//mu = 0.774863
	//sigma = 0.631326

	double mu=1.0, sigma=0.5;

	vector<double> trial=Metropolis_sampling(nsample,delta,0.0,mu,sigma);

	//If we have an acceptance rate of 0.50 we can continue

	//Exercise 08.1

	//with our trials sampling we can show that we are able to calculate by blocking method 
	//the expectation value of the hamiltonian over a trial wave function

	
	vector<double> e_trial(nsample,0);

	for(int i=0;i<nsample;i++){
		e_trial[i]=E_loc(trial[i],mu,sigma);
	}

	//Now do the blocking method
	vector<double> ave_trial=Blocking(nblk,e_trial);
	vector<double> sum_trial(nblk,0);
	vector<double> err_trial(nblk,0);

	DataBlocking(ave_trial,sum_trial,err_trial);

	//Then print the output

	ofstream outdata;
	outdata.open("example.dat");

	for(int i=0;i<nblk;i++){
		outdata<< i+1<<" "<<sum_trial[i]<<" "<<err_trial[i]<<endl;
	}

	outdata.close();

	
	//Exercise 08.2

	/*
	----------------------------------------------------------------------------------------------------------------
	// Now we want to minimize the expectation value by optimizing E_opt over the R^2 space of parameters (mu,sigma)
	// to do so we do a variational monte carlo, where metropolis sampling is used to calculate the integrals
	// and a gradient descent algorithm to optimize


	//first make an initial guess over the parameters

	best_mu=1;
	best_sigma=0.5;

	//...establish the parameters of the optimization
	nopt=1000;
	l_rate=0.01;

	//And now optimize

	Optimize(nopt,l_rate,best_mu, best_sigma);

	//Unfortunately, due to the nature of energy function, the gradient descent method is unstable
	----------------------------------------------------------------------------------------------------------------
	*/

	//Let's try with simulated annealing:
	//first make an initial guess over the parameters

	best_mu=1;
	best_sigma=0.5;

	//...establish the parameters of the optimization
	nopt=1000;
	double temp=10; //temperature of the annealing
	Simulated_annealing(nopt,temp,best_mu,best_sigma);

	cout<<"After the optimization the value of the parameters found are:"<<endl;
	cout<<"mu = "<<best_mu<<endl;
	cout<<"sigma = "<<best_sigma<<endl;
	cout<<"--------------------------------------------------------------"<<endl;


	//Now we can repeat the calculation of the exercise 08.1, but with new parameters and new deltas

	cout<<"Evaluation with optimal parameters"<<endl;

	cout<<"First 10000 steps useful to equilibrate"<<endl;

	double new_delta=1.1;

	vector<double> equilibrate=Metropolis_sampling(10000,new_delta,0.0,best_mu,best_sigma);

	cout<<"Then run the equilibrated algorithm"<<endl;

	vector<double> final=Metropolis_sampling(nsample,new_delta,equilibrate[9999],best_mu,best_sigma);


	vector<double> e_final(nsample,0);

	for(int i=0;i<nsample;i++){
		e_final[i]=E_loc(final[i],best_mu,best_sigma);
	}

	//Now do the blocking method
	vector<double> ave_final=Blocking(nblk,e_final);
	vector<double> sum_final(nblk,0);
	vector<double> err_final(nblk,0);

	DataBlocking(ave_final,sum_final,err_final);

	//Then print the output

	outdata.open("variational_energy_results.dat");

	for(int i=0;i<nblk;i++){
		outdata<< i+1<<" "<<sum_final[i]<<" "<<err_final[i]<<endl;
	}

	outdata.close();

	//In the end we can also save the data sampled, in order to see the sampled squared wavefunction
	//Let's create the histogram
	cout<<"Creating the histogram"<<endl;

	double range=6;
	int nbins=1000;
	double binsize=range/nbins;

	vector<double> hist(nbins,0);

	for(int j=0;j<nsample;j++){
		for(int k=0;k<nbins;k++){
			if(final[j]>(-3+k*binsize) && final[j]<(-3+(k+1)*binsize)){
				hist[k]=hist[k]+1.0;
			}
		}
	}

	outdata.open("Wavefunction.dat");

	for(int i=0;i<nbins;i++){
		outdata<<-3+i*binsize<<" "<<hist[i]<<endl;
	}

	outdata.close();
	
	return 0;
}


double wavefunction(double x , double mu, double sigma){
	return exp(-pow(x-mu,2)/(2*pow(sigma,2)))+exp(-pow(x+mu,2)/(2*pow(sigma,2)));
}
double squared_wave(double x,double mu, double sigma){
	return exp(-pow(x-mu,2)/pow(sigma,2))+exp(-pow(x+mu,2)/pow(sigma,2))+2*exp(-(pow(x,2)+pow(mu,2))/pow(sigma,2));
}

double E_loc(double x, double mu, double sigma){

	//double kinetic = -(1/(2*pow(sigma,2)))*(1-pow(x-mu,2)/pow(sigma,2))*exp(-pow(x-mu,2)/pow(sigma,2))-(1/(2*pow(sigma,2)))*(1-pow(x+mu,2)/pow(sigma,2))*exp(-pow(x+mu,2)/pow(sigma,2));
	double kinetic = -0.5*(exp(-pow(x+mu,2)/(2*pow(sigma,2)))*((1+exp(2*mu*x/pow(sigma,2)))*(pow(mu,2)-pow(sigma,2)+pow(x,2))-2*mu*x*(-1+exp(2*mu*x/pow(sigma,2))))/pow(sigma,4));
	double e = kinetic/wavefunction(x,mu,sigma) + Potential(x);

	return e; 
}


double Potential(double x){
	return pow(x,4)-2.5*pow(x,2);
}


//function useful for monte carlo


double Integrate(int nsample,double delta,double mu ,double sigma){
	//This function, useful for optimization methods, calculates only the E_loc after nsamble of monte carlo steps, without error
	vector<double> x(nsample,0);
	int acc=0,att=0;
	double r,p;
	double sum=0;

	x[0]=0; // starting point for the metropolis algorithm


	//Sample the point
	for(int i=0; i<nsample-1;i++){
		att++;
		xold=x[i];
		xnew=xold-delta+2*delta*rnd.Rannyu();

		p=squared_wave(xnew,mu,sigma)/squared_wave(xold,mu,sigma);
		r=rnd.Rannyu();

		if(p>=r){
			acc++;
			x[i+1]=xnew;
		}

		else if(p<r){
			x[i+1]=x[i];
		}

	}

	//Now calculate the function

	for(int j=0;j<nsample;j++){
		sum+=E_loc(x[j],mu,sigma);
	}

	sum = sum/nsample;

	return sum;
}


/*
void Optimize(double nopt,double l_rate,double guess_mu, double guess_sigma){
	//This function tries to implement a gradient descent algorithm by calculating the derivative of the target function
	// as an incremental ratio with h fixed
	double t_mu, t_sigma,o_mu=5,o_sigma=5;
	double h=0.1;
	double e,m_e,s_e;
	int count=0;

	double sample=10000;
	double th=0.000001;

	t_mu=  guess_sigma;
	t_sigma= guess_sigma;


	while(abs(t_mu-o_mu)>th && abs(t_sigma-o_sigma)>th && count < nopt){
		//Calculate function and derivatives
		e=Integrate(sample,delta,t_mu,t_sigma);
		m_e=Integrate(sample,delta,t_mu+h,t_sigma);
		s_e=Integrate(sample,delta,t_mu,t_sigma+h);

		//now update the parameters
		o_mu=t_mu;
		o_sigma=t_sigma;

		t_mu= t_mu -l_rate*(m_e-e)/h;
		t_sigma= t_sigma -l_rate*(s_e-e)/h;

		count++;

		cout<<count<<endl;
		cout<<"Current mu: "<<t_mu<<endl;
		cout<<"Current sigma: "<<t_sigma<<endl;
		cout<<"-----------------"<<endl;
	}

	cout<<"Best mu: "<<t_mu<<endl;
	cout<<"Best sigma: "<<t_sigma<<endl;
	//best_sigma=t_sigma;
	//best_mu=t_mu;

	return;
}
*/

void Simulated_annealing(double nopt, double temp,double guess_mu, double guess_sigma){
	double t_mu,t_sigma;
	double low_mu=0,low_sigma=0;
	double t_e;

	double gs=0;
	int sample=10000;
	int count=0;
	double h=0.01;

	t_mu=guess_mu;
	t_sigma=guess_sigma;

	for(int i=0;i<nopt;i++){
		t_e=Integrate(sample,delta,t_mu,t_sigma);
		//cout<<gs<<endl;
		if(t_e<gs){
			low_mu=t_mu;
			low_sigma=t_sigma;
			gs=t_e;
		}
		
		else if(t_e>gs){
			double r=rnd.Rannyu();
			double th=exp(-(double)count/temp);
			//cout<<"T: "<<th<<endl;
			//cout<<"R: "<<r<<endl;
			//cout<<"-------"<<endl;
			if(r<th){
				low_mu=t_mu;
				low_sigma=t_sigma;
				gs=t_e;
			}
		}
		
		count++;
		t_mu=abs(t_mu+rnd.Gauss(0,h));
		t_sigma=abs(t_sigma+rnd.Gauss(0,h));

	}

	cout<<"Approximation of ground state energy: "<<gs<<endl;

	best_sigma=low_sigma;
	best_mu=low_mu;
}

void Input(){
	ifstream ReadInput;

	cout<< "---------------------------------------------------------------"<<endl;
	cout<< "Variational Monte Carlo method to calculate ground state energy"<<endl;
	cout<< "---------------------------------------------------------------"<<endl;
	cout<<endl;

	//Read seed for random numbers
   int p1, p2;
   ifstream Primes("Primes");
   Primes >> p1 >> p2 ;
   Primes.close();

   ifstream input("seed.in");
   input >> seed[0] >> seed[1] >> seed[2] >> seed[3];
   rnd.SetRandom(seed,p1,p2);
   input.close();
  
	//Read input informations
  	ReadInput.open("input.dat");

  	ReadInput>>nblk;
  	ReadInput>>nsteps;

  	cout<<"Number of blocks: "<<nblk<<endl;
  	cout<<"Number of steps per block: "<<nsteps<<endl;
  	nsample=nblk*nsteps;
  	cout<<"Total sampling steps: "<<nsample<<endl;

  	ReadInput>>delta;

  	cout<<"Delta for metropolis algorithm: "<<delta<<endl;

  	cout<< "---------------------------------------------------------------"<<endl;
  	cout<<endl;

	return;
}

vector<double> Metropolis_sampling(double nsample,double delta, double x0,double mu, double sigma){
	
	vector<double> x(nsample,0);
	int acc=0,att=0;
	double r,p;

	x[0]=x0; // starting point for the metropolis algorithm


	//Sample the point
	for(int i=0; i<nsample-1;i++){
		att++;
		xold=x[i];
		xnew=xold-delta+2*delta*rnd.Rannyu();

		p=squared_wave(xnew,mu,sigma)/squared_wave(xold,mu,sigma);
		r=rnd.Rannyu();

		if(p>=r){
			acc++;
			x[i+1]=xnew;
		}

		else if(p<r){
			x[i+1]=x[i];
		}

	}

	cout<<"Metropolis algorithm ended"<<endl;
	cout<<"Acceptance rate = "<< (double)acc/(double)att<<endl;
	cout<<"--------------------------"<<endl;

	return x;
}




//Functions for data blocking


double error(vector<double> & av,vector<double> & av2, int n){  // Error function used to calculate
                                                                //statistical uncertainties
    if (n==0){
      return 0;
    }
    else{
      return sqrt((av2[n]-pow(av[n],2))/n);
    }
}

vector<double> Blocking(int nblocks,vector<double> & v){
 vector<double> b_v;
 int N= v.size();
 int L= N/nblocks;
 for(int i=0;i<nblocks;i++){
     double sum=0;
    for(int j=0;j<L;j++){
      int k=j+i*L;
      sum += v[k];
      }
    b_v.push_back(sum/L);
  }
  return b_v;
}

void DataBlocking(vector <double> value, vector <double> &sum_prog, vector <double> &err_prog){
    int N=sum_prog.size();
    vector <double> su2_prog(N,0);
    
    for(int i=0;i<N;i++){
        for(int j=0;j<i+1;j++){
            sum_prog[i]+=value[j];
            su2_prog[i]+=pow(value[j],2);
        }
        sum_prog[i]=sum_prog[i]/(i+1);      //Cumulative average
        su2_prog[i]=su2_prog[i]/(i+1);      //Cumulative square average
        err_prog[i] = error(sum_prog,su2_prog,i);
    }
}