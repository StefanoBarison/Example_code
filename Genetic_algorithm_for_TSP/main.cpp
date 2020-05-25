#include "tsp.h"

#include<iostream>
#include<fstream>
#include<vector>
#include<cmath>
#include<random>
#include<string>

using namespace std;


int main (){

	//Constants for the GA
	int n_cities=32;
	int n_ind=1000;
	int n_steps=1000;
	string type ("Elitism");
	ifstream indata;

	//initialize the map with n_cities random cities on a circumference and printing them

	indata.open("square.dat"); 

	map m(n_cities);
	//m.Circle_initialize();
	//m.Square_initialize();

	m.File_initialize(indata);
	m.Print_cities();
	m.Create_d_matrix();
	m.Print_d_matrix();

	//Create a population to perform the GA
	population p(n_ind,m);
	p.Initialize();
	p.Evaluate_all(m);

	//Vector to save the results

	vector<double> best_result;
	vector<double> best_ave;

	//Perform the GA
	for(int i=0;i<n_steps;i++){
		cout<<"Step "<<i+1<<endl<<endl;
		best_result.push_back(p.Get_best(m)->Get_lenght());
		best_ave.push_back(p.Get_best_average(m));
		p.Evolutive_step(m,type);

		if(i==n_steps-1){
			p.Print_best(m);
		}
	}


	//Save the data
	ofstream outdata1,outdata2;

	outdata1.open("results/test.dat");
	outdata2.open("results/test_average.dat");

	for(unsigned int i=0;i<best_result.size();i++){
		outdata1<<i<<","<<best_result[i]<<endl;
		outdata2<<i<<","<<best_ave[i]<<endl;
	}



	return 0;
	
}