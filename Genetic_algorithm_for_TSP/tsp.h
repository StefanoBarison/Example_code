#ifndef TSP
#define TSP

#include<iostream>
#include<fstream>

#define _USE_MATH_DEFINES
#include<cmath>

//Useful libraries

#include<vector>
#include<algorithm>
#include<iterator>
#include<random>
#include<string>

//Random generator
//extern std::random_device radm;
//extern std::mt19937_64 mt;

class city{

public:
	city(double x, double y){
		_x=x;
		_y=y;
	}

	city():city(0,0){}

	~city();

	double Get_x();
	double Get_y();
	void Print();
private:
	double _x;
	double _y;
};


class map{

public:

	map(){
		_ncities=0;
	}

	map(int n){
		_ncities=n;
	}
	~map(){
		_cities.clear();
		_distances.clear();
	}

	void Circle_initialize();
	void Square_initialize();
	void File_initialize(std::istream & indata);
	void Create_d_matrix();

	double Distance(int i, int j);
	double Get_size();

	void Print_cities();

	void Print_d_matrix();


private:
	int _ncities;
	std::vector<std::vector<double>> _distances;
	std::vector<city> _cities;

};

class individual{
	
public:
	
	//individual();
	individual(std::vector<int> data){
		_chromosome=data;
		_lenght=0.0;
	}

	~individual();

	bool operator ==(individual & ind2);
	bool operator <(individual & ind2);

	void Evaluate(map cities);

	double Get_lenght();

	void Print();

	void Print_lenght();

	std::vector<int> Get_genes();

	void Set_genes(std::vector<int> new_genes);

	bool Check();

	void Swap_mutate();

	void Push_back_mutate(int n);

	void Multi_swap_mutate(int n);

	void Uniform_swap_mutate(double p_u);

private:
	std::vector<int> _chromosome;
	double _lenght;

};


class population{

public:
	population(){
		_size=0;
		_ncities=0;
	}
	population(int n, map cities){
		_size=n;
		_ncities=cities.Get_size();
	}
	~population(){
		_pop.clear();
	}

	void Initialize();
	void Evaluate_all(map cities);

	individual* Get_individual(int i);
	double Get_size(); 

	void Print_best(map cities);
	individual* Get_best(map cities);
	double Get_best_average(map cities);

	void Sort(map cities); //A function to sort the population from the shortest path to the longest

	void Wheel_selection(map cities); //A function to perform the roulette wheel selection

	void Elitism(map cities,int elite); //A function to perform selection saving always the best candidates of a population,
							  // ie the best candidates are always passed to the next generation

	void Evolutive_step(map cities, std::string type);

private:
	int _size;
	int _ncities;
	std::vector<individual> _pop;

};



// Useful functions


void Crossover(individual * t1,individual * t2);


#endif