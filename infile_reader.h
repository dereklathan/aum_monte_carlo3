#ifndef INFILE_READER_H
#define INFILE_READER_H
#include <iostream>
#include <fstream>
#include <string>
using namespace std;

class Infile_reader{
	private:
		ifstream reader;
		string outfilename;
		int unitcellsize[3];
		float * percentdomainfill;
		int timesteps;
		int graph_interval;
		unsigned int particle_types;
		unsigned int * mass;
		string * particle_name;
		bool * fixed;
		int * strength;
		int boundary_conditions[3];
		double temperature;
		int num_sims;
		
	public:
		Infile_reader();
		//default constructor
		Infile_reader(string);
		//constructor. string is filename to be read
		void setData();
		//sets data based on input file
		string getoutfilename();
		//returns string for output vtf filename
		int* getunitcellsize();
		//returns array containing x,y,z params for unit cell respectively
		float* getpercentdomainfill();
		//returns percentage of domain to fill with atoms
		int gettimesteps();
		//returns number of timesteps 	
		int* get_boundary_conditions();
		//returns array of integers representing boundary conditions
		//for x,y,z axes respectively
		int get_graph_interval();
		//returns how many timesteps should pass before writing graph
		int get_particle_types();
		//how many particle types should be defined
		string* get_particle_name();
		//returns array of particle names
		unsigned int* get_mass();
		//returns array of particle masses;
		int* get_strength();
		//returns array of particle interaction strengths;
		bool* is_fixed();
		//returns bool array that tells if location of particle is fixed
		double get_temp();
		//returns temperature of cube
		int get_num_sims();
		//returns number of simulations to run

		
};
#endif
