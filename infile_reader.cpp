#include "infile_reader.h"
#include <iostream>
#include <fstream>
#include <string>
using namespace std;

Infile_reader::Infile_reader(string filename){
	reader.open(filename.c_str(), ifstream::in);	
}

void Infile_reader::setData(){
	nanoparticle_types=0;
	outfilename="";
	while(reader.peek()!='#')
		outfilename += (char)reader.get();
	while(reader.peek()!='\n')
		reader.get();
	reader >> num_sims;
	while(reader.peek()!='\n')
		reader.get();
	for(int c=0;c<3;c++)
		reader >> unitcellsize[c];
	while(reader.peek()!='\n')
		reader.get();
	reader >> temperature;
	while(reader.peek()!='\n')
		reader.get();
	reader >> timesteps;
	while(reader.peek()!='\n')
		reader.get();
	reader >> graph_interval;
	while(reader.peek()!='\n')
		reader.get();
	reader >> particle_types;
	particle_name = new string [particle_types];
	for(int c=0;c<particle_types;c++)
		particle_name[c]="";
	mass = new unsigned int [particle_types];
	percentdomainfill = new float [particle_types];
	strength = new int [particle_types];
	fixed = new bool [particle_types];
	n_particle = new bool [particle_types];
	n_particle_temp = new double [particle_types];
	while(reader.peek()!='\n')
		reader.get();
	reader >> boundary_conditions[0];
	while(reader.peek()!='\n')
		reader.get();
	reader >> boundary_conditions[1];
	while(reader.peek()!='\n')
		reader.get();
	reader >> boundary_conditions[2];
	while(reader.peek()!='\n')
		reader.get();
	reader.get();
	for(int c=0;c<particle_types;c++){
		while(reader.peek()!='\n')
			reader.get();
		reader.get();
		while(reader.peek()!='#')
			particle_name[c] += (char)reader.get();
		while(reader.peek()!='\n')
			reader.get();
		reader >> mass[c];
		while(reader.peek()!='\n')
			reader.get();
		reader >> percentdomainfill[c];
		while(reader.peek()!='\n')
			reader.get();
		reader >> strength[c];
		while(reader.peek()!='\n')
			reader.get();
		reader >> fixed[c];
		while(reader.peek()!='\n')
			reader.get();
		reader >> n_particle[c];
		if(n_particle[c]==true)
			nanoparticle_types++;
		while(reader.peek()!='\n')
			reader.get();
		reader >> n_particle_temp[c];
		while(reader.peek()!='\n')
			reader.get();
		reader.get();

	}
	reader.close();
}

int* Infile_reader::getunitcellsize(){
	int* pointer = unitcellsize;
	return pointer;
}

float* Infile_reader::getpercentdomainfill(){
	float* pointer = percentdomainfill;
	return pointer;
}

int Infile_reader::gettimesteps(){
	return timesteps;
}

string Infile_reader::getoutfilename(){return outfilename;}

int Infile_reader::get_graph_interval(){return graph_interval;}

int Infile_reader::get_particle_types(){return particle_types;}

int Infile_reader::get_nanoparticle_types(){return nanoparticle_types;}

string* Infile_reader::get_particle_name(){
	string* pointer = particle_name;
	return pointer;
}

unsigned int* Infile_reader::get_mass(){
	unsigned int* pointer = mass;
	return pointer;
}

int* Infile_reader::get_strength(){
	int* pointer = strength;
	return pointer;
}

bool* Infile_reader::is_fixed(){
	bool* pointer = fixed;
	return pointer;
}

bool* Infile_reader::is_nanoparticle(){
	bool* pointer = n_particle;
	return pointer;
}

double* Infile_reader::get_nanoparticle_temp(){
	double* pointer = n_particle_temp;
	return pointer;
}

int* Infile_reader::get_boundary_conditions(){
	int* pointer = boundary_conditions;
	return pointer;
}

double Infile_reader::get_temp(){return temperature;}

int Infile_reader::get_num_sims(){return num_sims;}

