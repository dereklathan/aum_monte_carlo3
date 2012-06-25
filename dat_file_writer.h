//this class is for writing *.dat graph files for xmgrace
#ifndef DAT_FILE_WRITER_H
#define DAT_FILR_WRITER_H
#include <iostream>
#include <fstream>
#include <string>
#include "cube.h"
#include "Atom.h"
#include <sstream>
#include <cstring>
using namespace std;
class dat_file_writer{
	private:
		ofstream outfile;
		string filename;
	public:
		dat_file_writer(string);
		//constructor. string is filename from infile
		void write_timestep_graph_x(Cube, int);
		//creates graph of number of particles vs layer of cube along x-axis int is timestep number
		void write_timestep_graph_y(Cube, int);
		//creates graph of number of particles vs layer of cube along y-axis int is timestep number
		void write_timestep_graph_z(Cube, int);
		//creates graph of number of particles vs layer of cube along z-axis int is timestep number
	

	
};
#endif
