#include "dat_file_writer.h"
using namespace std;

dat_file_writer::dat_file_writer(string outfilename){filename=outfilename;}

void dat_file_writer::write_timestep_graph_z(Cube cube, int timestep){
	int layer_atom_count;
	stringstream convert;
	string timestep_number_string;
	convert << timestep;
	timestep_number_string = convert.str();
	for(int c=timestep_number_string.size();c<5;c++)
		timestep_number_string.insert(0,"0");
	outfile.open((filename + "_" + timestep_number_string + "_z.dat").c_str(), ios::out);
	for(int c=0;c<cube.get_domain_z();c++){
		layer_atom_count=0;
		for(int d=0;d<cube.get_domain_x();d++){
			for(int e=0;e<cube.get_domain_y();e++){
				if(cube.get_atom(d, e, c).get_exists())
					layer_atom_count++;
			}
		}
		outfile << c << " " << layer_atom_count << endl;
	}
	outfile.close();
}
