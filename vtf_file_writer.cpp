#include "vtf_file_writer.h"
using namespace std;

vtf_file_writer::vtf_file_writer(string filename){
	outfile.open((filename + ".vtf").c_str(), ios::out);
}

void vtf_file_writer::close_file(){outfile.close();}

void vtf_file_writer::define_atom(Atom atom){
	outfile << "atom " << atom.get_index() << " radius 0.5 " << "type " << atom.get_name() << endl;
}

void vtf_file_writer::define_atoms(int m, int n, string type){
	outfile << "atom " << m << ":" << n << " radius 0.5 type " << type << endl;
}

void vtf_file_writer::write_timestep(Cube cube){
	outfile << "timestep indexed" << endl;
	outfile << "pbc " << cube.get_domain_x() << ".0 " << cube.get_domain_y() << ".0 " << cube.get_domain_z() << ".0" << endl;
	for(int c=0;c<cube.get_domain_x();c++){
		for(int d=0;d<cube.get_domain_y();d++){
			for(int e=0;e<cube.get_domain_z();e++){
				if(cube.get_atom(c,d,e).get_exists()){
					outfile << cube.get_atom(c,d,e).get_index() << " " << cube.get_atom(c,d,e).get_x_pos() << ".0 " << cube.get_atom(c,d,e).get_y_pos() << ".0 " << cube.get_atom(c,d,e).get_z_pos() << ".0" << endl;
				}
			}
		}
	}
	
}



