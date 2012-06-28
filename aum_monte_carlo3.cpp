//Author: Derek Lathan
/*
this simulation is similar to the previous with the exception of 
boundary conditions. periodic boundary conditions still apply to the
x and y axes, but the percentage to fill the domain only applies to 
the bottom layer of the z axis. if the particle tries to move below
the z axis then it will remain in the same location. if it moves up
then another particle will take its place keeping the fill percentage
of the bottom layer constant. if it move past the top of the z-axis,
the particle will escape. diffusion.
*/

#include <iostream>
#include <ctime>
#include <cstdlib>
#include "Atom.h"
#include "cube.h"
#include "infile_reader.h"
#include "vtf_file_writer.h"
#include "dat_file_writer.h"

int main(){
	srand(time(NULL));
	Infile_reader infile_reader("infile");
	infile_reader.setData();
	vtf_file_writer vtf_writer(infile_reader.getoutfilename());
	dat_file_writer dat_writer(infile_reader.getoutfilename());
	int graph_interval=0;
	Cube cube;
	cube.set_domain(infile_reader.getunitcellsize()[0], infile_reader.getunitcellsize()[1], infile_reader.getunitcellsize()[2]);
	cube.set_temp(infile_reader.get_temp());
	Atom atom;
	for(int c=0;c<infile_reader.get_particle_types();c++){
		atom.set_name(infile_reader.get_particle_name()[c]);
		atom.set_type_num(c);
		atom.set_strength(infile_reader.get_strength()[c]);
		atom.set_mass(infile_reader.get_mass()[c]);
		atom.set_fixed(infile_reader.is_fixed()[c]);
		for(int d=0;d<((float)cube.get_domain_x()*(float)cube.get_domain_y()*infile_reader.getpercentdomainfill()[c])/100;d++){
			while(!cube.insert_atom(atom, rand()%cube.get_domain_x(), rand()%cube.get_domain_y(), 0));
			//vtf_writer.define_atom(atom);
		}		
	}
	vtf_writer.define_atoms(0,10000,"particle");
/*	for(int c=0;c<cube.get_domain_x();c++){
		for(int d=0;d<cube.get_domain_y();d++){
			for(int e=0;e<cube.get_domain_z();e++){
				if(cube.get_atom(c,d,e).get_exists())
					cout << cube.get_atom(c,d,e).get_index() << " " << cube.get_atom(c,d,e).get_x_pos() << " " << cube.get_atom(c,d,e).get_y_pos() << " " << cube.get_atom(c,d,e).get_z_pos() << endl;
			}
		}
	}
*/	vtf_writer.write_timestep(cube);
	dat_writer.write_timestep_graph_z(cube,0);
	cube.set_interaction_factor(infile_reader.get_particle_types());
	for(int c=0;c<infile_reader.gettimesteps();c++){
		cube.advance_timestep_opentop();
		vtf_writer.write_timestep(cube);
		graph_interval++;
		if(graph_interval==infile_reader.get_graph_interval()){
			dat_writer.write_timestep_graph_z(cube, c+1);
			graph_interval=0;
		}
	}
	return 0;
}

