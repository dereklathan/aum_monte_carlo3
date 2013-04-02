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
#include <fstream>
#include <cstring>
#include "Atom.h"
#include "nanoparticle.h"
#include "cube.h"
#include "infile_reader.h"
#include "vtf_file_writer.h"
#include "dat_file_writer.h"

int main(int argc, char *argv[]){
	srand(time(NULL));
	Infile_reader infile_reader((string)argv[1]);
	infile_reader.setData();
	cout << infile_reader.get_nparticle_move_count() << endl;
	cout << infile_reader.getpercentdomainfill()[1] << endl;
	cout << infile_reader.get_nanoparticle_temp()[1] << endl;
	vtf_file_writer vtf_writer(infile_reader.getoutfilename());
	dat_file_writer dat_writer(infile_reader.getoutfilename());
	int graph_interval=0;
	int * atom_layer_count;
	int tcount;
	double z_cent_mass[infile_reader.gettimesteps()];
	double flux_in[infile_reader.gettimesteps()];
	double flux_out[infile_reader.gettimesteps()];
	Cube cube;
	cout << infile_reader.getunitcellsize()[0] << " " << infile_reader.getunitcellsize()[1] << " " <<  infile_reader.getunitcellsize()[2] << endl;
	cube.set_domain(infile_reader.getunitcellsize()[0], infile_reader.getunitcellsize()[1], infile_reader.getunitcellsize()[2]);
	cube.set_temp(infile_reader.get_temp());
	cube.set_nparticle_move_count(infile_reader.get_nparticle_move_count());
	Atom atom;
	Nanoparticle nanoparticle;
	for(int c=0;c<infile_reader.get_particle_types();c++){
		if(!infile_reader.is_nanoparticle()[c]){
			atom.set_name(infile_reader.get_particle_name()[c]);
			atom.set_type_num(c);
			atom.set_strength(infile_reader.get_strength()[c]);
			atom.set_mass(infile_reader.get_mass()[c]);
			atom.set_fixed(infile_reader.is_fixed()[c]);
			for(double d=0;d<((double)cube.get_domain_x()*(double)cube.get_domain_y()*infile_reader.getpercentdomainfill()[c])/100;d++){
				atom.set_replaceable(true);
				while(!cube.insert_atom(atom, rand()%cube.get_domain_x(), rand()%cube.get_domain_y(), 0));
			}		
		}
		else{
			nanoparticle.set_T(infile_reader.get_nanoparticle_temp()[c]);
			for(double d=0;d<((double)cube.get_domain_x()*(double)cube.get_domain_y()*(double)cube.get_domain_z()*infile_reader.getpercentdomainfill()[c])/100;d++)
				cube.insert_nanoparticle(nanoparticle);
		}
	}
	vtf_writer.define_atoms(0,10000,"particle");
	vtf_writer.define_nanoparticles(0,10000,"particle");
	vtf_writer.write_timestep(cube);
	dat_writer.write_timestep_graph_z(cube,0);
	cube.set_interaction_factor(infile_reader.get_particle_types()-infile_reader.get_nanoparticle_types());
//	int c=0;
//	while(!cube.obc_eq()){
	for(int c=0;c<infile_reader.gettimesteps();c++){
		z_cent_mass[c]=0;
		cout << "timestep " << c+1 << endl;
		cube.advance_timestep_opentop();
		atom_layer_count=dat_writer.get_timestep_array_z(cube);
		flux_in[c]=cube.get_flux_in();
		flux_out[c]=cube.get_flux_out();
		for(int d=0;d<cube.get_domain_z();d++){
			z_cent_mass[c]+=d*atom_layer_count[d];
		}
		cout << "total population: " << cube.get_population() << endl;
		z_cent_mass[c]/=(double)cube.get_population();
		vtf_writer.write_timestep(cube);
		graph_interval++;
		if(graph_interval==infile_reader.get_graph_interval()){
			dat_writer.write_timestep_graph_z(cube, c+1);
			graph_interval=0;
		}
		tcount=c+1;
	}

	ofstream z_cent_mass_file;
	ofstream flux_in_file;
	ofstream flux_out_file;
	flux_in_file.open((infile_reader.getoutfilename() + "_flux_in.dat").c_str(),ios::out);
	flux_out_file.open((infile_reader.getoutfilename() + "_flux_out.dat").c_str(),ios::out);
	z_cent_mass_file.open((infile_reader.getoutfilename() + "_z_cent_mass.dat").c_str(),ios::out);
	for(int c=0;c<tcount;c++){
		z_cent_mass_file << c+1 << " " << z_cent_mass[c] << endl;
		flux_in_file << c+1 << " " << flux_in[c] << endl;
		flux_out_file << c+1 << " " << flux_out[c] << endl;
	}
	z_cent_mass_file.close();
	flux_in_file.close();
	flux_out_file.close();
	return 0;
}

