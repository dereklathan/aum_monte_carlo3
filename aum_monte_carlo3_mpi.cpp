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
#include <cmath>
#include <fstream>
#include <cstdlib>
#include <string>
#include <cstring>
#include <sstream>
#include <mpi.h>
#include "Atom.h"
#include "cube.h"
#include "infile_reader.h"
#include "vtf_file_writer.h"
#include "dat_file_writer.h"

int main(int argc, char *argv[]){
	int rank, size, len, result, random;

	int * atom_layer_count;
	char hostname[MPI_MAX_PROCESSOR_NAME];
	Cube cube;
	MPI_Status status;
	MPI_Init(&argc, &argv);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &size);
	MPI_Get_processor_name(hostname, &len);
	if(rank==0){
		srand(time(NULL));
		for(int c=1;c<size;c++){
			random=rand();
			MPI_Send(&random, 1, MPI_INT, c, 98, MPI_COMM_WORLD);
		}
	}
	else{
		MPI_Recv(&random, 1, MPI_INT, 0, 98, MPI_COMM_WORLD, &status);
		srand(random);
	}
	cube.seed_random(rand());
	Infile_reader infile_reader(/*(string)argv[1]*/"infile");
	infile_reader.setData();
	cube.set_nparticle_move_count(infile_reader.get_nparticle_move_count());
	int counts_for_rank[infile_reader.gettimesteps()/infile_reader.get_graph_interval()][cube.get_domain_z()][infile_reader.get_num_sims()/size];
	double z_cent_mass[infile_reader.gettimesteps()][infile_reader.get_num_sims()/size];
	double flux_in[infile_reader.gettimesteps()][infile_reader.get_num_sims()/size];
	double flux_out[infile_reader.gettimesteps()][infile_reader.get_num_sims()/size];
	double int_factor[infile_reader.get_particle_types()-infile_reader.get_nanoparticle_types()][infile_reader.get_particle_types()-infile_reader.get_nanoparticle_types()];
	cube.set_domain(infile_reader.getunitcellsize()[0], infile_reader.getunitcellsize()[1], infile_reader.getunitcellsize()[2]);
	dat_file_writer dat_writer(infile_reader.getoutfilename());
	cube.set_temp(infile_reader.get_temp());
	Atom atom;
	Nanoparticle nanoparticle;
	if(rank==0){
		double ** temp_int_factor=cube.set_interaction_factor(infile_reader.get_particle_types()-infile_reader.get_nanoparticle_types());
		for(int c=0;c<infile_reader.get_particle_types()-infile_reader.get_nanoparticle_types();c++){
			for(int d=0;d<infile_reader.get_particle_types()-infile_reader.get_nanoparticle_types();d++)
				int_factor[c][d]=temp_int_factor[c][d];
		}
		for(int c=1;c<size;c++)
			MPI_Send(&int_factor, (infile_reader.get_particle_types()-infile_reader.get_nanoparticle_types())*(infile_reader.get_particle_types()-infile_reader.get_nanoparticle_types()), MPI_DOUBLE, c, 97, MPI_COMM_WORLD);
	}
	else{
		double ** temp_int_factor = new double*[infile_reader.get_particle_types()-infile_reader.get_nanoparticle_types()];
		for(int c=0;c<infile_reader.get_particle_types()-infile_reader.get_nanoparticle_types();c++)
			temp_int_factor[c] = new double[infile_reader.get_particle_types()-infile_reader.get_nanoparticle_types()];
		MPI_Recv(&int_factor, (infile_reader.get_particle_types()-infile_reader.get_nanoparticle_types())*(infile_reader.get_particle_types()-infile_reader.get_nanoparticle_types()), MPI_DOUBLE, 0, 97, MPI_COMM_WORLD, &status);
		for(int c=0;c<infile_reader.get_particle_types()-infile_reader.get_nanoparticle_types();c++){
			for(int d=0;d<infile_reader.get_particle_types()-infile_reader.get_nanoparticle_types();d++)
				temp_int_factor[c][d] = int_factor[c][d];
		}
		cube.set_interaction_factor(temp_int_factor);
	}
	for(int n=0;n<infile_reader.get_num_sims()/size;n++){
		for(int c=0;c<infile_reader.get_particle_types();c++){
			if(!infile_reader.is_nanoparticle()[c]){
				atom.set_name(infile_reader.get_particle_name()[c]);
				atom.set_type_num(c);
				atom.set_strength(infile_reader.get_strength()[c]);
				atom.set_mass(infile_reader.get_mass()[c]);
				atom.set_fixed(infile_reader.is_fixed()[c]);
				for(float d=0;d<((float)cube.get_domain_x()*(float)cube.get_domain_y()*infile_reader.getpercentdomainfill()[c])/100;d++){
					atom.set_replaceable(true);
					while(!cube.insert_atom(atom, rand()%cube.get_domain_x(), rand()%cube.get_domain_y(), 0));
				}		
			}
			else{
				nanoparticle.set_T(infile_reader.get_nanoparticle_temp()[c]);
				for(float d=0;d<((float)cube.get_domain_x()*(float)cube.get_domain_y()*(float)cube.get_domain_z()*infile_reader.getpercentdomainfill()[c])/100;d++)
				cube.insert_nanoparticle(nanoparticle);
			}
		}
		for(int c=0;c<infile_reader.gettimesteps();c++){
			z_cent_mass[c][n]=0;
			cout << rank << ":" << c << endl;
			cube.advance_timestep_opentop();
			flux_in[c][n]=(double)cube.get_flux_in();
			flux_out[c][n]=(double)cube.get_flux_out();
			atom_layer_count=dat_writer.get_timestep_array_z(cube);
			for(int d=0;d<cube.get_domain_z();d++)
				z_cent_mass[c][n]+=d*atom_layer_count[d];
			z_cent_mass[c][n]/=(double)cube.get_population();
			if((c+1)%infile_reader.get_graph_interval()==0){
				for(int d=0;d<cube.get_domain_z();d++)
					counts_for_rank[(c+1)/infile_reader.get_graph_interval()-1][d][n]=atom_layer_count[d];
			}	
		}
		
		cube.clear();
	}
	for(int c=0;c<infile_reader.get_num_sims()/size;c++){
		for(int d=0;d<cube.get_domain_z();d++){
			if(rank==1)
				cout << rank << ": " << counts_for_rank[d][c] << endl;
			}
	}	
	if(rank!=0){
		for(int c=0;c<infile_reader.get_num_sims()/size;c++){
			for(int d=0;d<cube.get_domain_z();d++){
				for(int e=0;e<infile_reader.gettimesteps()/infile_reader.get_graph_interval();e++)
					MPI_Send(&counts_for_rank[e][d][c],1, MPI_INT, 0, 99, MPI_COMM_WORLD);
			}
			for(int d=0;d<infile_reader.gettimesteps();d++){
				MPI_Send(&z_cent_mass[d][c], 1, MPI_DOUBLE, 0, 90, MPI_COMM_WORLD);
				MPI_Send(&flux_in[d][c], 1, MPI_DOUBLE, 0, 89, MPI_COMM_WORLD);
				MPI_Send(&flux_out[d][c], 1, MPI_DOUBLE, 0, 88, MPI_COMM_WORLD);
			}
		}
	}
	else{
		int all_counts[infile_reader.gettimesteps()/infile_reader.get_graph_interval()][cube.get_domain_z()][infile_reader.get_num_sims()];
		double all_z_cent_mass[infile_reader.gettimesteps()][infile_reader.get_num_sims()];
		double all_flux_in[infile_reader.gettimesteps()][infile_reader.get_num_sims()];
		double all_flux_out[infile_reader.gettimesteps()][infile_reader.get_num_sims()];
		double avg_z_cent_mass[infile_reader.gettimesteps()];
		double avg_flux_in[infile_reader.gettimesteps()];
		double avg_flux_out[infile_reader.gettimesteps()];
		double z_cent_mass_dev[infile_reader.gettimesteps()][infile_reader.get_num_sims()];
		double z_cent_mass_std_dev[infile_reader.gettimesteps()];
		double deviations[infile_reader.gettimesteps()/infile_reader.get_graph_interval()][cube.get_domain_z()][infile_reader.get_num_sims()];
		double standard_dev[infile_reader.gettimesteps()/infile_reader.get_graph_interval()][cube.get_domain_z()];
		double average[infile_reader.gettimesteps()/infile_reader.get_graph_interval()][cube.get_domain_z()];
		stringstream convert;
		string timestep_number_string;
		for(int c=0;c<infile_reader.get_num_sims()/size;c++){
			for(int d=0;d<infile_reader.gettimesteps();d++){
				all_z_cent_mass[d][c]=z_cent_mass[d][c];
				all_flux_in[d][c]=flux_in[d][c];
				all_flux_out[d][c]=flux_out[d][c];
			}
			for(int d=0;d<cube.get_domain_z();d++){
				for(int e=0;e<infile_reader.gettimesteps()/infile_reader.get_graph_interval();e++)
					all_counts[e][d][c]=counts_for_rank[e][d][c];
			}
		}
		for(int c=1;c<size;c++){
			for(int d=0;d<infile_reader.get_num_sims()/size;d++){
				for(int e=0;e<cube.get_domain_z();e++){
					for(int b=0;b<infile_reader.gettimesteps()/infile_reader.get_graph_interval();b++)
						MPI_Recv(&counts_for_rank[b][e][d], 1, MPI_INT, c, 99, MPI_COMM_WORLD, &status);
				}
				for(int e=0;e<infile_reader.gettimesteps();e++){
					MPI_Recv(&z_cent_mass[e][d], 1, MPI_DOUBLE, c, 90, MPI_COMM_WORLD, &status);
					MPI_Recv(&flux_in[e][d], 1, MPI_DOUBLE, c, 89, MPI_COMM_WORLD, &status);
					MPI_Recv(&flux_out[e][d], 1, MPI_DOUBLE, c, 88, MPI_COMM_WORLD, &status);
				}
			}
			int r=0;
			for(int e=infile_reader.get_num_sims()/size*c;e<(infile_reader.get_num_sims()/size)*(c+1);e++){
				for(int d=0;d<infile_reader.gettimesteps();d++){
					all_z_cent_mass[d][e]=z_cent_mass[d][r];
					all_flux_in[d][e]=flux_in[d][r];
					all_flux_out[d][e]=flux_out[d][r];
				}
				for(int d=0;d<cube.get_domain_z();d++){
					for(int b=0;b<infile_reader.gettimesteps()/infile_reader.get_graph_interval();b++)
						all_counts[b][d][e]=counts_for_rank[b][d][r];
				}
				r++;
			}
		}
/*
		for(int c=0;c<infile_reader.get_num_sims();c++){
			for(int d=0;d<cube.get_domain_z();d++)
				cout << all_counts[d][c] << " ";
			cout << endl;
		}
*/
		for(int e=0;e<infile_reader.gettimesteps()/infile_reader.get_graph_interval();e++){
			for(int c=0;c<cube.get_domain_z();c++){
				average[e][c]=0;
				standard_dev[e][c]=0;
				for(int d=0;d<infile_reader.get_num_sims();d++){
					average[e][c]+=(double)all_counts[e][c][d];
					deviations[e][c][d]=all_counts[e][c][d];
				}
				average[e][c]/=(double)infile_reader.get_num_sims();
			}
		}

		for(int e=0;e<infile_reader.gettimesteps()/infile_reader.get_graph_interval();e++){
			for(int c=0;c<cube.get_domain_z();c++){
				for(int d=0;d<infile_reader.get_num_sims();d++){
					deviations[e][c][d]-=average[e][c];
					deviations[e][c][d]=pow(deviations[e][c][d],2);
					standard_dev[e][c]+=deviations[e][c][d];
				}
				standard_dev[e][c]/=infile_reader.get_num_sims()-1;
				standard_dev[e][c]=sqrt(standard_dev[e][c]);
			}
		}

		for(int c=0;c<infile_reader.gettimesteps();c++){
			avg_z_cent_mass[c]=0;
			avg_flux_in[c]=0;
			avg_flux_out[c]=0;
			z_cent_mass_std_dev[c]=0;
			for(int d=0;d<infile_reader.get_num_sims();d++){
				avg_z_cent_mass[c]+=all_z_cent_mass[c][d];
				avg_flux_in[c]+=all_flux_in[c][d];
				avg_flux_out[c]+=all_flux_out[c][d];
				z_cent_mass_dev[c][d]=all_z_cent_mass[c][d];
			}	
			avg_z_cent_mass[c]/=(double)infile_reader.get_num_sims();
			avg_flux_in[c]/=(double)infile_reader.get_num_sims();
			avg_flux_out[c]/=(double)infile_reader.get_num_sims();
		}

		for(int c=0;c<infile_reader.gettimesteps();c++){
			for(int d=0;d<infile_reader.get_num_sims();d++){
				z_cent_mass_dev[c][d]-=avg_z_cent_mass[c];
				z_cent_mass_dev[c][d]=pow(z_cent_mass_dev[c][d],2);
				z_cent_mass_std_dev[c]+=z_cent_mass_dev[c][d];
			}
			z_cent_mass_std_dev[c]/=infile_reader.get_num_sims()-1;
			z_cent_mass_std_dev[c]=sqrt(z_cent_mass_std_dev[c]);
		}
/*
		for(int c=0;c<cube.get_domain_z();c++)
			cout << average[c] << " ";
		cout << endl;
		for(int c=0;c<cube.get_domain_z();c++)
			cout << standard_dev[c]/average[c]*100 << " ";
		cout << endl;
*/
		ofstream aver;
		ofstream dev;
		ofstream z_cent_mass_file;
		ofstream z_cent_mass_dev_file;
		ofstream flux_in_file;
		ofstream flux_out_file;
		z_cent_mass_file.open((infile_reader.getoutfilename() + "_z_cent_mass.dat").c_str(), ios::out);
		z_cent_mass_dev_file.open((infile_reader.getoutfilename() + "_z_cent_mass_dev.dat").c_str(), ios::out);
		flux_in_file.open((infile_reader.getoutfilename() + "_flux_in.dat").c_str(),ios::out);
		flux_out_file.open((infile_reader.getoutfilename() + "_flux_out.dat").c_str(),ios::out);
		for(int d=0;d<infile_reader.gettimesteps()/infile_reader.get_graph_interval();d++){
			convert << (d+1)*infile_reader.get_graph_interval();
			timestep_number_string=convert.str();
			aver.open((infile_reader.getoutfilename() + "_" + timestep_number_string + ".dat").c_str(), ios::out);
			dev.open((infile_reader.getoutfilename() + "_" + timestep_number_string + "_deviations.dat").c_str(), ios::out);
			for(int c=0;c<cube.get_domain_z();c++){
				aver << c << " " << average[d][c] << endl;
				dev << c << " " << standard_dev[d][c]/average[d][c]*100 << endl;
			}
			aver.close();
			dev.close();
			convert.str("");
		}
		for(int c=0;c<infile_reader.gettimesteps();c++){
			z_cent_mass_file << c+1 << " " << avg_z_cent_mass[c] << endl;
			z_cent_mass_dev_file << c+1 << " " << z_cent_mass_std_dev[c]/avg_z_cent_mass[c]*100 << endl;
			flux_in_file << c+1 << " " << avg_flux_in[c] << endl;
			flux_out_file << c+1 << " " << avg_flux_out[c] << endl;
		}
		aver.close();
		dev.close();
		z_cent_mass_file.close();
		z_cent_mass_dev_file.close();
		flux_in_file.close();
		flux_out_file.close();

	}
	MPI_Finalize();
	return 0;
}

