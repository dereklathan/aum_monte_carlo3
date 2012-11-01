#ifndef CUBE_H
#define CUBE_H
#include <iostream>
#include <fstream>
#include <ctime>
#include <cstdlib>
#include <cmath>
#include "Atom.h"
#include "nanoparticle.h"
#include <vector>

class Cube{
	private:
		Atom *** atomlocation;
		Nanoparticle *** nparticlelocation;
		double population;
		int left_to_attempt;
		int domain_x, domain_y, domain_z;
		bool all_attempted();
		double temperature;
		double T;
		double ** interaction_factor;
		int fixedcount;
		double flux_in, flux_out, E1, E2, dE;
		void move_nanoparticles();
		void calc_rms(vector<Atom>);
		vector<double> x_rms;
		vector<double> y_rms;
		vector<double> z_rms;

		//moves all nanoparticles

	public:		
		Cube();
		//constructor	
		void clear();
		//clears the cube
		void set_domain(int, int, int);
		//sets size of cube
		int get_domain_x();
		//returns length of x axis
		int get_domain_y();
		//returns length of y axis
		int get_domain_z();
		//returns length of z axis
		int get_volume();
		//returns volume of cube
		void set_temp(double t);
		//set temperature of cube
		double ** set_interaction_factor(int);
		void set_interaction_factor(double**);
		double get_temp();
		//returns temperature of cube
		Atom insert_atom(Atom&);
		//inserts atom into random location in cube and returns it
		Nanoparticle insert_nanoparticle(Nanoparticle&);
		//inserts nanoparticle into random location in cube and returns it
		bool insert_atom(Atom&, int, int, int);
		//attempts to insert atom into specified location if possible. returns true if successful. false otherwise.
		bool insert_nanoparticle(Nanoparticle&, int, int, int);
		//attempts to insert nanoparticle into specified location if possible. returns true if successful. false otherwise.
		Atom get_atom(int, int, int);
		//returns atom from specified coordinates
		void set_occupy_space(int, int, int, bool);
		//given x,y,z coordinates respectively, sets space occupied(true) or vacant(false)
		bool get_occupy_space(int, int, int);
		//checks if a space is occupied given x,y,z coordinates respectively
		double get_population();
		//returns number of atoms in cube
		void advance_timestep_pbc();
		double calculate_pot_energy_pbc();
		//calculates potential energy of system for periodic
		//boundary conditions
		double calculate_pot_energy_pbc(int,int,int,int,int,int);
		//calculates potential energy of system for periodic boundary 
		//conditions and returns it. first 3 integers are coordinates 
		//for particle. second 3 are where the particle will be when 
		//energy is calculated.
		bool can_still_move();
		double calculate_pot_energy_opentop();
		//calculates potential energy of system for periodic boundary
		//conditions on x and y axes and open end on z axis.
		double calculate_pot_energy_opentop(int, int, int, int, int, int);
		//calculates potential energy after particle at x1,y1,z1 moves
		//to x2,y2,z2
		double calculate_dE_opentop(int, int, int, int, int, int);
		//calculates difference in energy when particle at x1,y1,z1
		//moves to x2,y2,z2
		void advance_timestep_opentop();
		bool obc_eq();
		//returns true if in equilibrium for opentop boundary conditions
		void seed_random(int);
		//seed random number generator with int
		vector<string> toAppend();
		//return vector of newly created atoms for defining in vtf file
		double get_flux_in();
		//number of particles that entered domain after previous timestep;
		double get_flux_out();
		//number of particles that left the domain after previous timestep;
		double get_x_rms(int);
		double get_y_rms(int);
		double get_z_rms(int);

};
#endif
