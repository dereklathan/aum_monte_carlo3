#ifndef CUBE_H
#define CUBE_H
#include <iostream>
#include <ctime>
#include <cstdlib>
#include <cmath>
#include "Atom.h"
#include <vector>

class Cube{
	private:
		Atom *** atomlocation;
		int population;
		int left_to_attempt;
		int domain_x, domain_y, domain_z;
		bool all_attempted();
		int temperature;
		double ** interaction_factor;
		int fixedcount;

	public:
		bool attempted[];
		int num_threads;		
		Cube();
		//constructor		
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
		void set_interaction_factor(int);
		double get_temp();
		//returns temperature of cube
		Atom insert_atom(Atom&);
		//inserts atom into random location in cube and returns it
		bool insert_atom(Atom&, int, int, int);
		//attempts to insert atom into specified location if possible. returns true if successful. false otherwise.
		Atom get_atom(int, int, int);
		//returns atom from specified coordinates
		void set_occupy_space(int, int, int, bool);
		//given x,y,z coordinates respectively, sets space occupied(true) or vacant(false)
		bool get_occupy_space(int, int, int);
		//checks if a space is occupied given x,y,z coordinates respectively
		int get_population();
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

};
#endif
