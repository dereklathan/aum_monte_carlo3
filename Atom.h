#ifndef ATOM_H
#define ATOM_H
#include <iostream>
#include <string>
using namespace std;

class Atom{
	private:
		unsigned int x_pos, y_pos, z_pos;
		string name;
		unsigned int mass;
		unsigned int grav_pot_energy;
		unsigned int internal_pot_energy;
		bool attempted;
		bool exists;
		bool fixed;
		double strength;
		int index;
		int atom_type_num;
		bool replaceable;

	public:
		Atom();
		//constructor
		void set_name(string);
		//sets name of atom
		string get_name();
		//returns name of atom
		void set_mass(unsigned int);
		//sets mass of atom
		unsigned int get_mass();
		//gets mass of atom
		void set_fixed(bool);
		//sets atom to fixed position
		bool is_fixed();
		//returns true if fixed position. false otherwise
		void set_x_pos(unsigned int);
		//sets x coordinate of atom center
		void set_y_pos(unsigned int);
		//sets y coordinate of atom center
		void set_z_pos(unsigned int);
		//sets z coordinate of atom center
		unsigned int get_x_pos();
		//gets x coordinate of atom center
		unsigned int get_y_pos();
		//gets y coordinate of atom center
		unsigned int get_z_pos();
		//gets z coordinate of atom center
		void set_grav_pot_energy();
		//calculates grav_pot_energy based on mass and z_pos
		unsigned int get_grav_pot_energy();
		//returns grav_pot_energy
		void set_strength(double);
		//sets interaction strength factor
		double get_strength();
		//returns interaction strength factor
		void set_exists(bool);
		//if exists should be set to true, false otherwise
		bool get_exists();
		//checks if the particle exists
		void set_attempted(bool);
		//if attempt to move has been made for timestep, should be true. false otherwise
		bool get_attempted();
		//returns true if attempt has been made, false otherwise
		void set_index(int);
		//set unique integer for tracking
		int get_index();
		//returns unique integer for tracking
		void set_type_num(int);
		//set type number for interaction factor matrix
		int get_type_num();
		//returns type number for interaction factor matrix
		void set_replaceable(bool);
		bool get_replaceable();
		
};

#endif
