#include "Atom.h"
using namespace std;

Atom::Atom(){
	exists = false;
	attempted = false;
}

void Atom::set_name(string atomname){
	name=atomname;
}

string Atom::get_name(){return name;}

void Atom::set_mass(unsigned int atom_mass){mass=atom_mass;}

unsigned int Atom::get_mass(){return mass;}

void Atom::set_fixed(bool is_fixed){fixed=is_fixed;}

bool Atom::is_fixed(){return fixed;}

void Atom::set_x_pos(unsigned int x){
	x_pos=x;
}

void Atom::set_y_pos(unsigned int y){
	y_pos=y;
}

void Atom::set_z_pos(unsigned int z){
	z_pos=z;
}

unsigned int Atom::get_x_pos(){
	return x_pos;
}

unsigned int Atom::get_y_pos(){
	return y_pos;
}

unsigned int Atom::get_z_pos(){
	return z_pos;
}

void Atom::set_exists(bool does_exist){
	exists = does_exist;
}

bool Atom::get_exists(){
	if(exists) return true;
	else return false;
}

void Atom::set_grav_pot_energy(){
	grav_pot_energy = mass * get_z_pos();
}

unsigned int Atom::get_grav_pot_energy(){
	return grav_pot_energy;
}

void Atom::set_attempted(bool hasattempted){
	attempted=hasattempted;
}

bool Atom::get_attempted(){
	return attempted;
}

void Atom::set_strength(double val){
	strength=val;
}

double Atom::get_strength(){
	return strength;
}

void Atom::set_index(int i){index=i;}

int Atom::get_index(){return index;}

void Atom::set_type_num(int t){atom_type_num=t;}

int Atom::get_type_num(){return atom_type_num;}




