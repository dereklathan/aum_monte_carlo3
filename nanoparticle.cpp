#include "nanoparticle.h"
using namespace std;

	Nanoparticle::Nanoparticle(){
		exists=false;
	}
	void Nanoparticle::set_exists(bool does_exist){
		exists=does_exist;
	}
	bool Nanoparticle::get_exists(){
		return exists;
	}
	void Nanoparticle::set_T(double t){
		T=t;
	}
	double Nanoparticle::get_T(){
		return T;
	}
	void Nanoparticle::set_x_pos(int x){
		x_pos=x;
	}
	void Nanoparticle::set_y_pos(int y){
		y_pos=y;
	}
	void Nanoparticle::set_z_pos(int z){
		z_pos=z;
	}
	int Nanoparticle::get_x_pos(){
		return x_pos;
	}
	int Nanoparticle::get_y_pos(){
		return y_pos;
	}
	int Nanoparticle::get_z_pos(){
		return z_pos;
	}
