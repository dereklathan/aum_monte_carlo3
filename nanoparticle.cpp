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
