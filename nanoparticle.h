#ifndef NANOPARTICLE_H
#define NANOPARTICLE_H
using namespace std;

class Nanoparticle{
	private:
		bool exists;
		double T;
	public:
		Nanoparticle();
		//constructor
		void set_exists(bool);
		bool get_exists();
		void set_T(double);
		double get_T();				
};
#endif
