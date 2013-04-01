#ifndef NANOPARTICLE_H
#define NANOPARTICLE_H
using namespace std;

class Nanoparticle{
	private:
		bool exists;
		double T;
		int x_pos, y_pos, z_pos;
		double index;
	public:
		Nanoparticle();
		//constructor
		void set_exists(bool);
		bool get_exists();
		void set_T(double);
		double get_T();
		void set_x_pos(int);
		void set_y_pos(int);
		void set_z_pos(int);
		int get_x_pos();
		int get_y_pos();
		int get_z_pos();
		void set_index(double);
		double get_index();		
};
#endif
