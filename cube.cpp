#include "cube.h"
using namespace std;

Cube::Cube(){
	srand(time(NULL));
	population=0;
	fixedcount=0;
	flux_in=0;
	flux_out=0;
}

void Cube::set_domain(int x, int y, int z){
	domain_x=x;
	domain_y=y;
	domain_z=z;
	atomlocation = new Atom**[x];
	nparticlelocation = new Nanoparticle**[x];
	for(int c=0;c<x;c++){
		atomlocation[c] = new Atom*[y];
		nparticlelocation[c] = new Nanoparticle*[y];
	}
	for(int c=0;c<x;c++){
		for(int d=0;d<=y;d++){
			atomlocation[c][d] = new Atom[z];
			nparticlelocation[c][d] = new Nanoparticle[z];
		}
	}
}

void Cube::clear(){
	population=0;
	fixedcount=0;
	for(int c=0;c<domain_x;c++){
		for(int d=0;d<domain_y;d++){
			for(int e=0;e<domain_z;e++){
				atomlocation[c][d][e].set_exists(false);
				nparticlelocation[c][d][e].set_exists(false);
			}
		}
	}
}

int Cube::get_domain_x(){
	return domain_x;
}

int Cube::get_domain_y(){
	return domain_y;
}

int Cube::get_domain_z(){
	return domain_z;
}

int Cube::get_volume(){
	return (domain_x*domain_y*domain_z);
}

Atom Cube::insert_atom(Atom &atom){
	int xval, yval, zval;
	xval=rand()%get_domain_x();
	yval=rand()%get_domain_y();
	zval=rand()%get_domain_z();
	while(atomlocation[xval][yval][zval].get_exists()){
		xval=rand()%get_domain_x();
		yval=rand()%get_domain_y();
		zval=rand()%get_domain_z();
	}
	atom.set_exists(true);
	atom.set_x_pos(xval);
	atom.set_y_pos(yval);
	atom.set_z_pos(zval);
	if(atom.is_fixed()){
		atom.set_attempted(true);
		fixedcount++;
	}
	else
		atom.set_attempted(false);
	atom.set_index(population);
	population++;
	atomlocation[xval][yval][zval]=atom;
	return atom;
}

Nanoparticle Cube::insert_nanoparticle(Nanoparticle &nparticle){
	int xval, yval, zval;
	xval=rand()%get_domain_x();
	yval=rand()%get_domain_y();
	zval=rand()%get_domain_z();
	while(nparticlelocation[xval][yval][zval].get_exists()){
		xval=rand()%get_domain_x();
		yval=rand()%get_domain_y();
		zval=rand()%get_domain_z();
	}
	nparticle.set_x_pos(xval);
	nparticle.set_y_pos(yval);
	nparticle.set_z_pos(zval);
	nparticle.set_exists(true);
	nparticlelocation[xval][yval][zval]=nparticle;
	return nparticle;
}

bool Cube::insert_atom(Atom &atom, int x, int y, int z){
	if(atomlocation[x][y][z].get_exists())
		return false;
	else{
		atom.set_exists(true);
		atom.set_x_pos(x);
		atom.set_y_pos(y);
		atom.set_z_pos(z);
		if(atom.is_fixed())
			atom.set_attempted(true);
		else
			atom.set_attempted(false);
		atom.set_index(population);
		population++;
		atomlocation[x][y][z]=atom;
		return true;
	}
}

Atom Cube::get_atom(int x, int y, int z){
		return atomlocation[x][y][z];
}

void Cube::set_occupy_space(int x, int y, int z, bool occupied){
	atomlocation[x][y][z].set_exists(occupied);
}

bool Cube::get_occupy_space(int x, int y, int z){
	if(atomlocation[x][y][z].get_exists())
		return true;
	else
		return false;
}

double Cube::get_population(){return population;}

bool Cube::all_attempted(){
	if(left_to_attempt==0)
		return true;
	else
		return false;

}

void Cube::advance_timestep_pbc(){
	int x1, y1, z1, x2, y2, z2;
	int rand_axis, rand_dir, dir;
	left_to_attempt=population-fixedcount;
	vector<Atom> unattempted;
	for(int c=0;c<domain_x;c++){
		for(int d=0;d<domain_y;d++){
			for(int e=0;e<domain_z;e++){
				if(atomlocation[c][d][e].get_exists() && !atomlocation[c][d][e].is_fixed())
					unattempted.push_back(atomlocation[c][d][e]);
			}
		}
	}		
					
	while(!all_attempted()){
		int index = rand()%unattempted.size();
		x1=unattempted[index].get_x_pos();
		y1=unattempted[index].get_y_pos();
		z1=unattempted[index].get_z_pos();
		unattempted[index]=unattempted[unattempted.size()-1];
		unattempted.pop_back();
		rand_axis=rand()%3;
		rand_dir=rand()%2;
		x2=x1;
		y2=y1;
		z2=z1;
		if(rand_dir==0)	dir=1;
		else dir=-1;
		if(rand_axis==0) x2+=dir;
		else if(rand_axis==1) y2+=dir;
		else z2+=dir;
		if(x2<0) x2=domain_x-1;
		else if(x2>=domain_x) x2=0;
		else if(y2<0) y2=domain_y-1;
		else if(y2>=domain_y) y2=0;
		else if(z2<0) z2=domain_z-1;
		else if(z2>=domain_z) z2=0;
		left_to_attempt--;
		if(!atomlocation[x2][y2][z2].get_exists() && ((calculate_pot_energy_pbc(x1,y1,z1,x2,y2,z2)<=calculate_pot_energy_pbc()) || can_still_move())){
		
				atomlocation[x2][y2][z2]=atomlocation[x1][y1][z1];
				atomlocation[x1][y1][z1].set_exists(false);
				atomlocation[x2][y2][z2].set_exists(true);
				atomlocation[x2][y2][z2].set_x_pos(x2);
				atomlocation[x2][y2][z2].set_y_pos(y2);
				atomlocation[x2][y2][z2].set_z_pos(z2);
		}			
	}

}
void Cube::advance_timestep_opentop(){
	int movecount=0;
	flux_in=0;
	flux_out=0;
	int x1, x2, y1, y2, z1, z2, index;
	int rand_axis, rand_dir, dir;
	vector<Atom> unattempted;
	for(int c=0;c<domain_x;c++){
		for(int d=0;d<domain_y;d++){
			for(int e=0;e<domain_z;e++){
				if(atomlocation[c][d][e].get_exists() && !atomlocation[c][d][e].is_fixed())
					unattempted.push_back(atomlocation[c][d][e]);
			}
		}
	}
	while(unattempted.size()!=0){
		index=rand()%unattempted.size();
		x1=unattempted[index].get_x_pos();
		x2=x1;
		y1=unattempted[index].get_y_pos();
		y2=y1;
		z1=unattempted[index].get_z_pos();
		z2=z1;
		rand_axis = rand()%3;
		rand_dir = rand()%2;
		if(rand_dir==0) dir=1;
		else dir=-1;
		if(rand_axis==0) x2+=dir;
		else if(rand_axis==1) y2+=dir;
		else z2+=dir;
		if(x2>=domain_x) x2=0;
		else if(x2<0) x2=domain_x-1;
		if(y2>=domain_y) y2=0;
		else if(y2<0) y2=domain_y-1;
		if(z2>=0 && !atomlocation[x2][y2][z2].get_exists()){
			if(nparticlelocation[x1][y1][z1].get_exists()){
				T=nparticlelocation[x1][y1][z1].get_T();
			}
			else if(nparticlelocation[x2][y2][z2].get_exists())
				T=nparticlelocation[x2][y2][z2].get_T();
			else
				T=temperature;
			if(/*calculate_pot_energy_opentop(x1,y1,z1,x2,y2,z2)<=calculate_pot_energy_opentop()*/calculate_dE_opentop(x1,y1,z1,x2,y2,z2)<=0 || can_still_move()){
				movecount++;
				if(z2>=domain_z){
					atomlocation[x1][y1][z1].set_exists(false);
					population--;
					flux_out++;
				}
				else if(z1==0 && z2==1 && atomlocation[x1][y1][z1].get_replaceable()){
					atomlocation[x2][y2][z2]=atomlocation[x1][y1][z1];
					atomlocation[x2][y2][z2].set_replaceable(false);
					atomlocation[x2][y2][z2].set_exists(true);
					atomlocation[x2][y2][z2].set_x_pos(x2);
					atomlocation[x2][y2][z2].set_y_pos(y2);
					atomlocation[x2][y2][z2].set_z_pos(z2);
					atomlocation[x1][y1][z1].set_exists(false);
					flux_in++;
					insert_atom(atomlocation[x1][y1][z1], x1, y1, z1);
					atomlocation[x1][y1][z1].set_replaceable(true);
					atomlocation[x1][y1][z1].set_exists(true);
				}
				else{
					atomlocation[x2][y2][z2]=atomlocation[x1][y1][z1];
					atomlocation[x2][y2][z2].set_exists(true);
					if(z1==0 && z2==0 && atomlocation[x1][y1][z1].get_replaceable())
						atomlocation[x2][y2][z2].set_replaceable(true);
					else
						atomlocation[x2][y2][z2].set_replaceable(false);
					atomlocation[x1][y1][z1].set_exists(false);
					atomlocation[x2][y2][z2].set_x_pos(x2);
					atomlocation[x2][y2][z2].set_y_pos(y2);
					atomlocation[x2][y2][z2].set_z_pos(z2);
				}
			}	
		}
		unattempted[index]=unattempted[unattempted.size()-1];
		unattempted.pop_back();
	}
	move_nanoparticles();
	cout << movecount << " have moved\n";
}

double ** Cube::set_interaction_factor(int types){
	interaction_factor = new double*[types];
	for(int c=0;c<types;c++)
		interaction_factor[c] = new double[types];
	for(int c=0;c<types;c++){
		for(int d=0;d<types;d++)
			cin >> interaction_factor[c][d];
	}
	return interaction_factor;
}

void Cube::set_interaction_factor(double ** i_factor){
	interaction_factor=i_factor;
}

double Cube::calculate_pot_energy_pbc(){
	double grav_energy=0;
	double internal_energy=0;
	for(int c=0;c<domain_x;c++){
		for(int d=0;d<domain_y;d++){
			for(int e=0;e<domain_z;e++){
				if(atomlocation[c][d][e].get_exists()){
					grav_energy+=(e*atomlocation[c][d][e].get_mass());
					if(c==0){
						if(atomlocation[domain_x-1][d][e].get_exists())
							internal_energy+=(atomlocation[c][d][e].get_strength()*atomlocation[domain_x-1][d][e].get_strength()*interaction_factor[atomlocation[c][d][e].get_type_num()][atomlocation[domain_x-1][d][e].get_type_num()]);
						if(atomlocation[1][d][e].get_exists())
							internal_energy+=(atomlocation[c][d][e].get_strength()*atomlocation[1][d][e].get_strength()*interaction_factor[atomlocation[c][d][e].get_type_num()][atomlocation[1][d][e].get_type_num()]);					
					}
					else if(c==domain_x-1){
						if(atomlocation[0][d][e].get_exists())
							internal_energy+=(atomlocation[c][d][e].get_strength()*atomlocation[0][d][e].get_strength()*interaction_factor[atomlocation[c][d][e].get_type_num()][atomlocation[0][d][e].get_type_num()]);
						if(atomlocation[domain_x-2][d][e].get_exists())
							internal_energy+=(atomlocation[c][d][e].get_strength()*atomlocation[domain_x-2][d][e].get_strength()*interaction_factor[atomlocation[c][d][e].get_type_num()][atomlocation[domain_x-2][d][e].get_type_num()]);
					}
					else{
						if(atomlocation[c+1][d][e].get_exists())
							internal_energy+=(atomlocation[c][d][e].get_strength()*atomlocation[c+1][d][e].get_strength()*interaction_factor[atomlocation[c][d][e].get_type_num()][atomlocation[c+1][d][e].get_type_num()]);
						if(atomlocation[c-1][d][e].get_exists())
							internal_energy+=(atomlocation[c][d][e].get_strength()*atomlocation[c-1][d][e].get_strength()*interaction_factor[atomlocation[c][d][e].get_type_num()][atomlocation[c-1][d][e].get_type_num()]);
					}
					if(d==0){
						if(atomlocation[c][domain_y-1][e].get_exists())
							internal_energy+=(atomlocation[c][d][e].get_strength()*atomlocation[c][domain_y-1][e].get_strength()*interaction_factor[atomlocation[c][d][e].get_type_num()][atomlocation[c][domain_y-1][e].get_type_num()]);
						if(atomlocation[c][1][e].get_exists())
							internal_energy+=(atomlocation[c][d][e].get_strength()*atomlocation[c][1][e].get_strength()*interaction_factor[atomlocation[c][d][e].get_type_num()][atomlocation[c][1][e].get_type_num()]);
					}
					else if(d==domain_y-1){
						if(atomlocation[c][0][e].get_exists())
							internal_energy+=(atomlocation[c][d][e].get_strength()*atomlocation[c][0][e].get_strength()*interaction_factor[atomlocation[c][d][e].get_type_num()][atomlocation[c][0][e].get_type_num()]);
						if(atomlocation[c][domain_y-2][e].get_exists())
							internal_energy+=(atomlocation[c][d][e].get_strength()*atomlocation[c][domain_y-2][e].get_strength()*interaction_factor[atomlocation[c][d][e].get_type_num()][atomlocation[c][domain_y-2][e].get_type_num()]);
					}
					else{
						if(atomlocation[c][d+1][e].get_exists())
							internal_energy+=(atomlocation[c][d][e].get_strength()*atomlocation[c][d+1][e].get_strength()*interaction_factor[atomlocation[c][d][e].get_type_num()][atomlocation[c][d+1][e].get_type_num()]);
						if(atomlocation[c][d-1][e].get_exists())
							internal_energy+=(atomlocation[c][d][e].get_strength()*atomlocation[c][d-1][e].get_strength()*interaction_factor[atomlocation[c][d][e].get_type_num()][atomlocation[c][d-1][e].get_type_num()]);
					}
					if(e==0){
						if(atomlocation[c][d][domain_z-1].get_exists())
							internal_energy+=(atomlocation[c][d][e].get_strength()*atomlocation[c][d][domain_z-1].get_strength()*interaction_factor[atomlocation[c][d][e].get_type_num()][atomlocation[c][d][domain_z-1].get_type_num()]);
						if(atomlocation[c][d][1].get_exists())
							internal_energy+=(atomlocation[c][d][e].get_strength()*atomlocation[c][d][1].get_strength()*interaction_factor[atomlocation[c][d][e].get_type_num()][atomlocation[c][d][1].get_type_num()]);
					}
					else if(e==domain_z-1){
						if(atomlocation[c][d][0].get_exists())
							internal_energy+=(atomlocation[c][d][e].get_strength()*atomlocation[c][d][0].get_strength()*interaction_factor[atomlocation[c][d][e].get_type_num()][atomlocation[c][d][0].get_type_num()]);
						if(atomlocation[c][d][domain_z-2].get_exists())
							internal_energy+=(atomlocation[c][d][e].get_strength()*atomlocation[c][d][domain_z-2].get_strength()*interaction_factor[atomlocation[c][d][e].get_type_num()][atomlocation[c][d][domain_z-2].get_type_num()]);
					}
					else{
						if(atomlocation[c][d][e+1].get_exists())
							internal_energy+=(atomlocation[c][d][e].get_strength()*atomlocation[c][d][e+1].get_strength()*interaction_factor[atomlocation[c][d][e].get_type_num()][atomlocation[c][d][e+1].get_type_num()]);
						if(atomlocation[c][d][e-1].get_exists())
							internal_energy+=(atomlocation[c][d][e].get_strength()*atomlocation[c][d][e-1].get_strength()*interaction_factor[atomlocation[c][d][e].get_type_num()][atomlocation[c][d][e-1].get_type_num()]);
					}
				}
			}
		}
	}
	return grav_energy + (internal_energy/2);
}

double Cube::calculate_pot_energy_pbc(int x1, int y1, int z1, int x2, int y2, int z2){
	double grav_energy=0;
	double internal_energy=0;
	Atom temp[domain_x][domain_y][domain_z];
	for(int c=0;c<domain_x;c++){
		for(int d=0;d<domain_y;d++){
			for(int e=0;e<domain_z;e++)
				temp[c][d][e]=atomlocation[c][d][e];
		}
	}
	if(true/*x1>=0 && x1<domain_x && x2>=0 && x2<domain_x && y1>=0 && y1<domain_y && y2>=0 && y2<domain_y && z1>=0 && z1<domain_z && z2>=0 && z2<domain_z && temp[x1][y1][z1].get_exists() && !temp[x2][y2][z2].get_exists() && !temp[x2][y2][z2].is_fixed()*/){
		temp[x2][y2][z2]=temp[x1][y1][z1];
		temp[x2][y2][z2].set_x_pos(x2);
		temp[x2][y2][z2].set_y_pos(y2);
		temp[x2][y2][z2].set_z_pos(z2);
		temp[x1][y1][z1].set_exists(false);
		for(int c=0;c<domain_x;c++){
			for(int d=0;d<domain_y;d++){
				for(int e=0;e<domain_z;e++){
					if(temp[c][d][e].get_exists()){
						grav_energy+=(e*temp[c][d][e].get_mass());
						if(c==0){
							if(temp[domain_x-1][d][e].get_exists())
								internal_energy+=(temp[c][d][e].get_strength()*temp[domain_x-1][d][e].get_strength()*interaction_factor[temp[c][d][e].get_type_num()][temp[domain_x-1][d][e].get_type_num()]);
							if(temp[1][d][e].get_exists())
								internal_energy+=(temp[c][d][e].get_strength()*temp[1][d][e].get_strength()*interaction_factor[temp[c][d][e].get_type_num()][temp[1][d][e].get_type_num()]);					
						}
						else if(c==domain_x-1){
							if(temp[0][d][e].get_exists())
								internal_energy+=(temp[c][d][e].get_strength()*temp[0][d][e].get_strength()*interaction_factor[temp[c][d][e].get_type_num()][temp[0][d][e].get_type_num()]);
							if(temp[domain_x-2][d][e].get_exists())
								internal_energy+=(temp[c][d][e].get_strength()*temp[domain_x-2][d][e].get_strength()*interaction_factor[temp[c][d][e].get_type_num()][temp[domain_x-2][d][e].get_type_num()]);
						}
						else{
							if(temp[c+1][d][e].get_exists())
								internal_energy+=(temp[c][d][e].get_strength()*temp[c+1][d][e].get_strength()*interaction_factor[temp[c][d][e].get_type_num()][temp[c+1][d][e].get_type_num()]);
							if(temp[c-1][d][e].get_exists())
								internal_energy+=(temp[c][d][e].get_strength()*temp[c-1][d][e].get_strength()*interaction_factor[temp[c][d][e].get_type_num()][temp[c-1][d][e].get_type_num()]);
						}
						if(d==0){
							if(temp[c][domain_y-1][e].get_exists())
								internal_energy+=(temp[c][d][e].get_strength()*temp[c][domain_y-1][e].get_strength()*interaction_factor[temp[c][d][e].get_type_num()][temp[c][domain_y-1][e].get_type_num()]);
							if(temp[c][1][e].get_exists())
								internal_energy+=(temp[c][d][e].get_strength()*temp[c][1][e].get_strength()*interaction_factor[temp[c][d][e].get_type_num()][temp[c][1][e].get_type_num()]);
						}
						else if(d==domain_y-1){
							if(temp[c][0][e].get_exists())
								internal_energy+=(temp[c][d][e].get_strength()*temp[c][0][e].get_strength()*interaction_factor[temp[c][d][e].get_type_num()][temp[c][0][e].get_type_num()]);
							if(temp[c][domain_y-2][e].get_exists())
								internal_energy+=(temp[c][d][e].get_strength()*temp[c][domain_y-2][e].get_strength()*interaction_factor[temp[c][d][e].get_type_num()][temp[c][domain_y-2][e].get_type_num()]);
						}
						else{
							if(temp[c][d+1][e].get_exists())
								internal_energy+=(temp[c][d][e].get_strength()*temp[c][d+1][e].get_strength()*interaction_factor[temp[c][d][e].get_type_num()][temp[c][d+1][e].get_type_num()]);
							if(temp[c][d-1][e].get_exists())
								internal_energy+=(temp[c][d][e].get_strength()*temp[c][d-1][e].get_strength()*interaction_factor[temp[c][d][e].get_type_num()][temp[c][d-1][e].get_type_num()]);
						}
						if(e==0){
							if(temp[c][d][domain_z-1].get_exists())
								internal_energy+=(temp[c][d][e].get_strength()*temp[c][d][domain_z-1].get_strength()*interaction_factor[temp[c][d][e].get_type_num()][temp[c][d][domain_z-1].get_type_num()]);
							if(temp[c][d][1].get_exists())
								internal_energy+=(temp[c][d][e].get_strength()*temp[c][d][1].get_strength()*interaction_factor[temp[c][d][e].get_type_num()][temp[c][d][1].get_type_num()]);
						}
						else if(e==domain_z-1){
							if(temp[c][d][0].get_exists())
								internal_energy+=(temp[c][d][e].get_strength()*temp[c][d][0].get_strength()*interaction_factor[temp[c][d][e].get_type_num()][temp[c][d][0].get_type_num()]);
							if(temp[c][d][domain_z-2].get_exists())
								internal_energy+=(temp[c][d][e].get_strength()*temp[c][d][domain_z-2].get_strength()*interaction_factor[temp[c][d][e].get_type_num()][temp[c][d][domain_z-2].get_type_num()]);
						}
						else{
							if(temp[c][d][e+1].get_exists())
								internal_energy+=(temp[c][d][e].get_strength()*temp[c][d][e+1].get_strength()*interaction_factor[temp[c][d][e].get_type_num()][temp[c][d][e+1].get_type_num()]);
							if(temp[c][d][e-1].get_exists())
								internal_energy+=(temp[c][d][e].get_strength()*temp[c][d][e-1].get_strength()*interaction_factor[temp[c][d][e].get_type_num()][temp[c][d][e-1].get_type_num()]);
						}
					}
				}
			}
		}
		return grav_energy + (internal_energy/2);
	}
	else
		return 1.79769e+308;	
}

double Cube::calculate_pot_energy_opentop(){
	double grav_energy=0;
	double internal_energy=0;
	for(int c=0;c<domain_x;c++){
		for(int d=0;d<domain_y;d++){
			for(int e=0;e<domain_z;e++){
				if(atomlocation[c][d][e].get_exists()){
					grav_energy+=(e*atomlocation[c][d][e].get_mass());
					if(c==0){
						if(atomlocation[domain_x-1][d][e].get_exists())
							internal_energy+=(atomlocation[c][d][e].get_strength()*atomlocation[domain_x-1][d][e].get_strength()*interaction_factor[atomlocation[c][d][e].get_type_num()][atomlocation[domain_x-1][d][e].get_type_num()]);
						if(atomlocation[1][d][e].get_exists())
							internal_energy+=(atomlocation[c][d][e].get_strength()*atomlocation[1][d][e].get_strength()*interaction_factor[atomlocation[c][d][e].get_type_num()][atomlocation[1][d][e].get_type_num()]);					
					}
					else if(c==domain_x-1){
						if(atomlocation[0][d][e].get_exists())
							internal_energy+=(atomlocation[c][d][e].get_strength()*atomlocation[0][d][e].get_strength()*interaction_factor[atomlocation[c][d][e].get_type_num()][atomlocation[0][d][e].get_type_num()]);
						if(atomlocation[domain_x-2][d][e].get_exists())
							internal_energy+=(atomlocation[c][d][e].get_strength()*atomlocation[domain_x-2][d][e].get_strength()*interaction_factor[atomlocation[c][d][e].get_type_num()][atomlocation[domain_x-2][d][e].get_type_num()]);
					}
					else{
						if(atomlocation[c+1][d][e].get_exists())
							internal_energy+=(atomlocation[c][d][e].get_strength()*atomlocation[c+1][d][e].get_strength()*interaction_factor[atomlocation[c][d][e].get_type_num()][atomlocation[c+1][d][e].get_type_num()]);
						if(atomlocation[c-1][d][e].get_exists())
							internal_energy+=(atomlocation[c][d][e].get_strength()*atomlocation[c-1][d][e].get_strength()*interaction_factor[atomlocation[c][d][e].get_type_num()][atomlocation[c-1][d][e].get_type_num()]);
					}
					if(d==0){
						if(atomlocation[c][domain_y-1][e].get_exists())
							internal_energy+=(atomlocation[c][d][e].get_strength()*atomlocation[c][domain_y-1][e].get_strength()*interaction_factor[atomlocation[c][d][e].get_type_num()][atomlocation[c][domain_y-1][e].get_type_num()]);
						if(atomlocation[c][1][e].get_exists())
							internal_energy+=(atomlocation[c][d][e].get_strength()*atomlocation[c][1][e].get_strength()*interaction_factor[atomlocation[c][d][e].get_type_num()][atomlocation[c][1][e].get_type_num()]);
					}
					else if(d==domain_y-1){
						if(atomlocation[c][0][e].get_exists())
							internal_energy+=(atomlocation[c][d][e].get_strength()*atomlocation[c][0][e].get_strength()*interaction_factor[atomlocation[c][d][e].get_type_num()][atomlocation[c][0][e].get_type_num()]);
						if(atomlocation[c][domain_y-2][e].get_exists())
							internal_energy+=(atomlocation[c][d][e].get_strength()*atomlocation[c][domain_y-2][e].get_strength()*interaction_factor[atomlocation[c][d][e].get_type_num()][atomlocation[c][domain_y-2][e].get_type_num()]);
					}
					else{
						if(atomlocation[c][d+1][e].get_exists())
							internal_energy+=(atomlocation[c][d][e].get_strength()*atomlocation[c][d+1][e].get_strength()*interaction_factor[atomlocation[c][d][e].get_type_num()][atomlocation[c][d+1][e].get_type_num()]);
						if(atomlocation[c][d-1][e].get_exists())
							internal_energy+=(atomlocation[c][d][e].get_strength()*atomlocation[c][d-1][e].get_strength()*interaction_factor[atomlocation[c][d][e].get_type_num()][atomlocation[c][d-1][e].get_type_num()]);
					}
					if(e==0 && atomlocation[c][d][1].get_exists())
						internal_energy+=(atomlocation[c][d][e].get_strength()*atomlocation[c][d][1].get_strength()*interaction_factor[atomlocation[c][d][e].get_type_num()][atomlocation[c][d][1].get_type_num()]);
					else if(e==domain_z-1){
						if(atomlocation[c][d][domain_z-2].get_exists())
						internal_energy+=(atomlocation[c][d][e].get_strength()*atomlocation[c][d][domain_z-2].get_strength()*interaction_factor[atomlocation[c][d][e].get_type_num()][atomlocation[c][d][domain_z-2].get_type_num()]);
					}
					else{
						if(atomlocation[c][d][e+1].get_exists())
							internal_energy+=(atomlocation[c][d][e].get_strength()*atomlocation[c][d][e+1].get_strength()*interaction_factor[atomlocation[c][d][e].get_type_num()][atomlocation[c][d][e+1].get_type_num()]);
						if(atomlocation[c][d][e-1].get_exists())
							internal_energy+=(atomlocation[c][d][e].get_strength()*atomlocation[c][d][e-1].get_strength()*interaction_factor[atomlocation[c][d][e].get_type_num()][atomlocation[c][d][e-1].get_type_num()]);
					}
	
				}
			}
		}
	}
	E1=grav_energy + (internal_energy/2);
	dE = E2-E1;
	return grav_energy + (internal_energy/2);
}

double Cube::calculate_pot_energy_opentop(int x1, int y1, int z1, int x2, int y2, int z2){
	double grav_energy=0;
	double internal_energy=0;
	Atom temp[domain_x][domain_y][domain_z];
	for(int c=0;c<domain_x;c++){
		for(int d=0;d<domain_y;d++){
			for(int e=0;e<domain_z;e++)
				temp[c][d][e]=atomlocation[c][d][e];
		}
	}
	if(z2>=domain_z)
		temp[x1][y1][z1].set_exists(false);
	else if(z2<0 && (x1!=x2 || y1!=y2)){
		temp[x2][y2][z1]=temp[x1][y1][z1];
		temp[x2][y2][z1].set_x_pos(x2);
		temp[x2][y2][z1].set_y_pos(y2);
		temp[x1][y1][z1].set_exists(false);
	}
	else if(x1!=x2 || y1!=y2 || z1!=z2){
		temp[x2][y2][z2]=temp[x1][y1][z1];
		temp[x2][y2][z2].set_x_pos(x2);
		temp[x2][y2][z2].set_y_pos(y2);
		temp[x2][y2][z2].set_z_pos(z2);
		temp[x1][y1][z1].set_exists(false);
	}
	for(int c=0;c<domain_x;c++){
		for(int d=0;d<domain_y;d++){
			for(int e=0;e<domain_z;e++){
				if(temp[c][d][e].get_exists()){
					grav_energy+=(e*temp[c][d][e].get_mass());
					if(c==0){
						if(temp[domain_x-1][d][e].get_exists())
							internal_energy+=(temp[c][d][e].get_strength()*temp[domain_x-1][d][e].get_strength()*interaction_factor[temp[c][d][e].get_type_num()][temp[domain_x-1][d][e].get_type_num()]);
						if(temp[1][d][e].get_exists())
							internal_energy+=(temp[c][d][e].get_strength()*temp[1][d][e].get_strength()*interaction_factor[temp[c][d][e].get_type_num()][temp[1][d][e].get_type_num()]);					
					}
					else if(c==domain_x-1){
						if(temp[0][d][e].get_exists())
							internal_energy+=(temp[c][d][e].get_strength()*temp[0][d][e].get_strength()*interaction_factor[temp[c][d][e].get_type_num()][temp[0][d][e].get_type_num()]);
						if(temp[domain_x-2][d][e].get_exists())
							internal_energy+=(temp[c][d][e].get_strength()*temp[domain_x-2][d][e].get_strength()*interaction_factor[temp[c][d][e].get_type_num()][temp[domain_x-2][d][e].get_type_num()]);
					}
					else{
						if(temp[c+1][d][e].get_exists())
							internal_energy+=(temp[c][d][e].get_strength()*temp[c+1][d][e].get_strength()*interaction_factor[temp[c][d][e].get_type_num()][temp[c+1][d][e].get_type_num()]);
						if(temp[c-1][d][e].get_exists())
							internal_energy+=(temp[c][d][e].get_strength()*temp[c-1][d][e].get_strength()*interaction_factor[temp[c][d][e].get_type_num()][temp[c-1][d][e].get_type_num()]);
					}
					if(d==0){
						if(temp[c][domain_y-1][e].get_exists())
							internal_energy+=(temp[c][d][e].get_strength()*temp[c][domain_y-1][e].get_strength()*interaction_factor[temp[c][d][e].get_type_num()][temp[c][domain_y-1][e].get_type_num()]);
						if(temp[c][1][e].get_exists())
							internal_energy+=(temp[c][d][e].get_strength()*temp[c][1][e].get_strength()*interaction_factor[temp[c][d][e].get_type_num()][temp[c][1][e].get_type_num()]);
					}
					else if(d==domain_y-1){
						if(temp[c][0][e].get_exists())
							internal_energy+=(temp[c][d][e].get_strength()*temp[c][0][e].get_strength()*interaction_factor[temp[c][d][e].get_type_num()][temp[c][0][e].get_type_num()]);
						if(temp[c][domain_y-2][e].get_exists())
							internal_energy+=(temp[c][d][e].get_strength()*temp[c][domain_y-2][e].get_strength()*interaction_factor[temp[c][d][e].get_type_num()][temp[c][domain_y-2][e].get_type_num()]);
					}
					else{
						if(temp[c][d+1][e].get_exists())
							internal_energy+=(temp[c][d][e].get_strength()*temp[c][d+1][e].get_strength()*interaction_factor[temp[c][d][e].get_type_num()][temp[c][d+1][e].get_type_num()]);
						if(temp[c][d-1][e].get_exists())
							internal_energy+=(temp[c][d][e].get_strength()*temp[c][d-1][e].get_strength()*interaction_factor[temp[c][d][e].get_type_num()][temp[c][d-1][e].get_type_num()]);
					}
					if(e==0 && temp[c][d][1].get_exists())
						internal_energy+=(temp[c][d][e].get_strength()*temp[c][d][1].get_strength()*interaction_factor[temp[c][d][e].get_type_num()][temp[c][d][1].get_type_num()]);
					else if(e==domain_z-1){
						if(temp[c][d][domain_z-2].get_exists())
						internal_energy+=(temp[c][d][e].get_strength()*temp[c][d][domain_z-2].get_strength()*interaction_factor[temp[c][d][e].get_type_num()][temp[c][d][domain_z-2].get_type_num()]);
					}
					else{
						if(temp[c][d][e+1].get_exists())
							internal_energy+=(temp[c][d][e].get_strength()*temp[c][d][e+1].get_strength()*interaction_factor[temp[c][d][e].get_type_num()][temp[c][d][e+1].get_type_num()]);
						if(temp[c][d][e-1].get_exists())
							internal_energy+=(temp[c][d][e].get_strength()*temp[c][d][e-1].get_strength()*interaction_factor[temp[c][d][e].get_type_num()][temp[c][d][e-1].get_type_num()]);
					}
	
				}
			}
		}
	}
	E2=grav_energy + (internal_energy/2);
	return grav_energy + (internal_energy/2);
}

double Cube::calculate_dE_opentop(int x1,int y1,int z1,int x2, int y2, int z2){
	double grav_energy=0;
	double internal_energy=0;
	Atom temp[5][5][5];
	int x_index, y_index, dx, dy, dz;
	for(int c=-2;c<3;c++){
		if(x1+c<0)
			x_index=domain_x+x1+c;
		else if(x1+c>domain_x-1)
			x_index=domain_x-1-x1+c;
		else
			x_index=x1+c;

		for(int d=-2;d<3;d++){
			if(y1+d<0)
				y_index=domain_y+y1+d;
			else if(y1+d>domain_y-1)
				y_index=domain_y-1-y1+d;
			else
				y_index=y1+d;
			//cout << x1 << " " << c << " " << x_index << " " << y1 << " " << d << " " << y_index << endl;
			for(int e=-2;e<3;e++){
				if(z1+e<0 || z1+e>domain_z-1)
					temp[c+2][d+2][e+2].set_exists(false);
				else{
					temp[c+2][d+2][e+2]=atomlocation[x_index][y_index][z1+e];
					temp[c+2][d+2][e+2].set_exists(true);
				}
			}
		}
	}
	for(int c=0;c<5;c++){
		for(int d=0;d<5;d++){
			for(int e=0;e<5;e++){
				if(temp[c][d][e].get_exists()){
					grav_energy+=e*temp[c][d][e].get_mass();
					if(c>0 && temp[c-1][d][e].get_exists())
						internal_energy+=(temp[c][d][e].get_strength()*temp[c-1][d][e].get_strength()*interaction_factor[temp[c][d][e].get_type_num()][temp[c-1][d][e].get_type_num()]);
					if(c<4 && temp[c+1][d][e].get_exists())
						internal_energy+=(temp[c][d][e].get_strength()*temp[c+1][d][e].get_strength()*interaction_factor[temp[c][d][e].get_type_num()][temp[c+1][d][e].get_type_num()]);
					if(d>0 && temp[c][d-1][e].get_exists())
						internal_energy+=(temp[c][d][e].get_strength()*temp[c][d-1][e].get_strength()*interaction_factor[temp[c][d][e].get_type_num()][temp[c][d-1][e].get_type_num()]);
					if(d<4 && temp[c][d+1][e].get_exists())
						internal_energy+=(temp[c][d][e].get_strength()*temp[c][d+1][e].get_strength()*interaction_factor[temp[c][d][e].get_type_num()][temp[c][d+1][e].get_type_num()]);
					if(e>0 && temp[c][d][e-1].get_exists())
						internal_energy+=(temp[c][d][e].get_strength()*temp[c][d][e-1].get_strength()*interaction_factor[temp[c][d][e].get_type_num()][temp[c][d][e-1].get_type_num()]);
					if(e<4 && temp[c][d][e+1].get_exists())
						internal_energy+=(temp[c][d][e].get_strength()*temp[c][d][e+1].get_strength()*interaction_factor[temp[c][d][e].get_type_num()][temp[c][d][e+1].get_type_num()]);
				}
			}
		}
	}
	E1=grav_energy+(internal_energy/2);
	grav_energy=0;
	internal_energy=0;
	if(x1==0 && x2==domain_x-1)
		dx=-1;
	else if(x1==domain_x-1 && x2==0)
		dx=1;
	else
		dx=x2-x1;
	if(y1==0 && y2==domain_y-1)
		dy=-1;
	else if(y1==domain_y-1 && y2==0)
		dy=1;
	else
		dy=y2-y1;
	temp[2+dx][2+dy][2+z2-z1]=temp[2][2][2];
	if(z1!=0 && z2!=1)
		temp[2][2][2].set_exists(false);
	for(int c=0;c<5;c++){
		for(int d=0;d<5;d++){
			for(int e=0;e<5;e++){
				if(temp[c][d][e].get_exists()){
					grav_energy+=e*temp[c][d][e].get_mass();
					if(c>0 && temp[c-1][d][e].get_exists())
						internal_energy+=(temp[c][d][e].get_strength()*temp[c-1][d][e].get_strength()*interaction_factor[temp[c][d][e].get_type_num()][temp[c-1][d][e].get_type_num()]);
					if(c<4 && temp[c+1][d][e].get_exists())
						internal_energy+=(temp[c][d][e].get_strength()*temp[c+1][d][e].get_strength()*interaction_factor[temp[c][d][e].get_type_num()][temp[c+1][d][e].get_type_num()]);
					if(d>0 && temp[c][d-1][e].get_exists())
						internal_energy+=(temp[c][d][e].get_strength()*temp[c][d-1][e].get_strength()*interaction_factor[temp[c][d][e].get_type_num()][temp[c][d-1][e].get_type_num()]);
					if(d<4 && temp[c][d+1][e].get_exists())
						internal_energy+=(temp[c][d][e].get_strength()*temp[c][d+1][e].get_strength()*interaction_factor[temp[c][d][e].get_type_num()][temp[c][d+1][e].get_type_num()]);
					if(e>0 && temp[c][d][e-1].get_exists())
						internal_energy+=(temp[c][d][e].get_strength()*temp[c][d][e-1].get_strength()*interaction_factor[temp[c][d][e].get_type_num()][temp[c][d][e-1].get_type_num()]);
					if(e<4 && temp[c][d][e+1].get_exists())
						internal_energy+=(temp[c][d][e].get_strength()*temp[c][d][e+1].get_strength()*interaction_factor[temp[c][d][e].get_type_num()][temp[c][d][e+1].get_type_num()]);
				}
			}
		}
	}
	E2=grav_energy+(internal_energy/2);
	dE=E2-E1;
	return dE;
	
}


void Cube::set_temp(double t){temperature=t;}

double Cube::get_temp(){return temperature;}

bool Cube::can_still_move(){
	double r = (double)rand()/RAND_MAX;
	if(r<=pow((double)2.71828,(double)(0-dE)/T)){
		return true;
	}
	else
		return false;
}

void Cube::seed_random(int n){srand(n);}

bool Cube::obc_eq(){
	if(flux_in==0 || flux_out==0)
		return false;
	else if(fabs(flux_out-flux_in)/flux_in < 0.05){
		cout << flux_in << " " << flux_out << endl;
		cout << "q = " << (flux_out+flux_in)/2 << endl;
		return true;
	}
	else
		return false;
}

double Cube::get_flux_in(){return flux_in;}

double Cube::get_flux_out(){return flux_out;}

void Cube::move_nanoparticles(){
	int axis, dir, dest_x, dest_y, dest_z, index;
	vector<Nanoparticle> unmoved;
	for(int c=0;c<domain_x;c++){
		for(int d=0;d<domain_y;d++){
			for(int e=0;e<domain_z;e++){
				if(nparticlelocation[c][d][e].get_exists())
					unmoved.push_back(nparticlelocation[c][d][e]);
			}
		}
	}
	while(unmoved.size()!=0){
		index=rand()%unmoved.size();
		dest_x=unmoved[index].get_x_pos();
		dest_y=unmoved[index].get_y_pos();
		dest_z=unmoved[index].get_z_pos();
		axis=rand()%3;
		dir=rand()%2;
		if(dir==0)
			dir=-1;
		else
			dir=1;
		if(axis==0){
			dest_x+=dir;
			if(dest_x>=domain_x)
				dest_x=0;
			else if(dest_x<0)
				dest_x=domain_x-1;
		}
		else if(axis==1){
			dest_y+=dir;
			if(dest_y>=domain_y)
				dest_y=0;
			else if(dest_y<0)
				dest_y=domain_y-1;
		}
		else{
			dest_z+=dir;
			if(dest_z>=domain_z)
				dest_z=0;
			else if(dest_z<0)
				dest_z=domain_z-1;
		}
		if(!nparticlelocation[dest_x][dest_y][dest_z].get_exists()){
			nparticlelocation[unmoved[index].get_x_pos()][unmoved[index].get_y_pos()][unmoved[index].get_z_pos()].set_exists(false);
			unmoved[index].set_x_pos(dest_x);
			unmoved[index].set_y_pos(dest_y);
			unmoved[index].set_z_pos(dest_z);
			nparticlelocation[dest_x][dest_y][dest_z]=unmoved[index];
		}
		unmoved[index]=unmoved[unmoved.size()-1];
		unmoved.pop_back();
	}
								
}
