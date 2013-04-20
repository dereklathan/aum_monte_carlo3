#include <iostream>
#include <cmath>
#include <string>
#include <cstring>
#include <cstdlib>
using namespace std;

int main(int argc, char *argv[]){
	if(argc!=3){
		cout << "USAGE: ./probcalc T_base T_nanoparticle";
		return 1;
	}
	double T_base=atof(argv[1]);
	double T_np=atof(argv[2]);
	for(int c=1;c<6;c++)
		cout << c << " " << pow(2.71828,-c/T_base) << " " << pow(2.71828,-c/T_np) << " " << pow(2.71828,-c/T_np)/pow(2.71828,-c/T_base) << endl;
	return 0;
}
