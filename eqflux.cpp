#include <iostream>
#include <fstream>
#include <cstdlib>
#include <string>
using namespace std;

int main(int argc, char *argv[]){
	if(argc!=4){
		cout << "USAGE: ./eqflux FLUXIN.dat FLUXOUT.dat starting_timestep" << endl;
		return 1;
	}
	ifstream fluxin, fluxout;
	fluxin.open(argv[1], ifstream::in);
	fluxout.open(argv[2], ifstream::in);
	int start=atoi(argv[3]);
	int count=0;
	double inval, outval, dval, inavg, outavg;
	string line;
	inavg=0;
	outavg=0;
	for(int c=0;c<start;c++){
		getline(fluxin,line);
		getline(fluxout,line);
	}
	while(!fluxin.eof() && !fluxout.eof()){
		fluxin >> dval;
		fluxout >> dval;
		fluxin >> inval;
		fluxout >> outval;
		inavg+=inval;
		outavg+=outval;
		count++;
	}
	inavg/=count;
	outavg/=count;
	cout << (inavg+outavg)/2 << endl;
		
}
