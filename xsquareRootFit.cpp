#include <iostream>
#include <fstream>
#include <cmath>
#include <string>
#include <cstring>
#include <cstdlib>
using namespace std;

int main(int argc, char *argv[]){
	if(argc!=3){
		cout << "USAGE: ./xsquareRootFit FILENAME.dat TIMESTEPS" << endl;
		return 1;
	}
	ifstream data;
	data.open(argv[1], ifstream::in);
	double val1;
        double val2;
        double coef;
	double num=0;
        double denom=0;
	char * pEnd;
	int timesteps = atoi(argv[2]);
	for(int c=0;c<timesteps;c++){
		data >> val1;
		data >> val2;
		num += pow(val1,0.5)*val2;
                denom += pow(val1,1.0);
               
	}
        coef = num/denom;
	data.close();
	cout << coef << endl;
	return 0;
}
