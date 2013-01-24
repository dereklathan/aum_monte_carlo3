#include <iostream>
#include <fstream>
#include <cmath>
#include <string>
#include <cstring>
#include <cstdlib>
using namespace std;

int main(int argc, char *argv[]){
	if(argc!=5){
		cout << "USAGE: ./area_between_curves FILENAME.dat TIMESTEPS CONSTANT POWER" << endl;
		return 1;
	}
	ifstream data;
	data.open(argv[1], ifstream::in);
	double val;
	double area=0;
	char * pEnd;
	int timesteps = atoi(argv[2]);
	for(int c=0;c<timesteps;c++){
		data >> val;
		data >> val;
		area += abs(strtod(argv[3],NULL)*pow((double)c+1,strtod(argv[4],NULL))-val);
	}
	data.close();
	cout << area << endl;
	return 0;
}
