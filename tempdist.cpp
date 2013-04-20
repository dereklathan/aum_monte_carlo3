#include <iostream>
#include <fstream>
#include <cstdlib>
#include <string>
#include <cstring>
using namespace std;
int main(int argc, char *argv[]){
	ifstream data(argv[1], ifstream::in);
	double avg;
	double v1,v2,dumb;
	data >> dumb;
	data >> v1;
	data >> dumb;
	data >> v2;
	avg=v2-v1;
	while(!data.eof()){
		v1=v2;
		data >> dumb;
		data >> v2;
		avg+=(v2-v1);
		avg/=2;
	}
	data.close();
	cout << avg*2 << endl;
	return 0;
}

