work : aum_monte_carlo3.o cube.o Atom.o infile_reader.o vtf_file_writer.o dat_file_writer.o
	g++ -O3 aum_monte_carlo3.o cube.o Atom.o infile_reader.o vtf_file_writer.o dat_file_writer.o -o work

aum_monte_carlo_3.o : aum_monte_carlo3.cpp
	g++ -O3 -c aum_monte_carlo3.cpp

cube.o : cube.cpp cube.h
	g++ -O3 -c cube.cpp

Atom.o : Atom.cpp Atom.h
	g++ -O3 -c Atom.cpp

infile_reader.o : infile_reader.cpp infile_reader.h
	g++ -O3 -c infile_reader.cpp

vtf_file_writer.o : vtf_file_writer.cpp vtf_file_writer.h
	g++ -O3 -c vtf_file_writer.cpp

dat_file_writer.o : dat_file_writer.cpp dat_file_writer.h
	g++ -O3 -c dat_file_writer.cpp

clean : 
	rm *.o
