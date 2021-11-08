# makefile

version1: version1.o
	mpicc -std=c99 -lm version1.o -o  matrixMain
	mpirun -n 2 ./matrixMain

version1.o: version1.c matrixFunctions.c
	mpicc -std=c99 -c  version1.c

version2: version2.o
	mpicc -std=c99 version2.o -o version2
	mpirun -n 3 ./version2

version2.o: version2.c matrixFunctions.c
	mpicc -std=c99 -c version2.c

timeVersion1: timeVersion1.o
	mpicc -std=c99 timeVersion1.o -o timeVersion1
	# mpirun -n 2 ./timeVersion1 5 5 10
	mpirun -n 2 ./timeVersion1 13 13 10
timeVersion1.o: timeVersion1.c matrixFunctions.c
	mpicc -std=c99 -c timeVersion1.c

timeVersion2: timeVersion2.o
	mpicc -std=c99 timeVersion2.o -o timeVersion2
	mpirun -n 2 ./timeVersion2 5 5 10

timeVersion2.o: timeVersion2.c matrixFunctions.c
	mpicc -std=c99 -c timeVersion2.c
hex:
	hexdump -v -e '25/8 "%2f,"' -e '"\n"' matData.txt

clean:
	rm *.o
	rm time