CC = g++
CFLAGS = -Wall -O3 --std=c++11

main.exe : main.o random.o monteC.o MolDyn_NVE.o 
	$(CC) random.o main.o monteC.o MolDyn_NVE.o -o main.exe
main.o : main.cpp
	$(CC) -c main.cpp -o main.o $(CFLAGS)
random.o : random.cpp random.h
	$(CC) -c random.cpp -o random.o $(CFLAGS)
monteC.o : monteC.cpp monteC.h
	$(CC) -c monteC.cpp -o monteC.o $(CFLAGS)
MolDyn_NVE.o : MolDyn_NVE.cpp MolDyn_NVE.h
	$(CC) -c MolDyn_NVE.cpp -o MolDyn_NVE.o $(CFLAGS)	
clean :
	rm *.o main.exe #seed.out
