CC = g++
CFLAGS = -Wall -O3 --std=c++11

main.exe : main.o random.o TSP_sim_ann.o
	mpicxx random.o main.o TSP_sim_ann.o -o main.exe
main.o : main.cpp
	mpicxx -c main.cpp -o main.o $(CFLAGS)
random.o : random.cpp random.h
	$(CC) -c random.cpp -o random.o $(CFLAGS)
TSP_sim_ann.o : TSP_sim_ann.cpp TSP_sim_ann.h
	$(CC) -c TSP_sim_ann.cpp -o TSP_sim_ann.o $(CFLAGS)
clean :
	rm *.o main.exe #seed.out
