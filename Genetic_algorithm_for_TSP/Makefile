CC = g++
CFLAGS = -Wall -O3 --std=c++11

main.exe : main.o tsp.o 
	$(CC) main.o tsp.o -o main.exe
main.o : main.cpp 
	$(CC) -c main.cpp -o main.o $(CFLAGS)
tsp.o : tsp.cpp tsp.h
	$(CC) -c tsp.cpp -o tsp.o $(CFLAGS)
clean :
	rm *.o main.exe 
