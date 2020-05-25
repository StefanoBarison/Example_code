CC = g++
CFLAGS = -Wall -O3 --std=c++11

vmc.exe : vmc.o random.o 
	$(CC) random.o vmc.o -o vmc.exe
random.o : random.cpp random.h
	$(CC) -c random.cpp -o random.o $(CFLAGS)
vmc.o : vmc.cpp vmc.h
	$(CC) -c vmc.cpp -o vmc.o $(CFLAGS)
clean :
	rm *.o main.exe seed.out
