
lib = -fopenmp

all:
	make task1
	make task2
	make task3
	
task1:
	gcc -o $@ task1.c $(lib) 
task2:
	gcc -o  $@ task2.c -lm $(lib) 
task3:
	g++ -o $@ task3_1.cpp $(lib)
	g++ -o task3_guided task3.cpp $(lib) -DTYPEGUIDED
	g++ -o task3_static task3.cpp $(lib) -DTYPESTATIC
	g++ -o task3_dynamic task3.cpp $(lib) -DTYPEDYNAMIC

remove: all
	rm task1 task2 task3 task3_dynamic task3_guided task3_static
