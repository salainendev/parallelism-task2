all:
	gcc -o program task1.c -fopenmp -DFIRSTVAR -DTYPEGUIDED
	gcc -o program2 task2.c -fopenmp -DFIRSTVAR -DTYPEGUIDED
	g++ -o program3_1var_guided task3.cpp -fopenmp -DFIRSTVAR -DTYPEGUIDED
	g++ -o program3_1var_static task3.cpp -fopenmp -DFIRSTVAR -DTYPESTATIC
	g++ -o program3_1var_dynamic task3.cpp -fopenmp -DFIRSTVAR -DTYPEDYNAMIC
	g++ -o program3 task1.cpp -fopenmp -DSECONDVAR
task1:
	gcc -o program task1.c -fopenmp -DFIRSTVAR -DTYPEGUIDED
task2:
	gcc -o program2 task2.c -fopenmp -DFIRSTVAR -DTYPEGUIDED
task3:
	g++ -o program3_1var_guided task3.cpp -fopenmp -DFIRSTVAR -DTYPEGUIDED
	g++ -o program3_1var_static task3.cpp -fopenmp -DFIRSTVAR -DTYPESTATIC
	g++ -o program3_1var_dynamic task3.cpp -fopenmp -DFIRSTVAR -DTYPEDYNAMIC
remove:
	rm program
	rm program2
	rm program3
	rm program3_1var_dynamic
	rm program3_1var_guided
	rm program3_1var_static
