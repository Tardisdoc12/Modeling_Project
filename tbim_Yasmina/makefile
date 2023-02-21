%.o: %.cpp
	g++ -o $@ -c $<

model:tbim_main.o nano.o
	g++ nano.o tbim_main.o -o tbim_main
	rm *.o
	./tbim_main
	rm tbim_main
