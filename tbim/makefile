%.o: %.cpp
	g++ -o $@ -c $<

model:tbim_main.o
	g++ tbim_main.o -o tbim_main
	rm *.o
	./tbim_main
	rm tbim_main
