RCASimulator: RCASimulator.cpp
	g++ -std=c++11 -I/opt/X11/include -L/opt/X11/lib/ -lX11 -lm -lpthread -lpng -lgmpxx -lgmp -O2 -o RCASimulator RCASimulator.cpp 
