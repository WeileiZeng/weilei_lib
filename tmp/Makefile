main: libshared.so main.o
	$(CXX) -o main  main.o -L. -lshared

libshared.so: shared.cpp
	$(CXX) -fPIC -c shared.cpp -o shared.o
	$(CXX) -shared  -Wl,-soname,libshared.so -o libshared.so shared.o

clean:
	rm *.o *.so
