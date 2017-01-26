g++ -Wall -fPIC -O2 -c libIsing.cc
g++ -shared -o libIsing.so libIsing.o
rm libIsing.o
