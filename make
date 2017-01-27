g++ -Wall -fPIC -O2 -c libIsing.c
g++ -shared -o libIsing.so libIsing.o
rm libIsing.o
