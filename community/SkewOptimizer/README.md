Skew Optimizer
==============

Command line util to optimize the Murphy E value by adjusting the skew value. <br />

Contents
--------

* msieve_poly.cpp -- code from msieve.
* skewopt.cpp -- main program for skewopt.

Building on Linux
-----------------
Have the gmp library already installed, plus g++.


g++ -c msieve_poly.cpp <br />
g++ skewopt.cpp msieve_poly.o -lgmp -o skewopt <br />

Running
-------
./skewopt <br />
Gives the command line parameters required <br />
<br />

Example
-------
2*x^6 + 1, x - 118571099379011784113736688648896417641748464297615937576404566024103044751294464 <br />
<br />
./skewopt  -118571099379011784113736688648896417641748464297615937576404566024103044751294464 1 1 0 0 0 0 0 2 0 0 <br />
Best Skew: 1.29779 <br />
