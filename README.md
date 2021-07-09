# Fomite Simulation

Source code (C) for a stochastic simulation of fomite-mediated disease transmission.

For help/details contact andrew.di.battista@ultraleap.com

This code is available as open-source under Apache 2.0 license. (See LICENSE for details)

Please note the code is used for experimental purposes and is not a refined piece of software! Why was it written in C ? Stochastic simulating are computationaly intense; originaly developed using Matlab/Octave scripts, precompiled C code runs ~500 times faster! 



**How to Compile**

_Linux/MacOS_

- The project was originaly built with _gcc_ on a Linux platform; if not already preinstalled, try _sudo apt install build-essential_ to get _gcc_ and _make_. On MacOS, use homebrew _brew isntall gcc_ and _brew istall make_.

- From the command line, navigate to the folder where the source code has been downloaded and extracxted to. Run _make_.

_Windows_ (not tested)

- _gcc_ is also available for Windows e.g. minGW. However, the source code can also be imported into a standard VS project.



**How to run**

_Linux/MacOS_

- run as ./fomite myoutputcsvfile      

_note_ no .csv extension required

On Windows, run _fomite.exe myoutputcsvfile_ (or whatever the name of your executable file ends up being).

The output .csv file holds statisitical results from the simualtion (_R values_, the number of infections at each location etc.).



The default code will is setup to run 10000 realisations of an airport simualtion with 12000 people. The simulation will vary the proportion of touchscreens being replaced by touch-less alternatives.
Note: this may take a while to run on some computers. You can modify the simualtion by changing parameters in _fomite.c_ (see comments in the source code for details or contact me!)