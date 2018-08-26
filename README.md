# Gene regulatory network

Modeling gene regulatory networks requires to take into account the stochastic nature of
biochemical reactions. The presence of slow processes in gene expression, such as 
transcription and translation, reveal the necessity to include delays in the stochastic
algorithms. This introduces theoretical and computational complications because of the 
non-Markovian properties of the system. 

In this project I implemented a stochastic algorithm that includes delays using
[Anderson's algorithm] (https://aip.scitation.org/doi/abs/10.1063/1.2799998). This work
 was part of my M.Sc. in Biophysics thesis, that was supervised by 
Prof. M. Palassini at University of Barcelona.

## Explanation
Cyclical temporal patterns of gene expression in biological systems have attracted a great
interest. Experimental observations have been realized in different systems, like [the 
p53 transcription factor] (http://msb.embopress.org/content/2/1/2006.0033). Here I have 
modeled the p53-Mdm2 gene regulatory network,
whose core consists of a delayed feedback loop, in the following way:
  * Reaction 1: 0 -> P
  * Reaction 2: 0 -> M
  * Reaction 3: P -> P + M (delayed reaction)
  * Reaction 4: P -> 0
  * Reaction 5: M + P -> M
  * Reaction 6: M -> 0
P and M being the number of p53 and Mdm2 proteins, respectively. This reaction set is explicited 
in the input file 'p53\_data.dat' is the input file. It has been structured
in the following way: 1st column contains the reaction number (that matches the description
above), 2nd column corresponds to the number of reactants in the reaction (note that if we
have a generation of one specie, the number here needs to be 1 not 0), 3rd column corresponds
to the number of products in the reaction (analogous as in the previous case, if we have a degradation
of one specie, the number here needs to be 1 not 0), 4th and 5th columns contain the reactants
labels, 6th and 7th columns contain the products labels, and finally columns from 8 to 11
contain the stochiometric matrix.
'constants.dat' contains the reaction rates for each one of the network reactions (in 1/s):
column #1 corresponding to the reaction number, #2 is the reaction rate value, #3 contains a 0 if the 
reaction does not have a delay or 1 otherwise, and #4 is the delay of the reaction (in seconds).
'molecules.dat' contains the initial number of molecules for each specie (p53 and Mdm2 proteins).
'P0\_Hill.dat' contains the parameters of the Hill function that govern reaction 3 (that accounts for 
cooperativity and saturation in Mdm2 production activated by p53): 1st column contains the half-activation
threshold, and 2nd column contains the cooperativity (Hill exponent).

# Installation

* Clone this repo to your computer and cd into `gene-regulatory-network` directory.
* To compile the C code execute `make` (source code has been successfully compiled in 
Mac OS X 10.11.6 64-bit, it might require linking to math libraries when
compiled in Linux and/or adding `<malloc.h>` at the beginning of anderson\_Model2\_Hurdle.c).
* To run the program: `./execution.sh` (it requires ksh and gnuplot).

After execution, the program will create an output file containing the data for the number of
molecules of the two simulated species (p53 and Mdm2) over time (out\_p53.dat) and a plot
for these two time-series (timeseries.pdf).

![alt text](https://raw.githubusercontent.com/mdies/gene-regulatory-network/blob/master/timeseries.png)
