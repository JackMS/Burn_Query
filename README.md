Burn_Query
==========

A HDF5 Query Generation Tool designed for use with the Burn package and SNSPH, using 
USEEPP for interfacing with the HDF5 files. 

Created because I got sick of writing subtly different variations of the same
code to perform queries for a research project. When run with HDF5 files as 
its arguments will gather input from the user at the CLI to construct a query for
the first isotope listed above a mass fraction specified and then search for secondary
isotopes specified when that first is found. The canonical example is to search all
particles for 26Al above E6 abundance and report 28Si,29Si, and 30Si also above E6; 
which with this program would consist of 

./burn_query [files]

1 [add isotope]

13 [nN]

13 [nZ]

6  [10E-n]

1 [add isotope]

14 [nN]

14 [nZ]

6  [10E-n]

1 [add isotope]

15 [nN]

14 [nZ]

6  [10E-n]

1 [add isotope]

16 [nN]

14 [nZ]

6  [10E-n]

4 [run query]

y [to confirm correct query]

file-name [for output file]

...

...

...

Done!

A lot to be done to clean up the process of gathering data (all in one line isotope entries, for one)
and features still missing (such as wild-card mode for finding every isotope for a given element present 
above the cut-off) but it works in the present form "well enough". 

A work in progress, and I'm always welcome to commentary and feedback regarding
features to add or ways to improve the code.

To build the code simply use the makefile as per usual, requires SE to be installed 
(http://www.astro.keele.ac.uk/nugrid/work-packages/io-technologies/hdf5-useepp-format/installing-se) 
and defaults to looking for it in your HOME directory. Since this will run principally in 
remote systems I imagine you shouldn't need to change this, but the Makefile is very simple 
and can be readily modified as needed.

