Final Project Part 2
Connor Tluck

This program is meant to compute the potential energy difference between two cubes that are seperated however
are experiencing a changing charge desnity between them. The program uses two differnt integration techniques to do 
this then keeps track of the difference in time between each them of them. The integration takes place in 6 dimensions 
based off the dimensions of the cube. The integration was done using the vegas method from the GNU library as well as 
a monte carlo method of integration. This base code began early on in class during the compuation of pi. The potential 
difference in energy between the two cubes is found from the dipole approximation. 

Upon examination of the results that are printed it seems that the monte carlo integration runs about 20 times slower than the 
Vegas method. The accuracy also is far worse needing many more iterations. 
