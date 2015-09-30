#define solvation_min  0.000
#define solvation_max  0.300

#define bin 50 //for distribution

#define nstxout 10.0 // nstxout from the mdp file
#define dt 0.001  // dt from the mdp file
#define factor dt * nstxout 
#define ps 2.00  // relaxation desired in time (ps)
//#define monomers 20 // No. of monomers in the chain. Need for the normalization
#define monomers 1 

//#define stay 30 * factor   // Number of frames a water molecule should stay.
//Stay of water molecules. Determined from the residence time

#define SIZE 7784800 		// < 0.62 nm
#define SIZE_2 4308800  	// < 0.62 nm
