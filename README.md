# Local-Strain-Calculation

Find local strain in 3dimension sheared configurations : For details see PRE.88.062206(2013) for how to calculate local strain
   After calculation of local strain we will choose the cutoff of ed so that particle beyond ed_cut are 
   active particle for plastic rearrangement. Then we compute the cluster size distribution of those active particles.
   

   @Compilation :   gcc ed_cluster.c -L/Data/.gsl-2.2.1/lib/ -lgsl -lgslcblas -lm -o ed_cluster.out
   
   @Input  : i) Program needs two configurations file of avalanche, here those are: pxyz-*.restart-0 and pxyz-*.restart-1
             ii)During local strain calculation we will tessilation of the system to get the local volume. For that we use 
	    ./tessilate_bmlj_32K executable with filename-0 as a feed to obtain the tessilation file written as tessilate.

	    ./tessilate_bmlj_32K this executable is obtained from the fortran code 3Dtessilation_srikanth_bmlj.f (simple gfortran compilation), to change the system size one needs to change the N=32000 with proper size inside the code at a few places.

	    iii) N(number of particle).


   @Output  :  We choose five different ed_cut here, and corresponding cluster size distribution files will be generated with proper extension: both in bin and raw data.
   The distribution of the ed will also be saved with proper filename.
	    
