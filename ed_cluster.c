/**
   @Program : ed.c   
   @Author  : himangsu
   @Note    : Find local strain in 3dimension sheared configurations : For details see PRE.88.062206(2013) for how to calculate local strain
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
	    
**/
#include<stdio.h>
#include<math.h>
#include<stdlib.h>
#include<string.h>


#define totframe 1000

#define dim 3
#define dim1 (dim+1)
#define dim2 (dim*dim) 
#define dim3 (dim2*dim)
#define dim2_dim1 (dim2+dim)

#define binwidth 0.000001
#define binsize 1000
#define binof 1.2
#define tolerance 0.000001


#define edmax 100   // maximum of ed
#define no_of_edcut 6
//#define epsilon_cutoff 0.15  //epsilon cut off where power-law brakes: see below how we chose different ed_cut and compue the results at a shot
#define rcutoff 1.4        //First cordination shell of SiO2



int solve_linear_eq (double *a_data, double *b_data, double *x_data, long N);
int delta(int i, int j);

int main(int argc, char *argv[])
{     

  FILE *fp,*fp1,*fp2,*fp3,*fplst;

  static char FILE[256],FILE1[256],FILE2[256],FILE3[256];
  static char FILEpath[256],FILElst[256],FILEopen[256];
 
  
  static char filename[128],command[128],charac[10];

  static double L;
  static int N;  
  static long i,j,k,l;
  static double t;
  
  static double charge,x_cor,y_cor,z_cor,tild0,tild1;
  static double xr,yr,zr;
  static int iframe,id,itype,nx,ny,nz;
  static char char_dummy[100];
  static float dummy_float,epsilon_cutoff;

  
  static double Li,Lf,dr;

  static float findex;  
  static int index;  
  static double bincount[binsize],CDFbincount[binsize],tot_count;
  static double epsilon[dim][dim],epsilon_d,Trepsilon;
  static double epsilon_dev[dim][dim], epsilon_sq[dim][dim]; 

  static int atom[dim1], iatom;
  static double a_data[dim2*dim2],b_data[dim2];
  static double x_data[dim2];
  static double a_org[dim2*dim2];


  // list of edcut for which distribution to be evaluated
  static int iedcut;
  static double epscutoffarray[no_of_edcut];
  epscutoffarray[0]=0.12;
  epscutoffarray[1]=0.15;
  epscutoffarray[2]=0.17;
  epscutoffarray[3]=0.20;
  epscutoffarray[4]=0.22;
  epscutoffarray[5]=0.25;
  //Note that array-lengh should be changed accordingly for tot_ns[][length]

  
 
  if(argc>1){
    N=atof(argv[1]);
   }
  else{
    printf("Give N \n");
    return 0;
  }
  
 

  for(i=0;i<binsize;i++){
    CDFbincount[i]=0.0;
    bincount[i]=0.0;
  }

  //--------------------------------------------------------------------
  
  //=================defining various array ============================
  float **xyz_0,**xyz_1;
    xyz_0=malloc(N*sizeof(float *));
    xyz_1=malloc(N*sizeof(float *));
    for(j=0;j<N;j++){
      xyz_0[j]=malloc(3*sizeof(float));
      xyz_1[j]=malloc(3*sizeof(float));
    }
    for(j=0;j<N;j++)
      for(k=0;k<3;k++){
	xyz_0[j][k]=0.0;
	xyz_1[j][k]=0.0;
      }
   

    int *Activeparticle;
    Activeparticle=malloc(N*sizeof(int));

    int **ActiveparticleARRAY;
    ActiveparticleARRAY=malloc(N*sizeof(int *));
    for(j=0;j<N;j++)ActiveparticleARRAY[j]=malloc(no_of_edcut*sizeof(int));
    for(j=0;j<N;j++)for(k=0;k<6;k++)ActiveparticleARRAY[j][k]=0;

    int iplist,*PLIST,*LIST,*ID,*ns,**tot_ns;
    PLIST=malloc(N*sizeof(int));
    LIST =malloc(N*sizeof(int));
    ID   =malloc(N*sizeof(int));
    ns   =malloc((N+1)*sizeof(int));

    tot_ns   =malloc((N+1)*sizeof(int *));
    for(j=0;j<=N;j++)tot_ns[j]=malloc(no_of_edcut*sizeof(int));
    
    int iroot,inew,n,ncl,k1,k2,ntot;
    
    for(i=0;i<=N;i++)for(k=0;k<no_of_edcut;k++)tot_ns[i][k]=0;
    //===================================================================
   
      

      //===================Reading Initial file====================//
      sprintf(FILE,"pxyz-0.0002-aqs-1252.restart-0");
      fp=fopen(FILE,"r");
      printf("%s\n",FILE);      

      while (!feof(fp)){
	fscanf(fp,"%*[^\n]\n");
	fscanf(fp,"%d\n",&N);
	fscanf(fp,"%*[^\n]\n");
	fscanf(fp,"%*[^\n]\n");
	fscanf(fp,"%lf %lf \n",&Li,&Lf);
	fscanf(fp,"%*[^\n]\n");	fscanf(fp,"%*[^\n]\n");	fscanf(fp,"%*[^\n]\n");
	fscanf(fp,"%f %lf %f %s %s %s\n",&dummy_float, &tild0, &dummy_float,&char_dummy,&char_dummy,&char_dummy);
	for(i=0;i<4;i++)fscanf(fp,"%*[^\n]\n");
	for(i=0;i<N;i++){
	  fscanf(fp,"%d %d %lf %lf %lf %d %d %d\n",&id,&itype,&x_cor,&y_cor,&z_cor,&nx,&ny,&nz);
	  id-=1;
	  xyz_0[id][0]=x_cor;
	  xyz_0[id][1]=y_cor;
	  xyz_0[id][2]=z_cor;
	}
	
	for(i=0;i<=N;i++)fscanf(fp,"%*[^\n]\n");
      }
      fclose(fp);
      
      L=Lf-Li;
      

      //===================Reading final file====================//

      sprintf(FILE1,"pxyz-0.0002-aqs-1252.restart-1");
      fp1=fopen(FILE1,"r");
      
      while (!feof(fp1)){
	fscanf(fp1,"%*[^\n]\n");
	fscanf(fp1,"%d\n",&N);
	fscanf(fp1,"%*[^\n]\n");
	fscanf(fp1,"%*[^\n]\n");
	fscanf(fp1,"%lf %lf \n",&Li,&Lf);
	fscanf(fp1,"%*[^\n]\n");	
	fscanf(fp1,"%*[^\n]\n");
	fscanf(fp1,"%*[^\n]\n");
	fscanf(fp1,"%f %lf %f %s %s %s\n",&dummy_float, &tild1, &dummy_float,&char_dummy,&char_dummy,&char_dummy);
	for(i=0;i<4;i++)fscanf(fp1,"%*[^\n]\n");
	for(i=0;i<N;i++){
	  fscanf(fp1,"%d %d %lf %lf %lf %d %d %d\n",&id,&itype,&x_cor,&y_cor,&z_cor,&nx,&ny,&nz);
	  id-=1;
	  xyz_1[id][0]=x_cor;
	  xyz_1[id][1]=y_cor;
	  xyz_1[id][2]=z_cor;
	}
	for(i=0;i<=N;i++)fscanf(fp1,"%*[^\n]\n");
      }
      fclose(fp1);

      //L=Lf-Li;

      for(i=0;i<N;i++){
	Activeparticle[i]=0;
	for(k=0;k<no_of_edcut;k++)ActiveparticleARRAY[i][k]=0;
      }
	

      sprintf(command,"./tessilate_bmlj_32K %s > tlist",FILE);
      system(command);
      
      sprintf(FILE2,"tlist");
      fp2=fopen(FILE2,"r");
      
      
     //---------------------------------------------------------------
      //This part is newly added; For some configurations, tessilation
      //does not work properly, creating a file started with "avcount"
      //or "crmin" string; For such file with following operation ed
      //count neglected.
      fscanf(fp2,"%s",&charac);
      if(atof(charac) == 0){
      	//continue;
	goto here;
      }
      else rewind(fp2);
        while(!feof(fp2)){
	  fscanf(fp2,"%d %d %d %d",&atom[0],&atom[1],&atom[2],&atom[3]);
	  for(i=0;i<dim1;i++)atom[i]--;
	  for(i=0;i<dim2*dim2;i++) 
	    a_data[i]=0.0;
	  for(i=0;i<dim2;i++){
	    b_data[i]=0.0; 
	    x_data[i]=0.0;
	  } 
	  
	  
	  for(i=1;i<4;i++){
	    l=i-1;

	    xr=xyz_0[atom[i]][0]-xyz_0[atom[0]][0];
	    yr=xyz_0[atom[i]][1]-xyz_0[atom[0]][1];
	    zr=xyz_0[atom[i]][2]-xyz_0[atom[0]][2];
	    
	    xr=xr-tild0*round(zr/L);
	    xr=xr-L*round(xr/L);
	    yr=yr-L*round(yr/L);
	    zr=zr-L*round(zr/L);
	      	      
	      //assign a matrix
	      a_data[dim3*l]=xr;
	      a_data[dim3*l+1]=yr;
	      a_data[dim3*l+2]=zr;
	      
	      a_data[dim3*l+dim2_dim1]=xr;
	      a_data[dim3*l+dim2_dim1+1]=yr;
	      a_data[dim3*l+dim2_dim1+2]=zr;
	      
	      a_data[dim3*l+dim2_dim1*2]=xr;
	      a_data[dim3*l+dim2_dim1*2+1]=yr;
	      a_data[dim3*l+dim2_dim1*2+2]=zr;
	      
	      //assign b vector
	      xr=xyz_1[atom[i]][0]-xyz_1[atom[0]][0];
	      yr=xyz_1[atom[i]][1]-xyz_1[atom[0]][1];
	      zr=xyz_1[atom[i]][2]-xyz_1[atom[0]][2];
	      
	      xr=xr-tild1*round(zr/L);
	      xr=xr-L*round(xr/L);
	      yr=yr-L*round(yr/L);
	      zr=zr-L*round(zr/L);
	      
	      b_data[l*dim]=xr;
	      b_data[l*dim+1]=yr;
	      b_data[l*dim+2]=zr;
	      
	  }
	 

	  solve_linear_eq(a_data,b_data,x_data,dim2);
	  
	  
	  //======== Transform F to delu/delx=======
 	  for(i=0;i<dim;i++)
	    for(j=0;j<dim;j++)
	      x_data[i*dim+j]-=delta(i,j);


	  //======== Make symmetric epsilon[i][j]=======
	  for(i=0;i<dim;i++) 
	    for(j=0;j<dim;j++)
	      epsilon[i][j]=(x_data[i*dim+j]+x_data[j*dim+i])/2.0;
	  
	  
	  //==========   e_dev=epsilon-Tr(epsilon)I/d ====
	  Trepsilon=0.0;
	  for(i=0;i<dim;i++)Trepsilon+=epsilon[i][i];
	  for(i=0;i<dim;i++)
	    for(j=0;j<dim;j++)
	      epsilon_dev[i][j]=epsilon[i][j]-(1.0/dim)*Trepsilon*delta(i,j);    
	  
	  //=========  e_d=sqrt[ Tr(e_dev^2)/2 ]   =======
	  
	  for(i=0;i<dim;i++)
	    for(j=0;j<dim;j++){
	      epsilon_sq[i][j]=0.0;	  
	      for(k=0;k<dim;k++)	      
		epsilon_sq[i][j]+=epsilon_dev[i][k]*epsilon_dev[k][j];
	    }
	  Trepsilon=0.0;
	  for(i=0;i<dim;i++)Trepsilon+=epsilon_sq[i][i];
	  
	  epsilon_d=sqrt(Trepsilon/2.0);
	  
	  
	  
	  //printf("\ned=%e\n",epsilon_d);
	  if(epsilon_d>tolerance && epsilon_d<100){
	    findex=epsilon_d/binwidth;
	    index=log(findex*1.0)/log(binof);
	    bincount[index]++;
	    for(i=0;i<=index;i++)CDFbincount[i]++;
	    tot_count++;
	  }
	  //==========================================================
	  
	  for(iedcut=0;iedcut<no_of_edcut;iedcut++){
	    epsilon_cutoff=epscutoffarray[iedcut];

	    for(i=0;i<dim1;i++){
	      if(epsilon_d>epsilon_cutoff)
		ActiveparticleARRAY[atom[i]][iedcut]=1;
	    } 

	  }
	 

 
	}//loop over tlist: tessilation atom list file
      fclose(fp2);



     
      
      //============================================================================
      
      for(iedcut=0;iedcut<no_of_edcut;iedcut++){
	epsilon_cutoff=epscutoffarray[iedcut];

	printf("\nepsilon_cutoff=%f\n",epsilon_cutoff);
	
	for(i=0;i<N;i++)Activeparticle[i]=ActiveparticleARRAY[i][iedcut];

      
	iplist=0;ncl=0;
      for(i=0;i<N;i++) 
	if(Activeparticle[i]>0)
	PLIST[iplist++]=i;
      
      //printf("ipart=%d\n",iplist);
      // Let's count cluster for the frame!!!
      
      if(iplist==0)continue;
      
      ntot=0;
      for(i=0;i<N;i++)ID[i]=-1;
      for(i=0;i<=N;i++)ns[i]=0;
      for(i=0;i<iplist;i++){
	iroot=PLIST[i];
	if(ID[iroot]==-1){
	  ncl++;
	  ID[iroot]=ncl;
	  n=1;
	  LIST[n]=iroot;
	  k1=n;
	  k2=k1;
	  while(k1<=k2){
	    for(k=k1;k<=k2;k++){
	      iroot=LIST[k];
	      for(j=0;j<iplist;j++){
		inew=PLIST[j];
		if(ID[inew]==-1){
		  
		  xr=xyz_0[iroot][0]-xyz_0[inew][0];
		  yr=xyz_0[iroot][1]-xyz_0[inew][1];
		  zr=xyz_0[iroot][2]-xyz_0[inew][2];
		  
		  xr=xr-tild0*round(zr/L);
		  xr=xr-L*round(xr/L);
		  yr=yr-L*round(yr/L);
		  zr=zr-L*round(zr/L);
		  dr=sqrt(xr*xr+yr*yr+zr*zr);
		  
		  if(dr<rcutoff){
		    ID[inew]=ncl;
		    n++;
		    LIST[n]=inew;
		  }
		}
	      }//loop j
	    }//loop k
	    k1=k2+1;
	    k2=n;
	  }//while
	  ns[n]++;
	  printf("ncl=%d n=%d\n",ncl,n);
	  ntot+=n;
	}//if
      }//loop i

      
      for(i=0;i<=N;i++)
	if(ns[i]>0)tot_ns[i][iedcut]+=ns[i];
      




      }//iedcut for loop


      //============================================================================

 here:

    printf("tot_count=%e\n",tot_count);

      double b1,b2,bw;
      sprintf(FILE,"Distrubution_ed.dat");
      fp=fopen(FILE,"w");
      for(i=0;i<binsize;i++){
      	b1=pow(binof,i);
      	b2=pow(binof,i+1);
      	bw=sqrt(b1*b2);
      	if(bincount[i]>0)
      	  fprintf(fp,"%e %e\n",bw*binwidth,bincount[i]/((b2-b1)*binwidth*tot_count));
      }
      fclose(fp);
      
    
      for(iedcut=0;iedcut<no_of_edcut;iedcut++){
	epsilon_cutoff=epscutoffarray[iedcut];
	
	
	sprintf(FILE,"Dist_ns_edcut%3.2f.dat",epsilon_cutoff);
	fp=fopen(FILE,"w");
	for(i=0;i<=N;i++)if(tot_ns[i][iedcut]>0)fprintf(fp,"%d %d\n",i,tot_ns[i][iedcut]);
      fclose(fp);
      
      
      for(i=0;i<binsize;i++)bincount[i]=0.0;
      tot_count=0.0;
      for(i=1;i<=N;i++){
	index=log(i*1.0)/log(2.0);
	bincount[index]+=tot_ns[i][iedcut];
	tot_count+=tot_ns[i][iedcut];
      }
      
      sprintf(FILE,"bin_Dist_ns_edcut%3.2f.dat",epsilon_cutoff);
      fp=fopen(FILE,"w");
      for(i=0;i<binsize;i++){
	b1=pow(2.0,i);
	b2=pow(2.0,i+1);
	bw=sqrt(b1*b2);
	if(bincount[i]>0)
	  fprintf(fp,"%e %e\n",b1,bincount[i]/((b2-b1)*tot_count));
      }
      fclose(fp);
      
      }// iedcut

      
    free(xyz_0); 
    free(xyz_1); 
    
    sprintf(command, "rm tlist");
    system(command);
    
    return 0;
} //end main


int delta(int i, int j){
  int value;
  if(i==j) value=1;
  else value=0;
  return value;
}

	  
#include <gsl/gsl_linalg.h>
   
int solve_linear_eq (double *a_data, double *b_data, double *x_data, long N){
  
  gsl_matrix_view m = gsl_matrix_view_array (a_data, N, N);
  gsl_vector_view b = gsl_vector_view_array (b_data, N);
  gsl_vector *x = gsl_vector_alloc (N);
  int s;
  gsl_permutation * p = gsl_permutation_alloc (N);
  gsl_linalg_LU_decomp (&m.matrix, p, &s);
  gsl_linalg_LU_solve (&m.matrix, p, &b.vector, x);
  int i; 
  for (i = 0; i < N; ++i)x_data[i]=gsl_vector_get(x,i);
  gsl_vector_free(x);
  gsl_permutation_free(p);
  return 0; 
}
