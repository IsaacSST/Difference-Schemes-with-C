#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <lapacke.h>



//***********************MATRIX FUNCTIONS-taken from lecture notes ************************//
struct band_mat{
  long ncol;        /* Number of columns in band matrix            */
  long nbrows;      /* Number of rows (bands in original matrix)   */
  long nbands_up;   /* Number of bands above diagonal           */
  long nbands_low;  /* Number of bands below diagonal           */
  double *array;    /* Storage for the matrix in banded format  */
  /* Internal temporary storage for solving inverse problem */
  long nbrows_inv;  /* Number of rows of inverse matrix   */
  double *array_inv;/* Store the matrix decomposition if this is generated:                */
                    /* this is used to calculate the action of the inverse matrix.         */
                    /* (what is stored is not the inverse matrix but an equivalent object) */
  int *ipiv;        /* Additional inverse information         */
};
typedef struct band_mat band_mat;

/* Initialise a band matrix of a certain size, allocate memory,
   and set the parameters.  */ 
int init_band_mat(band_mat *bmat, long nbands_lower, long nbands_upper, long n_columns) {
  bmat->nbrows = nbands_lower + nbands_upper + 1;
  bmat->ncol   = n_columns;
  bmat->nbands_up = nbands_upper;
  bmat->nbands_low= nbands_lower;
  bmat->array      = (double *) malloc(sizeof(double)*bmat->nbrows*bmat->ncol);
  bmat->nbrows_inv = bmat->nbands_up*2 + bmat->nbands_low + 1;
  bmat->array_inv  = (double *) malloc(sizeof(double)*(bmat->nbrows+bmat->nbands_low)*bmat->ncol);
  bmat->ipiv       = (int *) malloc(sizeof(int)*bmat->ncol);
  if (bmat->array==NULL||bmat->array_inv==NULL) {
    return 0;
  }  
  /* Initialise array to zero */
  long i;
  for (i=0;i<bmat->nbrows*bmat->ncol;i++) {
    bmat->array[i] = 0.0;
  }
  return 1;
};


/* Get a pointer to a location in the band matrix, using
   the row and column indexes of the full matrix.           */
double *getp(band_mat *bmat, long row, long column) {
  int bandno = bmat->nbands_up + row - column;
  if(row<0 || column<0 || row>=bmat->ncol || column>=bmat->ncol ) {
    printf("Indexes out of bounds in getp: %ld %ld %ld \n",row,column,bmat->ncol);
    exit(1);
  }
  return &bmat->array[bmat->nbrows*column + bandno];
}

/* Retrun the value of a location in the band matrix, using
   the row and column indexes of the full matrix.           */
double getv(band_mat *bmat, long row, long column) {
  return *getp(bmat,row,column);
}

double setv(band_mat *bmat, long row, long column, double val) {
  *getp(bmat,row,column) = val;
  return val;
}

/* Solve the equation Ax = b for a matrix a stored in band format
   and x and b real arrays                                          */
int solve_Ax_eq_b(band_mat *bmat, double *x, double *b) {
  /* Copy bmat array into the temporary store */
  int i,bandno;
  for(i=0;i<bmat->ncol;i++) { 
    for (bandno=0;bandno<bmat->nbrows;bandno++) {
      bmat->array_inv[bmat->nbrows_inv*i+(bandno+bmat->nbands_low)] = bmat->array[bmat->nbrows*i+bandno];
    }
    x[i] = b[i];
  }

  long nrhs = 1;
  long ldab = bmat->nbands_low*2 + bmat->nbands_up + 1;
  int info = LAPACKE_dgbsv( LAPACK_COL_MAJOR, bmat->ncol, bmat->nbands_low, bmat->nbands_up, nrhs, bmat->array_inv, ldab, bmat->ipiv, x, bmat->ncol);
  return info;
}

int printmat(band_mat *bmat) {
  long i,j;
  for(i=0; i<bmat->ncol;i++) {
    for(j=0; j<bmat->nbrows; j++) {
       printf("%ld %ld %g \n",i,j,bmat->array[bmat->nbrows*i + j]);
    }
  }
  return 0;
}
//**********************************************************************************************//
//**********************************************************************************************//
void read_input(double *L, long *N, double *K ,double *gam, double *omA, double *omB, long *Iom);
void read_coefficients(double **mu, double **E, long N);
//reindexing function
long reindex(long j, long ncols);

int main(void){
   //output file
   FILE *out;
   out = fopen("output.txt","w");
   //parameters
   double L;
   long N;
   double K;
   double gam;
   double omA;
   double omB;
   long Iom;
   //reading values from  input file into variables
   read_input(&L, &N, &K, &gam, &omA, &omB, &Iom); 
   //coefficients
   double *mu = malloc(sizeof(double)*N);
   double *E = malloc(sizeof(double)*N);
   
   // reading values from coefficients file into arrays
   read_coefficients(&mu,&E,N);
   //grid spacing
   double dx = L/N;

   //band matrix
   band_mat bmat;
   long ncols = 2*N;
   long nbands_low = 4;
   long nbands_up = 4;
   init_band_mat(&bmat, nbands_low, nbands_up, ncols);
   double *x = malloc(sizeof(double)*ncols);
   double *b = malloc(sizeof(double)*ncols);
   
   //change in omega between each loop
   double dom = (omB - omA)/Iom;
   //initialising omega
   double om = omA;

   
   //arbitrary iteration variables
   long currentU;
   long currentV;
   long leftU;
   long leftV;
   long rightU;
   long rightV;
   long i;
   long j;
   
   //array to put coefficents in so code is easier to read
   double co[4];

   //frequency scan loop
   for (i = 0; i<Iom;i++){
      //setting coupling coefficent (this is a coefficient that is constant throughout the grid points)
      co[2] = om*gam;
      if (N == 1){
         //solution for N=1, from pen and paper method
         co[0] = 2*E[0]*cos(K*L)/(dx*dx) + om*om*mu[0] - 2*E[0]/(dx*dx);
         co[1] = co[0]*co[0] + co[2]*co[2];
         x[0] = (-1)*co[0]/co[1];
         x[1] = co[2]/co[1];
         
      }else{
       for (j=0;j<N;j++){
	  
	  //setting coefficients for putting into matrix
	  co[0] = E[j]/(dx*dx);
	  if (j < N-1){
	     co[1] = om*om*mu[j] - (E[j] + E[j+1])/(dx*dx);
	     co[3] = E[j+1]/(dx*dx);
	  }else{
             //dealing with right boundary coefficients using periodicity of E
	     co[1] = om*om*mu[j] - (E[j] + E[0])/(dx*dx);
	     co[3] = E[0]/(dx*dx);
	  }
          //matrix entries at left boundary
	  if (j == 0){
	     //reindexing
	     leftU = reindex(2*N-2,ncols);
	     leftV = reindex(2*N-1,ncols);
	     currentU = reindex(2*j,ncols);
	     currentV = reindex(2*j+1,ncols);
	     rightU  = reindex(2*j+2,ncols);
	     rightV  = reindex(2*j+3,ncols);
             
	     //----------real part-------------
	     b[reindex(2*j,ncols)] = (-1)*cos(0);
	     setv(&bmat,currentU,leftU,co[0]*cos(K*L));
             setv(&bmat,currentU,leftV,co[0]*sin(K*L));
	     setv(&bmat,currentU,currentU,co[1]);
	     setv(&bmat,currentU,currentV,(-1)*co[2]);
	     setv(&bmat,currentU,rightU,co[3]);
	  
	  
	     //-------------imaginary part------------
	     b[reindex(2*j+1,ncols)] = (-1)*sin(0);
             setv(&bmat,currentV,leftU,co[0]*(-1)*sin(K*L));
	     setv(&bmat,currentV,leftV,co[0]*cos(K*L));
	     setv(&bmat,currentV,currentV,co[1]);
	     setv(&bmat,currentV,currentU,co[2]);
	     setv(&bmat,currentV,rightV,co[3]);
          //matrix entries at right boundary
	  }else if (j == N-1){
	     //reindexing
	     leftU = reindex(2*j-2,ncols);
	     leftV = reindex(2*j-1,ncols);
	     currentU = reindex(2*j,ncols);
	     currentV = reindex(2*j+1,ncols);
	     rightU  = reindex(0,ncols);
	     rightV  = reindex(1,ncols);

	     //----------real part-------------
	     b[reindex(2*j,ncols)] = (-1)*cos(K*j*dx);
	     setv(&bmat,currentU,leftU,co[0]);
	     setv(&bmat,currentU,currentU,co[1]);
	     setv(&bmat,currentU,currentV,(-1)*co[2]);
	     setv(&bmat,currentU,rightU,co[3]*cos(K*L));
	     setv(&bmat,currentU,rightV,co[3]*(-1)*sin(K*L));
	  
	  
	     //-------------imaginary part------------
	     b[reindex(2*j+1,ncols)] = (-1)*sin(K*j*dx);
	     setv(&bmat,currentV,leftV,co[0]);
	     setv(&bmat,currentV,currentV,co[1]);
	     setv(&bmat,currentV,currentU,co[2]);
             setv(&bmat,currentV,rightU,co[3]*sin(K*L));
	     setv(&bmat,currentV,rightV,co[3]*cos(K*L));
          //everything in between the boundaries
	  }else{
	  
	    
	     //reindexing
	  
	     leftU = reindex(2*j-2,ncols);
	     leftV = reindex(2*j-1,ncols);
	     currentU = reindex(2*j,ncols);
	     currentV = reindex(2*j+1,ncols);
	     rightU  = reindex(2*j+2,ncols);
	     rightV  = reindex(2*j+3,ncols);
             
	     //----------real part-------------
	     b[reindex(2*j,ncols)] = (-1)*cos(K*j*dx);
	     setv(&bmat,currentU,leftU,co[0]);
	     setv(&bmat,currentU,currentU,co[1]);
	     setv(&bmat,currentU,currentV,(-1)*co[2]);
	     setv(&bmat,currentU,rightU,co[3]);
	  
	  
	     //-------------imaginary part------------
	     b[reindex(2*j+1,ncols)] = (-1)*sin(K*j*dx);
	     setv(&bmat,currentV,leftV,co[0]);
	     setv(&bmat,currentV,currentV,co[1]);
	     setv(&bmat,currentV,currentU,co[2]);
	     setv(&bmat,currentV,rightV,co[3]);
	  }
      }
     }
     //solving for cases apart from if N = 1, which is a special case handled already
     if (N != 1){
         solve_Ax_eq_b(&bmat, x, b);

     }
     //outputting to file
     for(int k=0; k<N; k++) {
		
        fprintf(out,"%.10lf %.10lf %.10lf %.10lf \n",om,k*dx,x[reindex(2*k,ncols)],x[reindex(2*k + 1,ncols)]);
     }
     om +=dom;
   }
   //freeing heap memory
   free(bmat.array);
   free(bmat.array_inv);
   free(bmat.ipiv);
   free(mu);
   free(E);
   free(x);
   free(b);
   fclose(out);

   return 0;
}

//reindexing function
long reindex(long j,long ncols){
   if (j<(ncols/2)){
      return 2*j;
   }else{
      //intuitive way of writing reindexing(for me)
      return ncols - 1 - 2*(j-ncols/2);
   }
}

//reading in coefficient functions
void read_coefficients(double **mu, double **E, long N) {
   FILE *infile;
   int i;
   double readmu;
   double readE;
   // opening file
   if(!(infile=fopen("coefficients.txt","r"))) {
       printf("Error opening file\n");
       exit(1);
   }
   //reading from file
   for (i = 0;i<N;i++){
       if (2!= fscanf(infile,"%lf %lf\n",&readmu,&readE)){
          printf("Error reading coefficients from file\n");
       }
       (*mu)[i] = readmu;
       (*E)[i] = readE;
   }
  

   fclose(infile);
   

}

//Read input function, changed so its much simpler than my previous assignment
void read_input(double *L, long *N, double *K,double *gam, double *omA, double *omB, long *Iom) {
   FILE *infile;
   // opening file
   if(!(infile=fopen("input.txt","r"))) {
       printf("Error opening file\n");
       exit(1);
   }
  
   if(1!=fscanf(infile,"%lf\n",L)) {
       printf("Error reading parameters from file\n");
       exit(1);
   }
   if(1!=fscanf(infile,"%ld\n",N)) {
       printf("Error reading parameters from file\n");
       exit(1);
   }
   if(1!=fscanf(infile,"%lf\n",K)) {
       printf("Error reading parameters from file\n");
       exit(1);
   }
   if(1!=fscanf(infile,"%lf\n",gam)) {
       printf("Error reading parameters from file\n");
       exit(1);
   }
   if(1!=fscanf(infile,"%lf\n",omA)) {
       printf("Error reading parameters from file\n");
       exit(1);
   }
   if(1!=fscanf(infile,"%lf\n",omB)) {
       printf("Error reading parameters from file\n");
       exit(1);
   }
   if(1!=fscanf(infile,"%ld\n",Iom)) {
       printf("Error reading parameters from file\n");
       exit(1);
   }
   fclose(infile);
}
