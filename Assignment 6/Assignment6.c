#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <lapacke.h>
//------------------Matrix routines copied from lecture 10 example code---------------------------
struct band_mat{
   long ncol;        /* Number of columns in band matrix            */
   long nbrows;      /* Number of rows (bands in original matrix)   */
   long nbands_up;   /* Number of bands above diagonal           */
   long nbands_low;  /* Number of bands below diagonal           */
   double *array;    /* Storage for the matrix in banded format  */
   /* Internal temporary storage for solving inverse problem */
   long nbrows_inv;  /* Number of rows of inverse matrix   */
   double *array_inv;/* Store the inverse if this is generated */
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

/* Finalise function: should free memory as required */
void finalise_band_mat(band_mat *bmat) {
   free(bmat->array);
   free(bmat->array_inv);
   free(bmat->ipiv);
}


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

void setv(band_mat *bmat, long row, long column, double val) {
   *getp(bmat,row,column) = val;
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




/* Return the 1D element index corresponding to a particular grid point.
   We can rewrite this function without changing the rest of the code if
   we want to change the grid numbering scheme!
   Output: long integer with the index of the point
   Input:
   long j:  The X grid point index
   long p:  The Y grid point index
   long P:  The number of Y points.
*/
long indx( long j, long p, long J) {
   return p*J + j;
}



//----------------------------------------------------------------------------------------------------
//----------------------------------------------------------------------------------------------------
void read_input(long *Nx, long *Ny, double *Lx,double *Ly, double *tf, double *lam, double *sig, double *kap, double *A0, double *A1, double *A2, double *tD);
long reindex(long p, long N);
long indexer(long j, long p, int switched, long Nx, long Ny);
int main(void){
   //output file
   FILE *out;
   out = fopen("output.txt","w");
   //parameters
   long Nx;
   long Ny;
   double Lx;
   double Ly;
   double tf;
   double lam;
   double sig;
   double kap;
   double A0;
   double A1;
   double A2;
   double tD;
   //reading values from  input file into variables
   read_input(&Nx, &Ny, &Lx, &Ly, &tf, &lam, &sig, &kap, &A0, &A1, &A2, &tD);
   //grid spacing
   double dx = Lx/Nx; //because x is periodic so don't need to include value x=Lx in grid as it is the same as x=0
   double dy = Ly/(Ny-1); //so Ly is included in the grid
   
   //timestep
   double dtin = 0.8342714479759265;////almost exactly the stability threshold for example input
   double dt = dtin;
    
   //initial condition coefficients A0 and A1 greatly affect stability:adjusting timestep accordingly
   double As = fabs(A0) + fabs(A1);
   double dtAs = (2.0)*(dtin/(As*As));
   
   double dtlam = dtin;
   double dtsig = dtin;
   double dtkap = dtin;
   //stable timesteps for each input parameter obtained by increasing each input in the example input one by one
   //-nb these are the only parameters that significantly affect stability
   if (fabs(lam) > 5.0){
      dtlam = (3.0)*dtin/(fabs(lam));
   }
   
   if (fabs(sig) > 20.0){
      dtsig = (fmin((2.0),fabs(sig) - 20.0))*dtin/(sqrt(fabs(sig)));
   }
   if (fabs(kap) > 8.0){
      dtkap = (5.0)*dtin/fabs(kap);
   }
   
   //taking the timestep to be the minimum of all the stable timesteps for each parameter, so it is guaranteed to be stable
   dt = fmin(fmin(dtAs,dtlam),fmin(dtsig,dtkap));
  
   //band matrices
   band_mat bmatx;
   band_mat bmaty;
   long ncols = Nx*Ny;
   //two matrices for time splitting scheme, much quicker to solve two matrices with few bands than one matrix with 2*Nx bands when Nx is large
   //two bands either side for x direction, as need to reindex x due to periodic boundary conditions
   init_band_mat(&bmatx, 2, 2, ncols);
   //y direction matrix is tridiagonal, as neumann boundary conditions don't require reindexing
   init_band_mat(&bmaty, 1, 1, ncols);

   //vectors for holding source terms
   double *bU = malloc(sizeof(double)*ncols);
   double *bV = malloc(sizeof(double)*ncols);
   //vectors for solving for U and V
   double *U = malloc(sizeof(double)*ncols);
   double *V = malloc(sizeof(double)*ncols);
   
   
   //iterands
   long i;
   long j;
   long k;
   long p;
   //storing value for pi
   double pi = acos(-1);
   
   //initial conditions and printing for time 0
   for (p = 0;p<Ny;p++){
      for (j = 0;j<Nx;j++){
         U[indexer(j,p,1,Nx,Ny)] = (A0 + A1*cos(pi*p*dy/Ly))*exp((-1)*(A2*A2*j*j*dx*dx/(Lx*Lx)));
         V[indexer(j,p,1,Nx,Ny)] = 0.0;
         fprintf(out,"%.10lf %.10lf %.10lf %.10lf \n",j*dx,p*dy,U[indexer(j,p,1,Nx,Ny)],V[indexer(j,p,1,Nx,Ny)]);
      }
   }

   //time variables
   double next_output_time = tD;
   double ctime = 0.0;
   double dt0;
   double dt1;
   //array to put coefficents in so code is easier to read
   double co[3];
   //setting common constant coefficients for putting into matrix
   co[0] = 1/(dx*dx);
   co[1] = 1/(dy*dy);
   
   //grid variables
   long leftx;
   long lefty;
   long pointprev;
   long point;
   long rightx;
   long righty;

   //output boolean
   int output;
   
   while (ctime<tf ){
      //only want to output at diagnostic timesteps so default output = 0
      output = 0;
	  
      dt0 = dt;
      //If we would go past the next output step, reduce the timestep.
      if (ctime+dt0>next_output_time) {
         dt0 = next_output_time - ctime;
         output = 1;
      }
      //output at final time
      if (ctime + dt0 > tf){
         dt0 = tf - ctime;
         output = 1;
      }
      
      //time split
      dt1 = dt0/2;
      //setting diagonal coefficient for x-direction matrix
      co[2] = (-1)*((1/dt1) + (2/(dx*dx)));
      
      //solving implicitly in x direction to get t+1/2 step
      for (j=0;j<Nx;j++){
         for (p=0;p<Ny;p++){
	     //indexing, taking into account boundary conditions
             point = indexer(j,p,0,Nx,Ny);
             pointprev = indexer(j,p,1,Nx,Ny);
             
             if (j == 0){
                leftx = indexer(Nx - 1,p,0,Nx,Ny);
             }else{
                leftx = indexer(j - 1,p,0,Nx,Ny);
             }

             if (j == Nx - 1){
                rightx = indexer(0,p,0,Nx,Ny);
             }else{
                rightx = indexer(j + 1,p,0,Nx,Ny);
             }
             
             
             //source terms
             bU[point] = sig*(V[pointprev]) + kap + (U[pointprev])*(U[pointprev])*(U[pointprev]) + (-1)*(lam + (1/dt1))*(U[pointprev]);
             bV[point] = (1 - (1/dt1))*(V[pointprev]) - U[pointprev];
             //matrix
             setv(&bmatx,point,point,co[2]);
             setv(&bmatx,point,leftx,co[0]);
             setv(&bmatx,point,rightx,co[0]);
           
         }
     }
     
     //solving for split timestep
     solve_Ax_eq_b(&bmatx, U, bU);
     solve_Ax_eq_b(&bmatx, V, bV);

     co[2] = (-1)*((1/dt1) + (2/(dy*dy)));
     //solving in y direction implicitly based on t+1/2 values;
     for (j=0;j<Nx;j++){
         for (p=0;p<Ny;p++){
	     //indexing, taking into account boundary conditions
             point = indexer(j,p,1,Nx,Ny);
             pointprev = indexer(j,p,0,Nx,Ny);
             if (p != 0){
                lefty = indexer(j,p - 1,1,Nx,Ny);
             }
             if (p != Ny - 1){
                righty = indexer(j,p + 1,1,Nx,Ny);
             }
             
             //source terms
             bU[point] = (-1)*(1/dt1)*(U[pointprev]);
             bV[point] = (-1)*(1/dt1)*(V[pointprev]);
             //matrix
             setv(&bmaty,point,point,co[2]);
             
             if (p == 0){
                setv(&bmaty,point,righty,2*co[1]);
             }else if (p == Ny -1){
                setv(&bmaty,point,lefty,2*co[1]);
             }else{
	        setv(&bmaty,point,lefty,co[1]);
                setv(&bmaty,point,righty,co[1]);
             }
         }
     }
     
     //solving for next timestep
     solve_Ax_eq_b(&bmaty, U, bU);
     solve_Ax_eq_b(&bmaty, V, bV);
	 
     //printing to file if at output timestep
     if (output == 1){
        for(i=0;i<Ny; i++) {
           for(k=0;k<Nx;k++){
              fprintf(out,"%.10lf %.10lf %.10lf %.10lf \n",k*dx,i*dy,U[indexer(k,i,1,Nx,Ny)],V[indexer(k,i,1,Nx,Ny)]);
           }
        }
        //incrementing output time
        next_output_time += tD;
     }

     //incrementing simulation time 
     ctime += dt0;
      
   }
   
   
   //freeing heap variables
   free(bU);
   free(bV);
   free(U);
   free(V);
   finalise_band_mat(&bmatx);
   finalise_band_mat(&bmaty);
   fclose(out);
   return 0;
}

//indexer that depends on which direction is being solved for (or has been solved for if trying to access the solution vector). If switched = 1 then y direction indexing scheme is used, if switched = 0 then its the x direction.
long indexer(long j, long p, int switched, long Nx, long Ny){
   if (switched){
      //indexing scheme so that matrix has Nx by Nx blocks of size Ny by Ny, i.e for a given x, the y-1,y,y+1 points are all next to each other, which allows for tridiagonal banded matrix when solving in y direction
      return indx(p,j,Ny);
   }else{
      //indexing scheme so that matrix has Ny by Ny blocks of size Nx by Nx for the same reason as with y direction. However x is also reindexed to deal with periodic boundary, leading to a 5-diagonal banded matrix for x direction. 
      return indx(reindex(j,Nx),p,Nx);
   }
}
//standard 'fold back' reindexer for x direction matrix to deal with periodic boundary. roughly the same as last assignment's reindexer but modified to deal with odd Nx
long reindex(long p,long N){
   //writing as N - N/2 deals with the case with odd N, since division of itegers in C rounds down to nearest whole integer i.e N/2 =/= (N - N/2)
   if (p< (N - (N/2))){
      return 2*p;
   }else{
      
      return 2*N - 1 - 2*p;
   }
}

//Read input function, mostly copied from my previous assignment, just with more variables
void read_input(long *Nx, long *Ny, double *Lx, double *Ly, double *tf, double *lam, double *sig, double *kap, double *A0, double *A1, double *A2, double *tD) {
  FILE *infile;
   // opening file
   if(!(infile=fopen("input.txt","r"))) {
       printf("Error opening file\n");
       exit(1);
   }
   //reading parameters
   if(1!=fscanf(infile,"%ld\n",Nx)) {
       printf("Error reading parameters from file\n");
       exit(1);
   }
   if(1!=fscanf(infile,"%ld\n",Ny)) {
       printf("Error reading parameters from file\n");
       exit(1);
   }
   if(1!=fscanf(infile,"%lf\n",Lx)) {
       printf("Error reading parameters from file\n");
       exit(1);
   }
   if(1!=fscanf(infile,"%lf\n",Ly)) {
       printf("Error reading parameters from file\n");
       exit(1);
   }
   if(1!=fscanf(infile,"%lf\n",tf)) {
       printf("Error reading parameters from file\n");
       exit(1);
   }
   if(1!=fscanf(infile,"%lf\n",lam)) {
       printf("Error reading parameters from file\n");
       exit(1);
   }
   if(1!=fscanf(infile,"%lf\n",sig)) {
       printf("Error reading parameters from file\n");
       exit(1);
   }
   if(1!=fscanf(infile,"%lf\n",kap)) {
       printf("Error reading parameters from file\n");
       exit(1);
   }
   if(1!=fscanf(infile,"%lf\n",A0)) {
       printf("Error reading parameters from file\n");
       exit(1);
   }
   if(1!=fscanf(infile,"%lf\n",A1)) {
       printf("Error reading parameters from file\n");
       exit(1);
   }
   if(1!=fscanf(infile,"%lf\n",A2)) {
       printf("Error reading parameters from file\n");
       exit(1);
   }
   if(1!=fscanf(infile,"%lf\n",tD)) {
       printf("Error reading parameters from file\n");
       exit(1);
   }
   
   fclose(infile);
}
