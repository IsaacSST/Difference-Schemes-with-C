#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>

void read_input(double *L, int *N, double *t_f,double *t_D, double *K, double *A, double *B);

int main(void){
  //****Parameters*****
  int N;
  double L;
  double t_f;
  double t_D;
  double K;
  double A;
  double B;
  read_input(&L, &N, &t_f, &t_D, &K, &A, &B); //reading values from file into variables

  double dx = L/N; //Grid spacing
  double U0 = K + A; //exact solution guesses
  double V0 = B;     //
  double dt0 = 0.5/sqrt(U0*U0 + V0*V0); // stable initial timestep
  
  //****Grid Storage****
  double *U, *U_next; //U at current and next timestep
  double *V, *V_next; //V at current and next timestep
  
  //allocating memory according to size of N
  U       = malloc(N*sizeof(double));
  U_next  = malloc(N*sizeof(double));
  V       = malloc(N*sizeof(double));
  V_next  = malloc(N*sizeof(double));
  
  if (U==NULL||U_next==NULL||V==NULL||V_next==NULL) {
    printf("Memory allocation failed\n");
    return 1;
  }

  int j;
  double x;

  //initialisation
  for(j=0;j<N;j++) {
    x = j*dx;
    U[j]  = K + A*cos(2*acos(-1.0)*x/L); //acos(-1.0) = pi
    V[j]  = B*sin(2*acos(-1.0)*x/L);
  }
  
  FILE *out;
  out = fopen("output.txt","w");
  
  
  // Output at start of simulation (copied from assignment 3)
  double ctime = 0.0;
  for (j=0; j<N; j++ ) {
    x = j*dx;
    printf("%g %g %g %g\n",ctime,x,U[j],V[j]);
    fprintf(out, "%g %g %g %g\n",ctime,x,U[j],V[j]);
  }
  double next_output_time = t_D;


  
  //loop over timesteps 
  while (ctime<t_f){
    //adjusting exact solution guess for each iteration
    U0 = U[1];
    V0 = V[1];
    //adjusting timestep with each iteration for stability
    dt0 = 0.5/sqrt(U0*U0 + V0*V0);
 
    
    int output = 0;
    // If we would go past the next output step, reduce the timestep.
    if (ctime+dt0>next_output_time) {
      dt0 = next_output_time - ctime;
      output = 1;
    }
    

    
    for (j=0; j<N; j++) {
      x = j*dx;
      
      double U2x; //second order x derivatives
      double V2x;
      double common = 1 - 2*(U0*U0) - 2*(V0*V0); // common term in the equations, expressed as variable to aid clarity.

      // conditions for spacial derivative to account for periodic boundaries
      if (j == N-1){
        U2x = (U[0] + U[j-1] - 2*U[j])/(dx*dx);
        V2x = (V[0] + V[j-1] - 2*V[j])/(dx*dx);
      }else if (j == 0){
        U2x = (U[j+1] + U[N-1] - 2*U[j])/(dx*dx);
        V2x = (V[j+1] + V[N-1] - 2*V[j])/(dx*dx);
      }else{
        U2x = (U[j+1] + U[j-1] - 2*U[j])/(dx*dx);
        V2x = (V[j+1] + V[j-1] - 2*V[j])/(dx*dx);
      }
      //scheme
      U_next[j] = U[j] + dt0*(U2x + U[j]*(2*common - 8*U0*U0) - 8*U0*V0*V[j]);
      V_next[j] = V[j] + dt0*(V2x + V[j]*(2*common - 8*V0*V0) - 8*U0*V0*U[j]);
     }
    

    
	   
    // Efficiently copy next values at timestep to U,V arrays (copied from assignment 3)
    double *tmp;
    tmp = V_next;
    V_next = V; 
    V = tmp;
    tmp = U_next;
    U_next = U; 
    U = tmp;

    // Increment time modified from assignment 3  
    ctime += dt0;
    if (output) {
      for (j=0; j<N; j++ ) {
	x = j*dx;
        printf("%g %g %g %g\n",ctime,x,U[j],V[j]);
        fprintf(out, "%g %g %g %g\n",ctime,x,U[j],V[j]);
      }
      next_output_time += t_D;
    }
  }
  //freeing pointers
  free(U);
  free(U_next);
  free(V);
  free(V_next);
  fclose(out);
  return 0;
}

//Read input function
void read_input(double *L, int *N, double *t_f,double *t_D, double *K, double *A, double *B) {
   FILE *infile;
   int placeholder = 0;
   int currentlength;
   int i;
   int j;
   char readstring[100];
   char storevalues[1000];
   // opening file
   if(!(infile=fopen("input.txt","r"))) {
       printf("Error opening file\n");
       exit(1);
   }
  
   // for loop takes input values from seperate lines in input file and puts them in a single string together, in order of appearance in file and seperated by spaces
   for (i = 0;i<7;i++){
     if( fgets(readstring,1000, infile)!=NULL ) { //ith line of file is read to a string (fgets automatically iterates through lines every time it is called until file is closed)
       currentlength = strlen(readstring); //length of input line in file (this includes the new line character)
       
       for (j = 0;j< currentlength-1;j++){
         storevalues[j + placeholder] = readstring[j]; //current input copied into holding string in appropriate position, excluding new line character at the end
       }
       
       if (i == 6){
         storevalues[placeholder + currentlength -1] = '\0'; //adding string terminator when all inputs are in the holding string
       }else{
         storevalues[placeholder + currentlength -1] = ' '; //adding space between values
       }
       placeholder = placeholder + currentlength; // place shifted to index after the space after previous input
       
     }
   }
   // reading the inputs from the holding string into the variables
   if(7!=sscanf(storevalues,"%lf %d %lf %lf %lf %lf %lf",L,N,t_f,t_D,K,A,B)) {
       printf("Error reading parameters from file\n");
       exit(1);
   }
   fclose(infile);
}
