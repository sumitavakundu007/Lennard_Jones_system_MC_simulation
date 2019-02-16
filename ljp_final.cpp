// Monte-Carlo Simulation for simulating the Lennard-Jones System and calculation of Radial-Distribution Function
// written by Sumitava Kundu

#include<iostream>
#include<math.h>
#include<cmath>
#include<cstdlib>
#include<ctime>
#include<fstream>
#include <stdio.h>
#include <iomanip>
#include<sstream>
#include<string>

using namespace std;
#define _USE_MATH_DEFINES

void read_coord(char *file_name, int n, double *x, double *y, double *z);

void input_file(int n, double *epsilon, double *rho, int *binNo, double *kb, double *T, double *delta, double *sigma, double *rc, double *x, double *y, double *z){   // input file where the code reads the parameter from
  int row = 9, col = 2;   // there are 9 rows and 2 columns in the "parameter.dat" file
  double A[row][col];
  string line;
  char str[67];
  char *file_name = new char;   // defining file_name as a pointer
  ifstream infile;
  infile.open("parameter.dat");   // opens the data file named "parameter.dat"
  if (!infile) {   // if code is unable to find the data file
        cout << "Unable to open file";
        exit(1); // terminate with error
  }
  int counter = 1;
  while((counter < row) && (infile.good())){
    for (int i = 0; i < col; i++){
      if (i == 0){
        infile >> str;   // reads the string at column 1 of each row
      }
      else {
        infile >> A[counter][i];   // read and store the parameters of each row 2nd column in array A
      }
    }
    //cout << A[counter][1] << endl;
    counter++;
  }
  //stores the parameter in the variables respectively
  *epsilon = A[1][1];
  *rho = A[2][1];
  *binNo = A[3][1];
  *kb = A[4][1];
  *T = A[5][1];
  *delta = A[6][1];
  *sigma = A[7][1];
  *rc = A[8][1];

  counter = 9;   // the last row
  for (int i = 0; i < col; i++){
    if (i == 0){
      infile >> str;   // read and store the string of last row 1st column
    }
    else {
      infile >> file_name;   // store the last row 2nd column element in a string called file_name
      read_coord(file_name, n, x, y, z);   // calling the function which reads the co-ordinate data file in the xyz format
    }
  }

  infile.close();
}

//======================================================================

void read_coord(char *file_name, int n, double *x, double *y, double *z){
  double C[n+3][4];
  char str[67];
  string line;
  ifstream infile;
  infile.open(file_name);   // open the file
  if (!infile) {
        cout << "Unable to open file" << endl;
        exit(1); // terminate with error
  }
  int counter = 0;
  while((counter < 2) && (infile.good())){   // reading the first two lines of xyz file
      getline (infile, line);
      counter++;
  }
  counter = 2;
  while(counter < (n+2) && (infile.good())){   // reading the rest of lines from xyz file
      for (int j = 0; j < 4; j++){
        if (j == 0){
          infile >> str;   // reads the string of each row 1st column
          //cout << str << "   \t   ";
        }
        else {
          infile >> C[counter][j];   // reads the 2nd, 3rd,  and 4th column of each row
        }
      }

      int j = 3;
      // incorporating the values inside x, y, z
      x[counter - 2] = C[counter][j-2];   // stores the elements of 2nd column of each row in the variable "x"
      y[counter - 2] = C[counter][j-1];   // stores the elements of 3rd column of each row in the variable "y"
      z[counter - 2] = C[counter][j];   // stores the elements of 4th column of each row in the variable "z"
      //cout << x[counter - 2] << "   \t   " << y[counter - 2] << "   \t   " << z[counter - 2] << endl;
      counter++;
  }
  infile.close();
}

//=====================================================================

double periodic_boundary_condition(double pair_distance, double box_length){   // Periodic boundary condition implemented
 if (pair_distance < - box_length/2){
  pair_distance = pair_distance + box_length;
 }
 else if (pair_distance > box_length/2){
  pair_distance = pair_distance - box_length;
 }
 return pair_distance;
}

//====================================================================

double system_total_energy(int n, double box_length, double epsilon, double sigma, double rc, double *x, double *y, double *z){   // function for calculating the total system energy
  double energy_total = 0.0, E_stepwise;
  for (int i = 0; i < (n-1); i++){
    double W_stepwise = 0.0;
    for (int j = i+1; j < n; j++){
        double x_diff = periodic_boundary_condition((x[i] - x[j]), box_length);
        double y_diff = periodic_boundary_condition((y[i] - y[j]), box_length);
        double z_diff = periodic_boundary_condition((z[i] - z[j]), box_length);
        double R2 = (x_diff*x_diff + y_diff*y_diff + z_diff*z_diff);   // calculating the distance between a pair

        if (R2 <= rc*rc)   // condition for checking the pair distance within the cut-off distance
        {
          E_stepwise = 4*epsilon*(pow(sigma, 12)/pow(R2, 6) - pow(sigma, 6)/pow(R2, 3));  // calculation of Lennard-Jones Potential
        }
        else {
          E_stepwise = 0;
        }
        W_stepwise = W_stepwise + E_stepwise;   // calculates the total energy for one selected particle
    }
    energy_total = energy_total + W_stepwise;   // calculates the total energy of the system
  }
  return energy_total;
}
//==========================================================================

void rdf(int n, int binNo, double binSize, double box_length, double *count, double *x, double *y, double *z){
  double dR;
  for (int i = 0; i < (n-1); i++){
    for (int j = i+1; j < n; j++){
      double dx = periodic_boundary_condition((x[i] - x[j]), box_length);
      double dy = periodic_boundary_condition((y[i] - y[j]), box_length);
      double dz = periodic_boundary_condition((z[i] - z[j]), box_length);
      dR = pow((dx*dx + dy*dy + dz*dz), 0.5);   // calculating the magnitude of each pair-wise distances
      double initial = 0.0;
      for (int k = 0; k < binNo; k++){
        if (dR <= box_length/2 && dR > initial && dR <= (initial + binSize)){
          count[k] = count[k] + 1.0;
        }
        initial = initial + binSize;
      }
    }
  }
}

//==========================================================================

int main(){
  //declearing the variables
  int n = 999, mcMAX = 1000, m = 100, binNo;
  // 'acc' is a counter which will increase with particle acceptance. Finally I will divide 'acc' by maximum number of Monte-Carlo steps to get
  // the average acceptance ratio
  // 'delta' is used to get the random moves of x, y, z
  double K, box_length, binSize, acc = 0, epsilon, rho, kb, T, delta, sigma, rc;
  // declearing x, y, z as pointer array
  // Because I need to have x[0], y[0], z[0], ......... x[N], y[N], z[N] at a time to calculate the energy before displacement
  double *x = new double[n];
  double *y = new double[n];
  double *z = new double[n];

  srand(time(0));
  // Reading the input parameters from 'parameter.dat' file
  input_file(n, &epsilon, &rho, &binNo, &kb, &T, &delta, &sigma, &rc, x, y, z);
  // each  parameter is taken as a pointer and bringing all of them including the values of x, y, z from the func 'input_file'
  // This will have x[0], y[0], z[0], ......... x[N], y[N], z[N]
  // Now I am sending all these values of x, y, z  along with the other variables to the energy calculation function to calculate the total system energy
  double *count = new double[binNo];   // declearing 'count' as pointer array

  // Setting the initial values of 'count', which I need during RDF calculation
  for (int i = 0; i < binNo; i++){
     count[i] = 0;
  }

  box_length =  cbrt(n/rho);   // determination of Box length, L
  binSize = (box_length*0.5)/double(binNo);   // binsize calculation

  //Opening the files to write
  ofstream outffile;
  outffile.open("coordinates.dat");
  ofstream outfile;
  outfile.open("energy_vs_mc_time_step_0.2_temp_0.723.dat");
  ofstream outFile;
  outFile.open("acceptance_ratio_0.2_temp_0.723.dat");
  ofstream outFfile;
  outFfile.open("rdf_0.2_temp_0.723.dat");

  // calculating the total system energy before monte carlo move
  double energy_before_disp = system_total_energy(n, box_length, epsilon, sigma, rc, x, y, z);   // calculating the total system energy before monte carlo move

  outfile << "0" << " \t   " << energy_before_disp << endl;   // print the initial total energy in the provided output file

  // Initializing Monte Carlo Loops
  for (int mc = 1; mc < mcMAX; mc++){
    int k = rand() % 1000;   // choosing a random particle which I will take as basis particle.
    // It will change in each MC loop
    double ranf = double(rand())/double(RAND_MAX);
    // Random number of dx, dy, dz implemented to give x, y, z random and different moves
    double dx = delta*(ranf - 0.5);
    double dy = delta*(ranf - 0.5);
    double dz = delta*(ranf - 0.5);

    // update the x, y, z co-ordinates after the move of the k-th particle
    x[k] = x[k] + dx;
    y[k] = y[k] + dy;
    z[k] = z[k] + dz;

    double energy_after_disp = system_total_energy(n, box_length, epsilon, sigma, rc, x, y, z);   // calculating the total system energy after monte carlo move

    // Applying 'Metropolis Method'
    if (energy_after_disp <= energy_before_disp){   // accept the move
      outfile << mc << "   \t    " << energy_after_disp << "   \t   " << k << "   \t   " << "1" << endl;
      // 'k' means which particle is being sampled in the current MC step
      // '1' means accepted and '0' means rejected
      acc = acc + 1.0;   // increasing 'acc' by 1
      // update the co-ordinates i.e. x, y, z values after adding with dx, dy, dz
      // Taking energy_after_disp into energy_before_disp to calculate next loop.
      // Because if the move is accepted then the energy_after_disp will be the energy_before_disp for the next loop
      energy_before_disp = energy_after_disp;
    }
    else {
      double random = double(rand())/double(RAND_MAX);   //generating a random no between 0 and 1
      double energy_diff = energy_after_disp - energy_before_disp;  // calculating energy differnce
      double P = exp(-(energy_diff)/(kb*T));   // Bolzman factor
      if (random < P){   // accept the move
        outfile << mc << "   \t    " << energy_after_disp << "   \t   " << k << "   \t   " << "1" << endl;   // writing in a provided output file
        acc = acc + random;  // increasing 'acc' by random number
        // Taking energy_after_disp into energy_before_disp to calculate next loop.
        energy_before_disp = energy_after_disp;
      }

      else {   // reject the move
        outfile << mc << "   \t    " << energy_before_disp << "   \t   " << k << "   \t   " << "0" << endl;
        // As the move is rejected then x, y,z values are being downgraded to the previous values before starting the current MC loop
        x[k] = x[k] - dx;
        y[k] = y[k] - dy;
        z[k] = z[k] - dz;

      }
    }
    // Here I am willing to write the positions of all particles in each 100 steps with a gap one after another
    if(mc % int(m) == 0){
      for (int i = 0; i < n; i++){
        outffile << x[i] << "   \t   " << y[i] << "   \t    " << z[i] << endl;   // printing the coordinates of all particles
      }
      outffile << "\n";
    }
    // Here I am willing to calculate RDF in each 100 steps
    if(mc % int(m) == 0){
      rdf(n, binNo, binSize, box_length, count, x, y, z);
    }
  }
  // print the average accepatnce ratio
  outFile << acc/double(mcMAX) << endl;

  double r[binNo], g[binNo];

  r[0] = 0.0;
  // Now I have count[1], .............. count[N] total values for all the MC steps. That's why I have to divide final g[i] with total number of MC steps
  for (int i = 0; i < binNo; i++){
    g[i] = count[i];
    g[i] = g[i]/(((4.0*M_PI)/3.0)*(pow((r[i]+binSize),3) - pow(r[i],3))*rho);  // Normalization of g(r)
    outFfile << r[i] << "   \t   " << (g[i]*2*m)/(n*mcMAX) << endl;
    r[i+1] = r[i] + binSize;
  }
  // closing all the files
  outFfile.close();
  outffile.close();
  outFile.close();
  outfile.close();
  return 0;
}
