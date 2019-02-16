// Radial Distribution Function code written by Sumitava Kundu

#include<iostream>
#include<math.h>
#include<cmath>
#include<cstdlib>
#include<ctime>
#include<iomanip>
#include<fstream>
#include<stdio.h>
#include<sstream>
#include<string>

using namespace std;

// Implementing periodic boundary condition
double Periodic(double u, double a){  // u = distance between two particles
                                      // a = box-length
 if (u <= -a/2){
  u = u + a;
 }
 else if (u > a/2){
  u = u - a;
 }
 return u;
}

int main(){
  int i, j, k, N = 999 ;
  double K, rho = 0.85, L, x[N], y[N], z[N], dx, dy, dz, dR;
  double initial, A[N+3][4];
  K = N/rho;
  L =  cbrt(K);
  srand(time(0));

  string line;
  char str[67];
  ifstream infile;
  infile.open("input_file.xyz");

  // Reading the inputs from input.xyz file
  int counter = 0;
  while((counter < 2) && (infile.good())){  // Redaing the first two lines
      getline (infile, line);
      counter++;
  }
  // Redaing the rest of the lines
  counter = 2;
  while((counter < (N+3)) && (infile.good())){
      for (j = 0; j < 4; j++){
        if (j == 0){   // reading the first column
          infile >> str;
          //cout << str << "   \t   ";
        }
        else {
          infile >> A[counter][j];   // reading the rest 3 columns
        }
      }
      j = 3;
      // reading the 3 columns per row
      x[counter - 2] = A[counter][j-2];
      y[counter - 2] = A[counter][j-1];
      z[counter - 2] = A[counter][j];
      //cout << x[counter - 2] << "   \t   " << y[counter - 2] << "   \t    " << z[counter - 2] << endl;
      counter++;
  }
  infile.close();

  // Creating the bins
  double binSize;
  int binNo = 100;
  double count[binNo], g[binNo];
  double r[binNo];
  binSize = (L*0.5)/double(binNo);   // calculating the bin-size

  // Setting up the initial values of 'count'
  for (i = 0; i < binNo; i++){
      count[i] = 0;
  }

  // Calculating the inter-particle distances
  for (i = 0; i < (N-1); i++){
    for (j = i+1; j < N; j++){
      dx = Periodic((x[i] - x[j]), L);
      dy = Periodic((y[i] - y[j]), L);
      dz = Periodic((z[i] - z[j]), L);
      dR = pow((pow(dx, 2) + pow(dy, 2) + pow(dz, 2)), 0.5);

      initial = 0.0;
      for (k = 0; k < binNo; k++){
        if (dR <= L/2 && dR > initial && dR <= (initial + binSize)){
          count[k] = count[k] + 1.0;
        }
        initial = initial + binSize;   // upgrading the bins
	                               // First initial to (initial+binsize)
				       // second (initial+binsize) to (initial+ 2 * binsize) in this way
      }
    }
  }

  r[0] = 0.0;
  for (i = 0; i < binNo; i++){
      g[i] = count[i];
      g[i] = g[i]/((4.0*3.14159/3.0)*(pow((r[i]+binSize),3) - pow(r[i],3))*rho);   // normalized value of g(r)
      cout << r[i] << "   \t   " << g[i]*2/N << endl;
      r[i+1] = r[i] + binSize;
  }


  return 0;
}
