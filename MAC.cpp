// Monte Carlo Option Price Simulator
// Creator: Mohammed Sharieff
// Just a note, this program needs to be run on a computer with a fast
// processor, it runs terribly on my Microsoft Surface 3

#define _USE_MATH_DEFINES
#include <math.h>
#include <iostream>
#include <iomanip>

using namespace std;

// Computes standard deviation of Option Simulation

double pVar(double * payoffs, double average, int sims){

	double std = 0, a = 0;

	for(int i = 0; i < sims; ++i){
		
		std += pow(payoffs[i] - average, 2);
	}

	std /= (double) sims;
	std = pow(std,0.5);

	return std;
}

// Probability Density Function

double PDF(double Z){

	double R =  0;
	double PI = M_PI;

	double A = 1 / pow(2 * PI, 0.5);
	double B = exp(-pow(Z,2)/2);

	R = A * B;


	return R;
}

// Calculates the Inverse Normal Distribution Function (=NORMINV() in Excel)

double NORMINV(double P){
   
   double R = 0, A = 0, B = 0, C = 0;
   double dX = 0, d3X = 0;
   double delta = 0.0001, error = 0.0005, E = 0;
   double coef = 0, prob = 0;
   int steps = 70, sims = 1000;
   
   B = -10;
   C = -3;
   
   
   for(int i = 0; i < sims; ++i){
      
      B = - 10;
      dX = (C - B) / (double) steps;
      d3X = (C - B) / (3 * (double) steps);
      prob = 0;
      
      for(int j = 0; j < steps + 1; ++j){
         
         if(j%2 == 0){
            coef = 2;
         } else {
            coef = 4;
         }
         if(j == 0 || j == steps){
            coef = 1;
         }
         
         prob += (PDF(B) * coef);
         B += dX;
         
      }
      
      prob *= d3X;
      
      E = abs(prob - P);

      
      if(E > error){
         C += 0.01;
      } else {
         return C;
      }
      
      
   }
  
   return A;
}

// Main method responsible for simulating option price

int main()
{
	double S = 0, K = 0, r = 0, q = 0, v = 0; 
	int t = 0, ttm = 0, sims = 0, type = 0;
	int min = 1, max = 9999;

	// ------------------------------------- USER INPUTS

	S = 105;      // Stock
	K = 100;      // Strike
	r = 0.05;     // Risk-Free Rate
	q = 0.02;     // Dividend
	v = 0.25;      // Volatility
	t = 60;       // Time Till Expiration
	ttm = 200;    // Time steps
	type = 1;     // Option Type (0) Call (1) Put
	sims = 100;   // Number of simulations

	// --------------------------------------------------


	double dX = (double) t / (double) ttm;
	double sDX = pow(dX, 0.5);
	double nS = 0, dS = 0, dW = 0, ranZ = 0;
	double sum = 0, pay = 0;
	double stdev = 0, con = 0;
	double * payoff = new double[sims];
	
   
   // Conducts the option price simulation
   
   cout << endl << "Conductiong Simulation\n" << endl;

   for(int i = 0; i < sims; ++i){
      
      nS = S;
      
      for(int j = 0; j < ttm; ++j){
         
         ranZ = (double)(rand() % (max - min + 1) + min) / 10000;
         dW = NORMINV(ranZ);
         nS = nS + (q * nS * dX) + (v * nS * sDX * dW);

      }
      
      if(type == 0){
         pay = nS - K;
		 if(pay > 0){
			 payoff[i] = pay;
		 } else {
			 payoff[i] = 0;
		 }

      } else {
         pay = K - nS;
		 if(pay > 0){
			 payoff[i] = pay;
		 } else {
			 payoff[i] = 0;
		 }
      }
      
      payoff[i] = payoff[i] * exp(-r * t);
      sum += payoff[i];
	  
	  cout << (i + 1) << "...";
   }
   
   

   sum /= (double) sims;
   stdev = pVar(payoff, sum, sims);
   stdev /= sqrt((double)sims);

   con = 1.96 * stdev;

   double minC = sum - con;
   double maxC = sum + con;

   if(minC < 0){
	   minC = 0;
   }


   // =================================================================================================== RESULTS

   cout << "\n\nOPTIONS PRICE MONTE CARLO SIMULATION\n" << endl;
   cout << "Average Simulated Option Price: $" << sum << endl;
   cout << "Simulation Standard Error: " << stdev << endl;
   cout << "We are 95% confident the option price is between : ($" << minC << ",$" << maxC << ")" << endl;
  
   // ===========================================================================================================

   
	return 0;
}
