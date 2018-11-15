#include <iostream>
#include <fstream>
#include <cstdint>
#include <random>
#include <vector>
#include <string>
#include <math.h>

#include <omp.h>

#include "State.hpp"

#define NUM_THREADS 4

int nps = 1e3;
std::mt19937_64 rng;
std::uniform_real_distribution<double> distribution(0.0,1.0);

struct material{
  std::string name = "Graphite";
  double Z = 6.0;
  double M = 12.011;
  double rho = 2.266;
  double Qmin = 7.80e-05;
};


template <typename T>
std::vector<T> linspace(T a, T b, size_t N) {
  T h = (b - a) / static_cast<T>(N-1);
  std::vector<T> xs(N);
  typename std::vector<T>::iterator x;
  T val;
  for (x = xs.begin(), val = a; x != xs.end(); ++x, val += h)
    *x = val;
  return xs;
} 

int main()
{
  auto XS = [&] ( const auto& particle, const auto& material ) {
	      auto amp = 0.1536 * pow(particle.Z, 2) * material.Z * material.rho / material.M;
	      auto xs = amp * ( 1 / particle.beta2 ) * ( ( 1/material.Qmin - 1/particle.Qmax )
							 - ( particle.beta2 / particle.Qmax )
							 * log( particle.Qmax/material.Qmin ) );
	      return xs; };
  
  auto pusher = [&] ( auto& particle, const auto& material, auto& gen){
		  particle.get_qmax();
		  auto xs = XS( particle, material );
		  auto delta =  ( 1 / xs ) * log( 1 / distribution( gen ) );
		  particle.position += delta;
		  return particle; };
  
  auto dFDQ = [&] (const auto& particle, const auto& material, auto& Q ) {
		auto xs = ( 1/material.Qmin - 1/particle.Qmax ) - (particle.beta2/particle.Qmax)
		  * log(particle.Qmax/material.Qmin);
		return ( 1/(Q*Q) - (particle.beta2/particle.Qmax)*(1/Q) ) / xs; };
  
  auto FQ = [&] ( const auto& particle, const auto& material, auto& Q, const auto& rand ) {
	      auto xs = ( 1/material.Qmin - 1/particle.Qmax ) - (particle.beta2/particle.Qmax)
		* log(particle.Qmax/material.Qmin);
	      return (( 1/material.Qmin - 1/Q ) - (particle.beta2/particle.Qmax)
		      * log( Q / material.Qmin ) ) / xs - rand; };
  
  auto Newton = [&] (const auto& particle, const auto& material, const auto& rand ) {
		  auto eps = 1e-3;
		  auto   Q = 1e-4;
		  auto fValue = FQ( particle, material, Q, rand );
		  auto itCount = 0;
		  while( abs( fValue ) > eps ) {
		    Q -= fValue / dFDQ( particle, material, Q );
		    fValue = FQ( particle, material, Q, rand );
		    ++ itCount; 
		  }
		  return Q; };
  
  auto DecE = [&] ( auto& particle, const auto& material, auto& gen ) {
		particle.get_qmax();
		auto delta = Newton( particle, material, distribution( gen ) );
	        particle.energy -= delta;
		return particle; };
  
  auto eTallier = [] ( const auto& bins, auto& state, const auto& p_energy, auto& tally ){
		    for( auto i = state.tal; i < tally.size(); ++i ){
		      if( state.position < bins[i] ) {
			tally[i] += p_energy - state.energy;
			break;
		      }
		      ++state.tal;
		    }
		    return tally;
		  };

  auto endTallier = [] ( const auto& bins, auto& state, auto& tally ){
		      for( auto i = state.tal; i < tally.size(); ++i ){
			if( state.position < bins[i] ) {
			  tally[i] += 1;
			  break;
			}
		      }
		      return tally;
		    };

  auto exitTallier = [] ( const auto& bins, auto& state, auto& tally ){
		      for( auto i = 0; i < tally.size(); ++i ){
			if( state.energy < bins[i] ) {
			  tally[i] += 1;
			  break;
			}
		      }
		      return tally;
		    };

  auto normalizer = [] ( const auto& nps, auto& tally, const auto& bins ){
		      for( auto i = 0; i < tally.size(); ++i ){
			tally[i] = tally[i] / (nps);
		      }
		      return tally;
		    };

  auto printer = [] ( const auto& tally, const auto& bins ) {
		   std::ofstream myfile ("output.txt");
		   if (myfile.is_open())
		     {
		       for( auto i = 0; i < tally.size(); ++i ){
			 auto pos = ( bins[i+1] + bins[i] ) / 2;
			 myfile << pos << " " <<  tally[i] << "\n";
		       }
		     }
		   myfile.close();
		 };

  auto physics = [&] ( auto& particle, const auto& mat ) {
		   do{ 
		     //auto p_energy = particle.energy;
		     // Position Mover
		     pusher( particle, mat, rng );
		     
		     // Energy Decrementer
		     DecE( particle, mat, rng );
	
		     // Energy Tally
		     //#pragma critical	  
		     // eTallier( bins, particle, p_energy, tally );
		   } while( particle.position <= 28.0 && particle.energy >= 1e-3 ); 
		 };
		   
  
  // Run problem -------------------------------------------------------------------
  double start_time, run_time;
  start_time = omp_get_wtime();

  rng.seed(123456789);

  auto bins = linspace( 0.0, 75.0, 200 );
  std::vector<double> tally( bins.size() - 1, 0.0 );
  
  // Loop -------------------------------------------------------------------------
#pragma omp parallel
  {
    //    int id = omp_get_thread_num();
    //    int numthreads = omp_get_num_threads();
    auto num = nps / 10;

#pragma omp for schedule( dynamic, 100 )
    for( int i = 0; i < nps; ++i ){
      state particle;
      material mat;
      physics( particle, mat );

      if( i % num == 0 ){
	std::cout << "nps = " << i << std::endl;
      }
#pragma omp critical
      exitTallier( bins, particle, tally );
    }
  }
  // -------------------------------------------------------------------------------

  //normalizer( nps, tally, bins );
  printer( tally, bins );
  std::cout << "nps = " << nps << std::endl;
  run_time = omp_get_wtime() - start_time;
  std::cout << "---------------------------------" << std::endl;
  std::cout << "Total time: " << run_time << std::endl;
  return 0;
}
