#include <iostream>
#include <fstream>
#include <cstdint>
#include <random>
#include <vector>
#include <string>
#include <math.h>

#include "State.hpp"

int nps = 5e2;
std::mt19937_64 rng;
std::uniform_real_distribution<double> distribution(0.0,1.0);

struct material{
  //std::string name = "Tungsten";
  double Z = 6.0;
  double M = 12.0107;
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
  
  auto pusher = [&] ( auto& particle, const auto& material, const auto& r){
		  particle.get_qmax();
		  auto xs = XS( particle, material );
		  auto delta =  ( 1 / xs ) * log( 1 / r );
		  particle.position += delta;
		  return particle; };
  
  auto decrementer = [&] ( auto& particle, const auto& material, const auto& r){
		       particle.get_qmax();
		       auto delta = pow( ( ( ( 1-r )/material.Qmin )+( r/particle.Qmax ) ),-1 );
		       particle.ene -= delta;
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
		  //std::cout << Q << " ";
		  return Q; };
  
  auto DecE = [&] ( auto& particle, const auto& material, const auto& r ) {
		particle.get_qmax();
		auto delta = Newton( particle, material, r );
		return particle.energy -= delta; };
  
  
  auto add_hist = [] ( const auto& particle, auto& history ){
		    history.push_back( particle );
		    return history; };

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
			  tally[i] += state.energy;
			  break;
			}
			++state.tal;
		      }
		      return tally;
		    };

  auto normalizer = [] ( const auto& nps, auto& tally, const auto& bins ){
		      state p;
		      for( auto i = 0; i < tally.size(); ++i ){
			tally[i] = tally[i] / (nps);
		      }
		      return tally;
		    };

  auto loop = [&] ( const auto& nps, const auto& bins, auto& tally ) {

		for( int i = 0; i < nps; ++i ){
		  state particle;
		  material mat;

		  auto p = nps / 10;
		  if( i % p == 0 ){
		    std::cout << "nps = " << i << std::endl;
		  }
		  
		  while(true){
		    auto p_energy = particle.energy;
		    //add_hist( particle, history );
		    if( particle.position <= 30.0 && particle.energy >= 3e-1 ){
		      // Position Mover
		      double r1 = distribution( rng );
		      pusher( particle, mat, r1 );
		      
		      // Energy Decrementer
		      double r2 = distribution( rng );
		      DecE( particle, mat, r2 );

		      // Energy Tally
		      eTallier( bins, particle, p_energy, tally );
		    }
		    else {
		      endTallier( bins, particle, tally );
		      break;
		    }
		  }
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

  
  // Run problem -------------------------------------------------------------------
  auto start = std::chrono::system_clock::now();

  state p1;
  std::vector<decltype(p1)> history;
  rng.seed(123456789);

  auto bins = linspace( 0.0, 30.0, 1000 );
  std::vector<double> tally( bins.size() - 1, 0.0 );
  
  loop( nps, bins, tally );
  normalizer( nps, tally, bins );
  printer( tally, bins );

  auto end = std::chrono::system_clock::now();
  std::chrono::duration<double> elapsed_seconds = end-start;
  std::cout << "---------------------------------" << std::endl;
  std::cout << "Total time: " << elapsed_seconds.count() << std::endl;
  return 0;
}
