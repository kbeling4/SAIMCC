#include <iostream>
#include <fstream>
#include <cstdint>
#include <random>
#include <vector>
#include <string>
#include <ctime>
#include <math.h>

#include <stlab/concurrency/default_executor.hpp>
#include <stlab/concurrency/future.hpp>

int nps = 1e4;
int num_groups = 4;
std::mt19937_64 rng;
std::uniform_real_distribution<double> distribution(0.0,1.0);

class state {
  double mass = 938.27231;
  double gamma2;

public:
  double Z = 1.0;
  double pos =    0;
  double ene = 1700;
  double beta2;
  double Qmax;

  auto qmax() {
    this->gamma2 = pow( ( this->ene + this->mass ) / this->mass, 2 );
    this->beta2  = 1 - ( 1 / this->gamma2  );
    this->Qmax = 1.022 * this->beta2 * this->gamma2;
  }
};

class timer {
  clock_t start;
  clock_t end;
public:
  std::chrono::duration<double> elapsed_seconds;

  auto get_start() {
    this->start = std::clock();
  }
  auto get_end() {
    this->end = std::clock();
  }
  double  total_time() {
    return double(this->end - this->start) / CLOCKS_PER_SEC;
  }
};

class material{
public:
  double Z = 74.0;
  double M = 183.84;
  double rho = 19.3;
  double Qmin = 6e-04;
};

int main()
{    
  auto XS = [&] ( const auto& particle, const auto& material ) {
	      auto amp = 0.1536 * pow(particle.Z, 2) * material.Z * material.rho / material.M;
	      auto xs = amp * ( 1 / particle.beta2 ) * ( ( 1/material.Qmin - 1/particle.Qmax )
							 - ( particle.beta2/particle.Qmax )
							 * log( particle.Qmax/material.Qmin ) );
	      return xs; };
  
  auto pusher = [&] ( auto& particle, const auto& material, const auto& r){
		  particle.qmax();
		  auto xs = XS( particle, material );
		  auto delta =  ( 1 / xs ) * log( 1 / r );
		  particle.pos += delta;
		  return particle; };
  
  auto decrementer = [&] ( auto& particle, const auto& material, const auto& r){
		       particle.qmax();
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
		  auto   Q = 1e-3;
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
		particle.qmax();
		auto delta = Newton( particle, material, r );
		return particle.ene -= delta; };
  
  
  auto add_hist = [] ( const auto& particle, auto& history ){
		    history.push_back( particle );
		    return history; };

  auto printer = [] (const auto& history) {
		   std::ofstream myfile ("output.txt");
		   if (myfile.is_open())
		     {
		       for( auto t : history ){
			 myfile << "( " << t.ene << ", " << t.pos << " ) ";
		       }
		     }
		   myfile.close(); };

  state p1;
  std::vector<decltype(p1)> history;
  timer t;
  auto num_part = nps / num_groups;
  auto looper = [&] () {
		  t.get_start();
		  rng.seed(123456789);

		  for( int i = 0; i < num_part; ++i){
		    state particle;
		    material mat;

		    auto p = num_part / 10;
		    if( i % p == 0 ){
		      std::cout << "nps = " << i << std::endl;
		    }
		    
		    while(true){
		      if( particle.pos <= 0.5 && particle.ene >= 1e-2 ){
			// Position Mover
			double r1 = distribution( rng );
			pusher( particle, mat, r1 );
			
			// Energy Decrementer
			double r2 = distribution( rng );
			DecE( particle, mat, r2 );
		      }
		      else {
			add_hist( particle, history );
			break;
		      }
		    }
		  }
		  return history;
		};

  // Run Problem ----------------------------------
  looper();

  // Output
  printer( history );

  t.get_end();
  std::cout << "---------------------------------" << std::endl;
  std::cout << "Total time: " << t.total_time() << std::endl;
  return 0;
}
