#include <boost/random.hpp>
#include <iostream>
#include <cstdint>
#include <random>
#include <vector>
#include <string>
#include <math.h>

int nps = 1e2;
std::mt19937_64 rng;
std::uniform_real_distribution<double> distribution(0.0,1.0);

struct state{
  std::string par = "proton";
  double pos =    0;
  double ene = 1700;
  double Z = 1.0;
  double mass = 938.2723;
};

struct material{
  std::string par = "Tungsten";
  double Z = 74.0;
  double M = 183.84;
  double rho = 19.3;
  double Qmin = 6e-04;
};

auto Qmax = [] ( const auto& particle ) {
	      auto gamma = ( particle.ene + particle.mass ) / particle.mass;
	      auto beta  = 1 - ( 1 / gamma );
	      auto qmax = 1.022 * beta * gamma;
	      std::vector<double> quant = { gamma, beta, qmax };
	      return quant;
	    };

auto XS = [&] ( const auto& particle, const auto& material ) {
	    auto amp = 0.1536 * pow(particle.Z, 2) * material.Z * material.rho / material.M;
	    auto quant = Qmax( particle );
	    auto xs = amp * ( 1 / quant[1] ) * ( ( 1 / material.Qmin - 1 / quant[2] )
					      - ( quant[1] / quant[2] )
					      * log( quant[2] / material.Qmin ) );
	    return xs;
	  };

auto pusher = [&] ( auto& particle, const auto& material, const auto& r){
		auto xs = XS( particle, material );
		auto delta =  ( 1 / xs ) * log( 1 / r );
		particle.pos += delta;
		return particle;
	      };

auto decrementer = [&] (auto& particle, const auto& material, const auto& r){
		     auto quant = Qmax( particle );
		     auto delta = pow( ( ( ( 1 - r )/material.Qmin ) + ( r / quant[2] ) ), -1 );
		     particle.ene -= delta;
		     return particle;
		   };

auto add_hist = [] ( const auto& particle, auto& history ){
		  history.push_back( particle );
		  return history;
		};

int main()
{
  rng.seed(123456789);
  state p1;
  std::vector<decltype(p1)> history;
  
  for( int i = 0; i < nps; ++i ){
    state particle;
    material mat;

    auto p = nps / 10;
    if( i % p == 0 ){
      std::cout << "nps = " << i << std::endl;
    }
    
    bool t = true;
    while(t){
      if( particle.pos <= 0.5 && particle.ene >= 200 ){
	// Position Mover
	double r1 = distribution( rng );
	pusher( particle, mat, r1 );
	
	// Energy Decrementer
	double r2 = distribution( rng );
	decrementer( particle, mat, r2 );
      }
      else {
	t = false;
	add_hist( particle, history );
      }
    }
  }

  for( auto t : history ){
    std::cout << t.ene << '\n';
  }
}
