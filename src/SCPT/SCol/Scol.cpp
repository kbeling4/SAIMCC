#include <boost/random.hpp>
#include <iostream>
#include <fstream>
#include <cstdint>
#include <random>
#include <vector>
#include <string>
#include <math.h>

int nps = 1e5;
std::mt19937_64 rng;
std::uniform_real_distribution<double> distribution(0.0,1.0);

struct state{
  std::string par = "proton";
  double pos =    0;
  double ene = 1700;
  double Z = 1.0;
  double mass = 938.27231;
  double gamma2;
  double beta2;
  double Qmax;
};

struct material{
  std::string par = "Tungsten";
  double Z = 74.0;
  double M = 183.84;
  double rho = 19.3;
  double Qmin = 6e-04;
};

auto Qmax = [] ( auto& particle ) {
	      particle.gamma2 = pow( ( particle.ene + particle.mass ) / particle.mass, 2 );
	      particle.beta2  = 1 - ( 1 / particle.gamma2  );
	      particle.Qmax = 1.022 * particle.beta2 * particle.gamma2;
	      return particle; };

auto XS = [&] ( auto& particle, const auto& material ) {
	    auto amp = 0.1536 * pow(particle.Z, 2) * material.Z * material.rho / material.M;
	    Qmax( particle );
	    auto xs = amp * ( 1 / particle.beta2 ) * ( ( 1 / material.Qmin - 1 / particle.Qmax )
					      - ( particle.beta2 / particle.Qmax )
					      * log( particle.Qmax / material.Qmin ) );
	    return xs; };

auto pusher = [&] ( auto& particle, const auto& material, const auto& r){
		auto xs = XS( particle, material );
		auto delta =  ( 1 / xs ) * log( 1 / r );
		particle.pos += delta;
		return particle; };

auto decrementer = [&] ( auto& particle, const auto& material, const auto& r){
		     Qmax( particle );
		     auto delta = pow( ( ( ( 1 - r )/material.Qmin )+( r/particle.Qmax ) ),-1 );
		     particle.ene -= delta;
		     return particle; };

auto dFDQ = [&] ( auto& particle, const auto& material, auto& Q ) {
	      Qmax( particle );
	      auto xs = ( 1/material.Qmin - 1/particle.Qmax ) - (particle.beta2/particle.Qmax)
		* log(particle.Qmax/material.Qmin);
	      return ( 1/(Q*Q) - (particle.beta2/particle.Qmax)*(1/Q) ) / xs; };

auto FQ = [&] ( auto& particle, const auto& material, auto& Q, const auto& rand ) {
	    Qmax( particle );
	    auto xs = ( 1/material.Qmin - 1/particle.Qmax ) - (particle.beta2/particle.Qmax)
	      * log(particle.Qmax/material.Qmin);
	    return (( 1/material.Qmin - 1/Q ) - (particle.beta2/particle.Qmax)
		    * log( Q / material.Qmin ) ) / xs - rand; };

auto Newton = [&] ( auto& particle, const auto& material, const auto& rand ) {
		auto eps = 1e-3;
		auto   Q = 1e-3;
		auto fValue = FQ( particle, material, Q, rand );
		auto itCount = 0;
		while( abs( fValue ) > eps ) {
		  Q -= fValue / dFDQ( particle, material, Q );
		  fValue = FQ( particle, material, Q, rand );
		  ++ itCount; 
		}
		//std::cout << Q << std::endl;
		return Q; };

auto DecE = [&] ( auto& particle, const auto& material, const auto& r ) {
	      auto delta = Newton( particle, material, r );
	      return particle.ene -= delta; };
		

auto add_hist = [] ( const auto& particle, auto& history ){
		  history.push_back( particle );
		  return history; };

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
      if( particle.pos <= 0.5 && particle.ene >= 1e-2 ){
	// Position Mover
	double r1 = distribution( rng );
	pusher( particle, mat, r1 );
	
	// Energy Decrementer
	double r2 = distribution( rng );
	DecE( particle, mat, r2 );
      }
      else {
	t = false;
	add_hist( particle, history );
      }
    }
  }

  std::ofstream myfile ("output.txt");
  if (myfile.is_open())
  {
    for( auto t : history ){
      myfile << t.ene << "\n";
    }
  }
  myfile.close();
  return 0;
}
