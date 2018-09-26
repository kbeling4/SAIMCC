#include <boost/random.hpp>
#include <iostream>
#include <cstdint>
#include <random>
#include <vector>

struct state{
  double pos =    0;
  double ene = 1700;
};

auto pusher( auto& particle, auto& r ){
  particle.pos = particle.pos + r;
  return particle;
}


int main()
{
  int nps = 1e2;
  std::mt19937_64 rng;
  rng.seed(123456789);
  std::uniform_real_distribution<double> distribution(0.0,1.0);

  for( int i = 0; i < nps; ++i ){
    state particle;
    bool t = true;
    while(t){
      if( particle.pos <= 2.0 && particle.ene >= 200 ){
	double r = distribution( rng );
	pusher( particle, r );
      }
      else {
	t = false;
	std::cout << particle.pos << std::endl;
      }
    }
  }
  
  // int nps = 1e2;
  // boost::uniform_real<> uni_dist(0,1);
  // std::vector<long double> range(nps);
  // for( auto i : range ){
  //   i = gen();
  //   std::cout << i << '\n';
  // }
}
