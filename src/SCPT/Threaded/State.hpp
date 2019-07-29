#include <math.h>
#include <string>

class state{

  double mass = 5633.51382;
  double gamma2;

public:
  std::string name = "Li-6";
  double position =    0;
  double energy;
  double Z = 3.0;
  double beta2;
  double Qmax;
  double tal = 0;

  auto get_qmax() {
    this->gamma2 = pow( ( this->energy + this->mass ) / this->mass, 2 );
    this->beta2  = 1 - ( 1 / this->gamma2  );
    this->Qmax = 1.022 * this->beta2 * this->gamma2; 
  }
};

