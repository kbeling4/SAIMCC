#include <chrono>

class Timer {
  clock_t    start;
  clock_t      end;
  std::chrono::duration<double> duration;
  
public:
  auto get_start() {
    auto s = std::chrono::system_clock::now();
    this->start = s;
  }

};
