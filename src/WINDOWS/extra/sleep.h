#include <thread>
#include <chrono>

inline void usleep(int x){
  std::this_thread::sleep_for(std::chrono::microseconds (x));
}