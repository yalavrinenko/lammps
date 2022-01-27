#include <chrono>
#include <thread>

inline void usleep(int x) {
  std::this_thread::sleep_for(std::chrono::microseconds(x));
}

//void usleep(int x);
