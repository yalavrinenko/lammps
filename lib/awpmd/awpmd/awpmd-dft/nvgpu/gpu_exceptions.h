//
// Created by yalavrinenko on 01.10.2019.
//

#ifndef DERIVS_GPU_API_H
#define DERIVS_GPU_API_H
#include <cuda_runtime_api.h>

namespace {
  thread_local std::string error_text;
}

class gpu_runtime_exception : public std::exception {
public:
  gpu_runtime_exception(char const *what, char const *descr, int code) :
      reason(what), m_code(code), cuda_description(descr) {
  }

  char const *what() const noexcept override {
    ::error_text = reason + ":" + std::to_string(m_code) + ":" + cuda_description;
    return ::error_text.c_str();
  }

  int code() const {
    return m_code;
  }

private:
  std::string reason;
  int m_code;
  std::string cuda_description;
};

#define SAFECALL(Operation, ErrorMessage) \
{\
    auto cuError = Operation; \
    if (cuError != cudaSuccess) {\
        Logger::Error(ErrorMessage, cudaGetErrorString(cuError), " Code: ", cuError);\
        throw gpu_runtime_exception(ErrorMessage, cudaGetErrorString(cuError), cuError); \
    }\
}


#endif //DERIVS_GPU_API_H
