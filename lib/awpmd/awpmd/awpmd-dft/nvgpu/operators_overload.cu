#include "operators_overload.hpp"

__host__ __device__ uint3 operator+ (const uint3 &a, const uint3 &b){
    return make_uint3(a.x + b.x,    a.y + b.y,    a.z + b.z);
}
__host__ __device__ uint3 operator* (const uint3 &a, const uint3 &b){
    return make_uint3(a.x * b.x,    a.y * b.y,    a.z * b.z);
}
__host__ __device__ uint3 operator* (const uint3 &a, const dim3 &b){
    return make_uint3(a.x * b.x,    a.y * b.y,    a.z * b.z);
}
__host__ __device__ uint3 operator* (const dim3 &a, const uint3 &b){
    return make_uint3(a.x * b.x,    a.y * b.y,    a.z * b.z);
}