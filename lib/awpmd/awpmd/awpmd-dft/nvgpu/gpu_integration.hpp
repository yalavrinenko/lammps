//
// Created by yalavrinenko on 09.10.2019.
//

#ifndef ADAPTIVE_TEST_GPU_INTEGRATION_HPP
#define ADAPTIVE_TEST_GPU_INTEGRATION_HPP

#include <vector>
#include <thrust/device_vector.h>
#include <thrust/execution_policy.h>

template <typename data_t>
struct memory_adapter{
protected:
  template <typename iterator_t>
  void allocate(iterator_t const &begin, iterator_t const &end) {
    memory_ = std::move(thrust::device_vector<data_t>{begin, end});
  }

  thrust::device_vector<data_t> memory_;

public:
  memory_adapter() = default;

  template <typename iterator_t>
  auto begin(iterator_t hbegin, iterator_t hend) -> decltype(memory_.begin()){
    if (memory_.empty())
      allocate(hbegin, hend);
    return memory_.begin();
  }

  template <typename iterator_t>
  auto end(iterator_t hbegin, iterator_t hend) -> decltype(memory_.end()){
    if (memory_.empty())
      allocate(hbegin, hend);
    return memory_.end();
  }

  auto raw_begin(std::vector<data_t> const &host_vector) -> decltype(thrust::raw_pointer_cast(&memory_[0])){
    if (memory_.empty())
      allocate(host_vector.begin(), host_vector.end());
    return thrust::raw_pointer_cast(&memory_[0]);
  }

  auto raw_end(std::vector<data_t> const &host_vector) -> decltype(thrust::raw_pointer_cast(&memory_[memory_.size()])){
    if (memory_.empty())
      allocate(host_vector.begin(), host_vector.end());
    return thrust::raw_pointer_cast(&memory_[memory_.size()]);
  }
};

template<typename output_t, typename arg_t>
struct reduction_functor{
  __device__ output_t operator() (output_t sum, arg_t value){
    return sum + value;
  }
};

template <typename T>
struct linear_index_to_row_index : public thrust::unary_function<T,T>
{
  T C; // number of columns
  __host__ __device__ explicit linear_index_to_row_index(T C) : C(C) {}

  __host__ __device__ T operator()(T i) { return i / C;  }
};

template <typename cell_t>
class Integrator_nvgpu {
public:
  template <typename output_t, typename PacketType, typename FunctionType, typename iterator_t>
  static output_t integrate(iterator_t cbegin, iterator_t cend, std::vector<PacketType> const &packets, FunctionType &ifunc) {
    thrust::device_vector<PacketType> dev_packets(packets.begin(), packets.end());
    thrust::device_vector<cell_t> dev_cells(cbegin, cend);
    thrust::device_vector<output_t> I(std::distance(cbegin, cend), output_t{});

    ifunc.packets = thrust::raw_pointer_cast(&dev_packets[0]);
    ifunc.count = packets.size();

    thrust::transform(thrust::device, dev_cells.begin(), dev_cells.end(), I.begin(), ifunc);

    return  thrust::reduce(thrust::device, I.begin(), I.end(), output_t{});
  }

  template <typename output_t, typename PacketType, typename FunctionType, typename iterator_t>
  static output_t integrate(iterator_t cbegin, iterator_t cend, std::vector<PacketType> const &packets, size_t batch_size, FunctionType &ifunc) {
    using datatype = typename output_t::value_type;

    thrust::device_vector<PacketType> dev_packets(packets.begin(), packets.end());
    thrust::device_vector<cell_t> dev_cells(cbegin, cend);
    thrust::device_vector<datatype> I(std::distance(cbegin, cend), datatype{});

    ifunc.packets = thrust::raw_pointer_cast(&dev_packets[0]);
    ifunc.count = packets.size();

    thrust::transform(thrust::device, dev_cells.begin(), dev_cells.end(), I.begin(), ifunc);

    auto const &N = packets.size();
    thrust::device_vector<unsigned int> out_keys(N);
    thrust::device_vector<datatype> out_value(N);

    auto key_begin = thrust::make_transform_iterator(thrust::counting_iterator<unsigned int>(0), linear_index_to_row_index<unsigned int>(batch_size));
    auto reduce_out = thrust::reduce_by_key(thrust::device, key_begin, key_begin + std::distance(cbegin, cend), I.begin(), out_keys.begin(), out_value.begin());

    std::vector<unsigned int> keys(N);
    std::vector<datatype> output(N);

    thrust::copy(out_keys.begin(), out_keys.end(), keys.begin());
    thrust::copy(out_value.begin(), out_value.end(), output.begin());
    return output;
  }

  template <typename FunctionType, typename iterator_t>
  static void for_each(iterator_t cbegin, iterator_t cend, memory_adapter<cell_t> &mem_adapter, FunctionType &func){
    thrust::for_each(thrust::device, mem_adapter.begin(cbegin, cend), mem_adapter.end(cbegin, cend), func);
  }

  template <typename output_t, typename FunctionType, typename iterator_t>
  static output_t integrate(iterator_t cbegin, iterator_t cend, memory_adapter<cell_t> &mem_adapter, FunctionType &ifunc){
    thrust::device_vector<output_t> I(std::distance(cbegin, cend), output_t{});
    thrust::transform(thrust::device, mem_adapter.begin(cbegin, cend), mem_adapter.end(cbegin, cend), I.begin(), ifunc);
    return  thrust::reduce(thrust::device, I.begin(), I.end(), output_t{});
//    return thrust::transform_reduce(thrust::device, mem_adapter.begin(cells), mem_adapter.end(cells), ifunc, output_t{}, reduction_functor<output_t, output_t>{});
  }

  template <typename output_t, typename FunctionType, typename iterator_t>
  static output_t integrate(iterator_t cbegin, iterator_t cend, memory_adapter<cell_t> &mem_adapter, size_t batch_size, FunctionType &ifunc){
    using datatype = typename output_t::value_type;

    thrust::device_vector<datatype> I(std::distance(cbegin, cend), output_t{});
    thrust::transform(thrust::device, mem_adapter.begin(cbegin, cend), mem_adapter.end(cbegin, cend), I.begin(), ifunc);

    auto N = std::distance(cbegin, cend) / batch_size;
    thrust::device_vector<unsigned int> out_keys(N);
    thrust::device_vector<datatype> out_value(N);

    auto key_begin = thrust::make_transform_iterator(thrust::counting_iterator<unsigned int>(0), linear_index_to_row_index<unsigned int>(batch_size));
    auto reduce_out = thrust::reduce_by_key(thrust::device, key_begin, key_begin + std::distance(cbegin, cend), I.begin(), out_keys.begin(), out_value.begin());

    std::vector<unsigned int> keys(N);
    std::vector<datatype> output(N);

    thrust::copy(out_keys.begin(), out_keys.end(), keys.begin());
    thrust::copy(out_value.begin(), out_value.end(), output.begin());
    return output;
  }
};

#endif // ADAPTIVE_TEST_GPU_INTEGRATION_HPP
