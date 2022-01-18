//
// Created by cheshire on 09.03.17.
//

#ifndef AWPMD_DFT_AWPMD_DFT_NVGPU_HPP
#define AWPMD_DFT_AWPMD_DFT_NVGPU_HPP

#include "../numeric/memorized_adaptive_mesh_cell.hpp"
#include "../numeric/adaptive_mesh_integrator.hpp"
#include "../awpmd-dft.hpp"
#include <tuple>

template <typename MeshCellType>
class Integrator_nvgpu;

template <typename cell_t, typename packet_t>
struct energy_functor_sic;

class XCEnergy_nvgpu : public XCEnergy {
public:
  XCEnergy_nvgpu(unsigned int packet_number, DFTConfig const &meshConfig, long gpu_id = 0);

  XCResult energy(std::vector<WavePacketInfo> const &wavepackets, bool calc_force) override;

  ~XCEnergy_nvgpu() override;

protected:
  void init_xc_approximation() override ;

  void device_allocation(size_t new_size);

  bool check_texture_object_support();

private:
  ElectronDensity<Mesh_nvgpu> rho;

  std::tuple<double, double>
  regular_mesh_integration(std::vector<std::vector<float>> &derivatives, const vector<WavePacketInfo> &wavepackets,
                           bool calc_force);

  std::tuple<double, double>
  adaptive_mesh_integration(std::vector<std::vector<float>> &derivatives, const vector<WavePacketInfo> &wavepackets,
                            bool calc_force);

  using MeshCell = AdaptiveMeshCell<float>;
  AdaptiveMeshIntegrator<MeshCell, class Integrator_nvgpu<MeshCell>> integrator_;

  using ForceMeshCell = LinkedMeshCell<float>;
  AdaptiveMeshIntegrator<ForceMeshCell, class Integrator_nvgpu<ForceMeshCell>> force_integrator_;

  using MeshhWorldTopo = AdaptiveMeshIntegrator<MeshCell, class Integrator_nvgpu<MeshCell>>::WorldTopology;

protected:
  void compute_forces(std::vector<std::vector<float>> &derivatives, const vector<WavePacketInfo> &wavepackets,
                      energy_functor_sic<MeshCell, PacketType> const &energy_eval);
};

#endif //AWPMD_DFT_AWPMD_DFT_NVGPU_HPP
