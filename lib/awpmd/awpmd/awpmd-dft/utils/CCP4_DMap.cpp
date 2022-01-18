//
// Created by yalavrinenko on 27.08.18.
//

#include "CCP4_DMap.h"
#include "../awpmd-dft.hpp"
#include "../DataTypes.hpp"

void ccp4_dmap_writer::write_header(Mesh const &m, DFTConfig const &config, std::string label) {
  ccp4_header header;

  header.NC = m.m_width;
  header.NR = m.m_height;
  header.NS = m.m_depth;
  std::tie(header.NX, header.NY, header.NZ) = std::tie(config.mesh_size.size.as_struct.x,
                                                       config.mesh_size.size.as_struct.y,
                                                       config.mesh_size.size.as_struct.z);
  double space_size[] = {1.0, 1.0, 1.0};
  for (auto i : {0, 1, 2})
    space_size[i] = config.mesh_fin.size.as_array[i] - config.mesh_start.size.as_array[i];

  std::tie(header.X, header.Y, header.Z) = std::tie(space_size[0], space_size[1], space_size[2]);

  header.NLABL = (int32_t) label.size();
  std::copy(label.begin(), label.end(), &header.LABEL[0]);

  float dmin = std::numeric_limits<float>::max(), dmax = std::numeric_limits<float>::min(), dmean = 0;
  m.apply([&dmin, &dmax, &dmean](float v) {
    if (dmin > v)
      dmin = v;
    if (dmax < v)
      dmax = v;
    dmean += v;
  });

  auto integral = dmean * 1.0 / (m.m_depth * m.m_height * m.m_width);

  dmean /= m.count();

  std::tie(header.AMIN, header.AMAX, header.AMEAN) = std::tie(dmin, dmax, dmean);

  m_stream.write((char const *) &header, sizeof(header));

  Logger::Warning(label, "info:\n\tMin:", dmin, "Max:", dmax, "Integral", integral);
}

void ccp4_dmap_writer::write(Mesh const &m, DFTConfig const &config, std::string label) {
  m_stream.open(m_path, std::fstream::out | std::fstream::binary);
  this->write_header(m, config, label);

  std::vector<float> buffer;
  buffer.reserve(m.count());
  m.apply([&buffer](float v) {
    buffer.emplace_back(v);
  });

  m_stream.write((char const *) buffer.data(), sizeof(float) * buffer.size());
  m_stream.close();
}
