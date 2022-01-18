//
// Created by cheshire on 20.04.17.
//
#include "SpaceMesh.hpp"
#include "../DataTypes.hpp"
#include <cstring>
#include <iostream>

unsigned int inline Mesh::address(unsigned int i, unsigned int j, unsigned int k) const{
    return k*(m_width * m_height) + j * m_width + i;
}

Mesh::Mesh(unsigned int width, unsigned int height, unsigned int depth): grid_size(width * height * depth), m_width(width), m_height(height),
                                                                         m_depth(depth)
{
    m_host_storage = new double[grid_size];
    m_device_storage = m_host_storage;

    Logger::Info("Create CPU Mesh. Size:", width,'x',height,'x',depth);
}

Mesh::Mesh(const Mesh& m):Mesh(m.m_width, m.m_height, m.m_depth) {}

void Mesh::clear(){
    memset(m_device_storage, 0, sizeof(float) * grid_size);
}

void Mesh::clear(int value){
    memset(m_device_storage, value, sizeof(float) * grid_size);
}

double& Mesh::at(unsigned int i, unsigned int j, unsigned int k){
    return m_device_storage[this->address(i,j,k)];
}

double& Mesh::ath(unsigned int i, unsigned int j, unsigned int k){
    return m_host_storage[this->address(i,j,k)];
}

double Mesh::ath(unsigned int i, unsigned int j, unsigned int k) const{
    return m_host_storage[this->address(i,j,k)];
}

void Mesh::copy_memory(int direction){
}

void Mesh::Free(){
    if (m_host_storage != nullptr)
        delete[] m_host_storage;
}

Mesh::~Mesh(){
}

double Mesh::reduction() {
  auto size = this->grid_size;
  do{
    if (size % 2 != 0)
      m_host_storage[size-2] += m_host_storage[size-1];

    size /= 2;

    for (auto i = 0u; i < size; ++i)
      m_host_storage[i] += m_host_storage[i + size];
  }
  while (size > 1);

  return m_host_storage[0];
}
