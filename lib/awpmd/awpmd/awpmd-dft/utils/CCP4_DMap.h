//
// Created by yalavrinenko on 27.08.18.
//

#ifndef MDUTILS_LIB_CCP4_DMAP_H
#define MDUTILS_LIB_CCP4_DMAP_H

#include "SpaceMesh.hpp"
#include <fstream>
#include <utility>

struct DFTConfig;
//description at http://www.ccp4.ac.uk/html/maplib.html#description
struct ccp4_header{
    int32_t NC, NR, NS;
    int32_t MODE = 2;
    int32_t NCSTART = 0;
    int32_t NRSTART = 0;
    int32_t NSSTART = 0;
    int32_t NX;
    int32_t NY;
    int32_t NZ;
    float X;
    float Y;
    float Z;
    float Alpha = 0, Beta = 0, Gamma = 0;
    int32_t MAPC = 1;
    int32_t MAPR = 2;
    int32_t MAPS = 3;
    float AMIN;
    float AMAX;
    float AMEAN;
    int32_t ISPG = 0;
    int32_t NSYMBT = 0;
    int32_t LSKFLG = 0;
    int32_t SKWMAT[9]={0};
    int32_t SKWTRN[3]={0};
    int32_t FUTURE[15]={0};
    char MAP[4] = "MAP";
    int32_t MACHST = 0x1111;
    int32_t ARMS = 0;
    int32_t NLABL;
    uint32_t LABEL[200];
};

class ccp4_dmap_writer{
public:
    explicit ccp4_dmap_writer(std::string path):
        m_path(std::move(path)){
    }

    void write(Mesh const& m, DFTConfig const &config, std::string label = "");
private:
    void write_header(Mesh const &m, DFTConfig const &config, std::string label = "");

    std::string m_path;
    std::ofstream m_stream;
};

#endif //MDUTILS_LIB_CCP4_DMAP_H
