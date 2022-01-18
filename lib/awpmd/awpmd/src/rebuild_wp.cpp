/*s***************************************************************************
 *
 *   Copyright (c), Ilya Valuev 2013        All Rights Reserved.
 *
 *   Author	: Ilya Valuev, MIPT, Moscow, Russia
 *
 *   Project	: AWPMD
 *
 *****************************************************************************/
  
/*s****************************************************************************
 * $Log: rebuild_wp.cpp,v $
 * Revision 1.1  2013/03/04 16:46:17  valuev
 * added rebuild_wp interface
 *
*******************************************************************************/


# include "rebuild_wp.h"


bool rebuild_wavepackets(AWPMD_split &system /*, add you parameters here */){
  int nel[2] = { system.get_electron_number(0), system.get_electron_number(1) };
  for(int spin = 0; spin < 2 ; spin++ ) {  // loop spins
    int wpid = 0; // wpid runs through all electrons of the given spin
    for(int el = 0; el < nel[spin] ; el++ ){
      int n_gaussians = system.get_split_number(spin, el); // detect how many wps belong to current electron
      WavePacket *wp = new WavePacket[n_gaussians]; // wave packet array for current electron, may use std:vector instead
      complex<double> *coeffs = new complex<double>[n_gaussians]; // amplitudes array for current electron
      int wpid0 = wpid; 
      for(int split = 0; split < n_gaussians; split++) { // filling array
        system.get_wavepacket(spin, wpid, coeffs[split], wp[split]); // get wave packet and amplitudes from the system
        wpid++;
      }
      
      // call procedure to rebuild wave function
      // apply_matching_pursuit(n_gaussians, wp, coeffs /* add more parameters here */)
      
      for(int split = 0; split < n_gaussians; split++) { // recompose the system: number of gaussians remains the same
        system.set_wavepacket(spin, wpid0, coeffs[split], wp[split]); // set wave packet and amplitudes 
        wpid0++;
      }
      delete [] wp;
      delete [] coeffs;
    } // next electron
  }// next spin
  return true;
}

