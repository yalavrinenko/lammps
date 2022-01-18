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
 * $Log: rebuild_wp.h,v $
 * Revision 1.1  2013/03/04 16:46:17  valuev
 * added rebuild_wp interface
 *
*******************************************************************************/

# ifndef REBUILD_WP_H
# define REBUILD_WP_H

# include "wpmd_split.h"

/// Interface function for Andreas to rebuid system with the Matching Pursuit method.
bool rebuild_wavepackets(AWPMD_split &system /*, add you parameters here */);


# endif