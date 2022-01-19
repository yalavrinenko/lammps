# ifndef WPMD_SPLIT_H
# define WPMD_SPLIT_H

/** @file wpmd_split.h
    @brief Representation of electrons by multiple wave packets within WPMD */

/*s****************************************************************************
 * $Log: wpmd_split.h,v $
 * Revision 1.60  2015/10/20 15:01:13  valuev
 * switched to width velocity
 *
 * Revision 1.59  2015/02/10 14:08:09  valuev
 * added step function
 *
 * Revision 1.58  2014/09/26 16:54:20  valuev
 * syncronizing with kintech svn
 *
 * Revision 1.57  2014/07/28 10:45:35  valuev
 * added partial WP integration
 *
 * Revision 1.56  2014/07/28 08:39:49  valuev
 * discrete fermi distributions
 *
 * Revision 1.55  2014/07/21 13:33:17  valuev
 * quantum MC, density for UHF case
 *
 * Revision 1.54  2014/07/18 18:49:21  valuev
 * quantum MC
 *
 * Revision 1.53  2014/07/17 11:57:32  valuev
 * more correct oscillator partition
 *
 * Revision 1.52  2014/07/16 10:32:11  valuev
 * quantum sums for AWPMC
 *
 * Revision 1.51  2014/07/09 15:27:51  valuev
 * tuning ii interaction for MC
 *
 * Revision 1.50  2014/07/08 17:07:21  valuev
 * prepared confined MC simulation with AWPMD nonsplit norm matrix
 *
 * Revision 1.49  2014/07/03 16:58:42  valuev
 * added overlap matrix det log measurement
 *
 * Revision 1.48  2014/07/03 13:32:02  valuev
 * fixed e-e energy for different spins, generalized 2-center term for vector operations
 *
 * Revision 1.47  2014/07/02 14:59:54  valuev
 * added forgotten yy term
 *
 * Revision 1.46  2014/07/02 09:36:16  valuev
 * added harmonic state projection onto many-electron (antisymmetrized) wave function
 *
 * Revision 1.45  2014/06/27 19:12:11  morozov
 * First version of projections to quantum harmonic oscillator eigenstates. Not fully working.
 *
 * Revision 1.44  2014/06/27 10:38:33  valuev
 * added second spin to state_io, checked norm_matrix, corrected y derivative for box
 *
 * Revision 1.43  2014/05/09 12:03:21  kazeev
 * Added directional potential multipliers into box.
 *
 * Revision 1.42  2014/04/21 18:42:41  kazeev
 * Added box_eterm_deriv. Untested.
 *
 * Revision 1.41  2013/10/17 12:32:01  morozov
 * Added AWPMD_split::el_density for calculation of |Psi(x)|^2
 *
 * Revision 1.40  2013/03/04 16:46:17  valuev
 * added rebuild_wp interface
 *
 * Revision 1.39  2012/10/04 11:42:05  valuev
 * first working version of constrainde dynamics
 *
 * Revision 1.38  2012/09/28 11:25:40  valuev
 * updated degeneracy constraints
 *
 * Revision 1.37  2012/09/20 05:38:31  valuev
 * optimizer (test version)
 *
 * Revision 1.36  2012/08/30 09:52:07  morozov
 * Interaction code is moved to a new module interact.cpp. awpmds.E_der is used instead of forces for minimization.
 *
 * Revision 1.35  2012/08/22 11:57:51  valuev
 * changes for KIM compliance
 *
 * Revision 1.34  2012/06/29 10:50:13  valuev
 * added linear constraints
 *
 * Revision 1.33  2012/06/15 16:59:00  morozov
 * OpenMP for electron-ion interaction and 2D psi dump. Added switching between dE and dE/dt criteria.
 *
 * Revision 1.32  2012/04/29 08:52:29  valuev
 * added power-conditional constraint
 *
 * Revision 1.31  2012/04/17 20:42:31  valuev
 * fixed Psi calculation, started adding e-i Debye screening
 *
 * Revision 1.30  2012/04/03 10:26:41  valuev
 * added variable constraints
 *
 * Revision 1.29  2012/03/23 15:50:57  valuev
 * added Psi calculation
 *
 * Revision 1.28  2012/01/27 10:54:03  valuev
 * fixed wrong term in the Hartree norm matrix
 *
 * Revision 1.27  2012/01/26 10:03:19  valuev
 * added free norm option
 *
 * Revision 1.26  2012/01/23 12:36:17  valuev
 * added wave function overlap calculation
 *
 * Revision 1.25  2012/01/17 06:00:07  valuev
 * neb workflow
 *
 * Revision 1.24  2011/10/12 18:58:48  valuev
 * added polar mode for C
 *
 * Revision 1.23  2011/10/12 16:39:02  valuev
 * restructured  interaction function
 *
 * Revision 1.22  2011/10/08 11:21:14  valuev
 * added external force storage and external force power calculation
 *
 * Revision 1.21  2011/07/27 19:40:21  morozov
 * Added determinant of norm-matrix for debug
 *
 * Revision 1.20  2011/07/22 14:24:55  valuev
 * forces with norm matrix
 *
 * Revision 1.19  2011/06/22 08:31:38  valuev
 * derivative fixes
 *
 * Revision 1.2  2011/06/11 16:53:55  valuev
 * sync with LAMMPS
 *
 * Revision 1.1  2011/06/10 17:15:07  morozov
 * First Windows project with the correct directory structure
 *
 * Revision 1.17  2011/06/09 22:55:08  valuev
 * norm matrices
 *
 * Revision 1.16  2011/06/07 19:58:42  valuev
 * corrected partitions
 *
 * Revision 1.15  2011/06/07 17:43:00  valuev
 * added Y derivatives
 *
 * Revision 1.14  2011/06/03 08:13:33  valuev
 * added partitions to account for ghost atoms
 *
 * Revision 1.13  2011/06/01 23:45:35  valuev
 * modified for LAMMPS compatibility
 *
 * Revision 1.12  2011/05/28 17:16:22  valuev
 * fixed template<>, some fixes to UHF
 *
 * Revision 1.11  2011/05/27 08:43:52  valuev
 * fixed split packet antisymmetrized version
 *
 * Revision 1.10  2011/05/25 05:23:43  valuev
 * fixed variable transformation for norm matrix
 *
 * Revision 1.9  2011/05/24 19:54:32  valuev
 * fixed sqmatrix::iterator
 *
 * Revision 1.8  2011/05/20 21:39:49  valuev
 * separated norm calculation
 *
 * Revision 1.7  2011/05/14 18:56:19  valuev
 * derivative for ee split interactions
 *
 * Revision 1.6  2011/05/05 08:56:02  valuev
 * working split wp version
 *
 * Revision 1.5  2011/05/04 16:48:52  valuev
 * fixed syntax
 *
 * Revision 1.4  2011/05/04 09:04:48  valuev
 * completed wp_split (except for ee forces)
 *
 * Revision 1.3  2011/04/22 09:54:24  valuev
 * working on split WP
 *
 * Revision 1.1  2011/04/20 08:43:09  valuev
 * started adding split packet WPMD
 *
 *******************************************************************************/

#include "wpmd.h"
#include <set>


///\en Term to be supplied to AWPMD_split::exp_value_2 to calculate 1-particle density matrix:
///    n(x1,x1')
class density_term1_t {
protected:
  Vector_3 x1, x2;
public:
  typedef cdouble result_t;

  ///\en Constructs the term with two coords
  density_term1_t(const Vector_3 &x1_, const Vector_3& x2_) :x1(x1_), x2(x2_) {
  }

  result_t operator()(const WavePacket &bra, const WavePacket &ket) const {
    return (conj(bra(x1))*ket(x2));
  }

  void add2(cdouble k, const WavePacket &bra, const WavePacket &ket) {}
};


///\en Term to be supplied to AWPMD_split::exp_value_2 to calculate 1-particle density matrix:
///    n(x1,x1')
class density_term2_t {
protected:
  Vector_3 x1, x2;
public:
  typedef cdouble result_t;

  ///\en Constructs the term with two coords
  density_term2_t(const Vector_3 &x1_, const Vector_3& x2_) :x1(x1_), x2(x2_) {
  }

  result_t operator()(const WavePacket &bra, const WavePacket &ket) const {
    cdouble res = (conj(bra(x1))*ket(x2));
    return res*res;
  }
  void add2(cdouble k, const WavePacket &bra, const WavePacket &ket) {}
};

class AWPMD_split: public AWPMD {
  friend class box_projection_term_sum;
  friend struct box_projection_term_ax;
protected:
  int s_add, spl_add;
  bool valid_norms; ///\en norms are valid for the given WP configuration
  bool split_wp; ///<\en (assigned automatically) if on, we really have more than one split per wp
public:

  std::pair<double, double> interaction_electron_kinetic(WavePacket const &packet, int spin, double *erforce, double *ervfroce) override;

  double interaction_ee_single(WavePacket const &packet_1,
                               WavePacket const &packet_2, double **eforce,
                               double **erforce, double cutoff) override;

  double interaction_ei_single(double const *x, double q, WavePacket const &packet,
                               int spin, double *f, double *eforce,
                               double *erforce, double cutoff) override;

  vector<Vector_2> split_c[2]; ///<\en split coefficients for electrons (c_re, c_im)  or (psi,phi) depending on the norm mode
  vector<int> nspl[2]; ///<\en number of wave packets for each electron (size is ne[i])

  vector<double> wf_norm[2]; ///<\en norms for each electron
  vector<double> wf_norm_der[2]; ///<\en norm derivative
  vector<cdouble> ovl_der[2]; ///<\en overlap derivative: \<psi|psi'\>


  vector< cdouble > Lh[2]; ///<\en Substitute for L in Hartree case: block matrices 1x(10*nspl[i])
  vector< sqmatrix<double> > Normh[2];  ///<\en Substitute for Norm in Hartree case: block matrices


  vector< vector<size_t> > constr_pivots[2]; ///\en Pivots for excluded/remaining variables as returned by LAPACK, first nconstr elements are indicies of excluded vars
  vector< vector<size_t> > rconstr_pivots[2]; ///\en Reverse pivots mapping initial var index to the place in constr_pivots array
  vector< recmatrix<double> > Jh[2];  ///<\en Constraint Jacobian for the Hartree case: Jij = dqi/dqj, where i is excluded variable, j is free variable

  ///\en If set to 1, uses amplitude-phase physical representation of
  ///    c: c=rho*exp(i*alpha), c[0]=rho, c[1]=alpha.  Otherwise
  ///    standard complex representation is used: c[0]=c_re, c[2]=c_im
  int c_polar_mode;

  ///\en Parameter controling how the electron norms are entering the
  ///    equations of motion.  NORMALIZE means that each single
  ///    electron wavefunction is normalized, i.e. divided by
  ///    sqrt(norm),\n FREE menans that no norm constraints are
  ///    applied,\n PROJECTED means projection constraint keeping the
  ///    forces orthogonal to the norm derivatives (NOT IMPLEMENTED)
  enum {NORMALIZE=0, FREE=1, PROJECTED=2} norm_mode;

  ///\en resizes all internal arrays according to new electrons
  ///(wavepackets) added
  virtual void resize(int flag);
public:
  double tmp;

  AWPMD_split() : s_add(0), spl_add(0), valid_norms(false), c_polar_mode(0), norm_mode(NORMALIZE),
                  split_wp(true) {
    vars_per_wp=10;
  }

  ///\en Copy operator.  Copies all structure information and
  ///    energies, does not copy matrices and temporary data.
  ///    Force/energy recalculation reqires the call to interaction()
  AWPMD_split &operator=(const AWPMD_split &other){
    if(this==&other)
      return *this;
    *((AWPMD *)this)=(AWPMD)other;
    s_add=other.s_add;
    spl_add=other.spl_add;
    valid_norms=other.valid_norms;
    c_polar_mode=other.c_polar_mode;
    norm_mode=other.norm_mode;

    for(int s=0;s<2;s++){
      split_c[s]=other.split_c[s];
      nspl[s]=other.nspl[s];
      wf_norm[s]=other.wf_norm[s];
      //wf_norm_der[s];
      //ovl_der[s];
      //E_der[s];
      //F_extra[s];
    }
    return
        *this;
  }

  ///\en Calculates complex overlap of the total wave function with
  ///    the other total wavefunction: conj(Psi)*other.Psi. Calculates
  ///    norms if necessary.
  cdouble operator*(AWPMD_split& other);

  ///\en Calculates total wave function at given space position \a x (has a meaning only for 1 electron).
  ///    In case of many electrons returns sqrt(el_density(x)).
  ///    Flags \a integ_flag are used to indicate along which axis the integration should be performed.
  ///    For example the call with integ_flag = iVector_3(0,1,1)  will return density at x[0] integrated over y and z.
  cdouble Psi(const Vector_3& x, const iVector_3 &integ_flag = iVector_3(0,0,0) );

  ///\en Calculates total electron density at given space position \a x.
  ///    Flags \a integ_flag are used to indicate along which axis the integration should be performed.
  ///    For example the call with integ_flag = iVector_3(0,1,1)  will return density at x[0] integrated over y and z.
  double el_density(const Vector_3& x, const iVector_3 &integ_flag = iVector_3(0,0,0));

  ///\en Prepares to setup a new system of particles using \ref
  ///    add_ion(), \ref add_electron() and \ref add_split().  There
  ///    is no need to call this function when using \ref
  ///    set_electrons() and \ref set_ions() to setup particles.
  virtual void reset(){
    for(int s=0;s<2;s++){
      split_c[s].clear();
      nspl[s].clear();
    }
    s_add=0;
    spl_add=0;
    AWPMD::reset();
  }



  ///\en Setup electrons: forms internal wave packet representations.
  ///    If PBCs are used the coords must be within the range [0, cell)
  ///    the \a splits array defines the number of wavepackets required for each electron
  ///    the data for splits should be placed in the corresponding data arrays
  ///    \a c array contains the splits mixing  coefficints
  ///    \a n is the number of electrons of a given spin component
  ///    Electron velocity v is multiplied by mass to obtain momentum.
  ///    Default mass (-1) means me.
  ///    Electronic charges q are -1 by default (when q=NULL), otherwise the charges are assigned for each split
  ///    \a pw_is_vel, if on, indicates that the \a pw array contains velocities (width_momentum[i]/mass) instead of momenta
  int set_electrons(int s, int nel, const Vector_3P x, const Vector_3P v, const double* w, const double* pw, const Vector_2 *c, const int *splits, double mass=-1, const double *q=NULL, bool pw_is_vel = false, double q0 = -1);

  ///\en Overrides AWPMD:set_electrons to correctly work for no-split case
  virtual int set_electrons(int spin, int n, const Vector_3P x, const Vector_3P v, const double* w, const double* pw, double mass=-1, const double *q=NULL, bool pw_is_vel = false, double q0=-1){
    if(n<=0)
      return 0;
    vector<int> splits(n,1);
    vector<Vector_2> c(n,Vector_2(1,0));
    /*vector<double> qv(n,q0);
    if(q){
      for(int i=0;i<n;i++)
        qv[i]=q[i];
    }*/
    return set_electrons(spin,n,x,v,w,pw,&c[0],&splits[0],mass,q, pw_is_vel,q0);
  }
  ///\en Initializes splits from two arrays
  int set_electrons_spl2(int s, int nel, const Vector_3P x, const Vector_3P v, const double* w, const double* pw, const double *c1, const double *c2, const int *splits, double mass = -1, const double *q = NULL, bool pw_is_vel = false, double q0 = -1) {
    int n = 0;
    for (int i = 0; i < nel; i++)
      n += splits[i];
    vector<Vector_2> c(n);
    for (int i = 0; i < n; i++) {
      c[i][0] = c1[i];
      c[i][1] = c2[i];
    }
    return set_electrons(s, nel, x, v, w, pw, &c[0], splits, mass, q, pw_is_vel, q0);
  }

  ///\en Returns the (modified after normalization) splits into two arrays
  int  get_splits_spl2(int s, double *c0, double *c1) const ;


  ///\en Starts adding new electron: continue with \ref add_split functions.
  int add_electron(int s){
    if(s < 0 || s > 1)
      return LOGERR(-1,fmt_iv("AWPMD_split.add_electron: invalid spin setting (%d)!",s),LINFO);
    calc_state&= (CALC_ELECTRONS_SET|CALC_IONS_SET);
    calc_state|= CALC_ELECTRONS_SET;
    s_add=s;
    spl_add=0;
    valid_norms=false;
    return ne[s_add];
  }

  ///\en Adds a new split to current electron.
  ///    May change the arguments according to the constraints set.
  ///    \return global id of the wavepacket (starting from 0 for each spin s)
  ///    Electron velocity v is multiplied by mass to obtain momentum.
  ///    Default mass (-1) means me.
  ///    The tags must be nonzero, >0 for the local particle, <0 for ghost particle.
  ///    Unique particle id  is abs(tag).
  ///    Default tag (0) means inserting the current particle id as local particle.
  int add_split(Vector_3 &x, Vector_3 &v, double &w, double &pw, Vector_2 &c, double mass=-1, double q=-1., int tag=0);


  ///\en Gets current electronic coordinates, and (optionally) number of wave packets for each electron.
  ///    Arguments \a c and \a splits may be NULL if their output is not needed.
  int get_electrons(int spin, Vector_3P x, Vector_3P v, double* w, double* pw, Vector_2 *c, int *splits=NULL, double mass=-1);


  ///\en Overrides \ref AWPMD::get_electrons() to correctly work for no-split case
  ///    Transforms the momenta to velocity v according to mass setting (-1 means me)
  virtual int get_electrons(int spin, Vector_3P x, Vector_3P v, double* w, double* pw, double mass=-1){
    return get_electrons(spin,x,v,w,pw,NULL,NULL,mass);
  }


  ///\en Prepares force arrays according to \a flag setting for interaction()
  virtual void clear_forces(int flagi,Vector_3P fi, Vector_3P fe_x,
                            Vector_3P fe_p, double *fe_w, double *fe_pw, Vector_2P fe_c);



  ///\en Same as \ref interaction(), but using Hartee factorization (no antisymmetrization)
  ///    !!! NOT IMPLEMENTED COMPLETELY !!!
  int interaction_hartree(int flag=0, Vector_3P fi=NULL, Vector_3P fe_x=NULL,
                          Vector_3P fe_p=NULL, double *fe_w=NULL, double *fe_pw=NULL, Vector_2P fe_c=NULL);

  ///\en Calculates interaction in the system of ni ions + electrons
  /// the electonic subsystem must be previously setup by set_electrons, ionic by set_ions
  /// 0x1   -- give back ion forces \n
  /// 0x2   -- add ion forces to the existing set \n
  /// 0x4   -- calculate electronic forces  \n
  /// 0x8   -- add electronic forces to the existing arrays \n
  /// 0x10  -- calculate internal electronic derivatives only: \n
  ///          will not update electronic force arrays, which may be NULL, \n
  ///          the forces may be obtained then using \ref get_el_forces() for all WPs \n
  ///          or separately for each WP using \ref get_wp_force()
  /// if PBCs are used the coords must be within the range [0, cell)
  /// 0x20  -- multiply forces by the inverted norm matrix
  virtual int interaction(int flag=0, Vector_3P fi=NULL, Vector_3P fe_x=NULL,
                          Vector_3P fe_p=NULL, double *fe_w=NULL, double *fe_pw=NULL,
                          Vector_2P fe_c=NULL);

  bool check_overlap(int s1, int c1, int c2) {
    if (approx == HARTREE && !split_wp)
      return true;
    return Oflg[s1](c1, c2);
  }

  ///\en Get electronic forcess in the arrays provided, using calculated internal representation
  ///    Valid flag settings are:\n
  ///    0x4   -- overwrite existing forces  \n
  ///    0x8   -- add electronic forces to the existing arrays \n
  void get_el_forces(int flag, Vector_3P fe_x,
                     Vector_3P fe_p, double *fe_w, double *fe_pw, Vector_2P fe_c);


  void get_wp_force(int s, int ispl, Vector_3P fe_x, Vector_3P fe_p, double *fe_w, double *fe_pw, Vector_2P fe_c, int normalized=1){
    WavePacket wk=wp[s][ispl];
    int indw1=8*ispl;
    int indn1=(nvar[s]/10)*8+2*ispl;
    double *ptr= normalized ? &(dq_dt[s][indw1]) : &(E_der[s][indw1]);
    double *ptrc= normalized ? &(dq_dt[s][indn1]) : &(E_der[s][indn1]);

    for(int j=0;j<3;j++){
      (*fe_x)[j]=ptr[j];
      (*fe_p)[j]=ptr[j+3];
    }
    *fe_w=ptr[6];
    *fe_pw=ptr[7];
    *fe_c=Vector_2(ptrc[0],ptrc[1]);

    //wk.int2phys_der< eq_minus_second >(E_der[s].begin()+indw1,(double *)fe_x,(double *)fe_p,fe_w,fe_pw,1./one_h);
    //*fe_c=-Vector_2(E_der[s][indn1],E_der[s][indn1+1]);
  }


  ///\en Calculates block norms and derivatives
  ///    \param flag is the same as for force calculation specification
  ///    \param normalize if set TRUE, normalizes coefficients so all electron norms are 1
  void calc_norms(int flag, bool normalize=false);

  ///\en Calcualtes the overlap between two electrons taking all split WPs into account.
  ///    Norms must be pre-calculated by \ref calc_norms().
  cdouble overlap(int ic1, int s1, int c1,int ic2, int s2, int c2);


  ///\en Calculates Norm matrix together with auxilliary matrices L and M
  ///    The result is saved in AWPMD_split::Norm[s] or AWPMD_split::Normh[s][el] for HARTREE case
  virtual void norm_matrix(int s);

  ///\en Performs LU-factorization of the Norm matrix
  ///    AWPMD::Norm[s] is replaced by the LU matrix
  ///    (or AWPMD_split::Normh[s][el] for HARTREE case)
  virtual int norm_factorize(int s);

  ///\en Inverts Norm matrix.
  ///    AWPMD::Norm[s] is replaced by the inverted matrix (or AWPMD_split::Normh[s][el] for HARTREE case)
  virtual int norm_invert(int s);

  ///\en Get the determinant of the norm-matrix for the particles with spin s
  virtual double norm_matrix_det(int s);

  ////\en Get the determinant logarithm of the norm-matrix for the particles with spin s
  virtual double norm_matrix_detl(int s) override;

  ///\en Multiplies forces (dHtot/dq_i) stored in f*_* by the inverse norm matrix.
  int calc_norm_forces();

  int calc_norm_forces_simplectic();

  ///\en Calculates the power of external force using forces multiplied by the inverse
  ///    norm matrix (AWPMD::Wext = \sum_i (dExt/dq_i)(dq_i/dt)).
  ///    Using the fact that dq_i/dt = inv(Norm)*(dHtot/d_qi) = -force_qi
  void calc_ext_power();


  ///\en Set constraint for specified parameter of given electron. Parameter id is
  ///    a number from 0 to 9: 0,1,2=x 3,4,5=p, 6=w, 7=pw, 8= cre (or |c|), 9=cim (or arg(c)).
  ///    Supported fix flags: 0=free (default), 1=fixed, 2=conditional.
  ///    \return 1 if OK or <0 if parameters are out of range.
  int set_var_constraint(int s, int electron_id, int wp_id, int param_id, int fix_flag=1);

  ///\en Get parameters of the given wavepacket.
  virtual int get_wavepacket(int spin, int wp_id, Vector_3 *x, Vector_3 *p, double *w, double *pw, Vector_2 *c=NULL) const {
    if(c)
      *c=split_c[spin][wp_id];
    return AWPMD::get_wavepacket(spin,wp_id,x,p,w,pw);
  }

  ///\en Set parameters of the given wavepacket.
  virtual int set_wavepacket(int spin, int wp_id, const Vector_3 *x, const Vector_3 *p, const double *w, const double *pw, const Vector_2 *c=NULL,double mass=-1){
    if(c)
      split_c[spin][wp_id]=*c;
    return AWPMD::set_wavepacket(spin,wp_id,x,p,w,pw,NULL,mass);
  }

  ///\en Get wave packet and its complex amplitude by spin and id
  int get_wavepacket(int spin, int wp_id, cdouble &c, WavePacket &wpack) const {
    wpack = wp[spin][wp_id];
    c=cdouble(split_c[spin][wp_id][0],split_c[spin][wp_id][1]);
    return 1;
  }

  ///\en Set wave packet and its complex amplitude by spin and id
  int set_wavepacket(int spin, int wp_id, const cdouble &c, const WavePacket &wpack){
    wp[spin][wp_id] = wpack;
    split_c[spin][wp_id][0] = c.real();
    split_c[spin][wp_id][1] = c.imag();
    return 1;
  }

  ///\en Makes timestep \a dt of electronic component: q-> q + (dq_dt)*dt for each variable.
  ///    If flag contains 0x10 uses internal variables, \a spin <0 means go through all spins,
  ///    the vaector \a dq_dt_ should contain generalized force (for all used spins), if NULL, they are taken from previous \ref interaction().
  virtual int step(double dt, int flag =0, int spin =-1,  const vector<double> *dq_dt_=NULL);

  int constraint_int_power(Vector_3P fe_x, Vector_3P fe_p, double *fe_w, double *fe_pw, Vector_2P fe_c, int release_all=0);

  enum DEG_ACTIONS { IMPOSE=0x1, RELEASE=0x2};
  int check_degeneracy_constraint(int action=IMPOSE, double prj_tolerance=-1., int max_constr_number=-1);

  ///\en Calculates overlap matrices O[s] and the inverses Y[s] (for antisymmetrized approximations)
  ///    \return exitcode <0 when matrix inversion fails, 1 on success
  int calc_overlaps();

protected:
  ///\en Calculates energy derivative contribution from the term
  ///    pref*Vp1p2*Op1p2/sqrt(norm_e1*norm_e2).  Here p1, p2 are wave
  ///    packet indicies composed of f.e. (s1,c1,k1), ic1 is the
  ///    precomputed array wp index dv_aj_conj, dv_bj_conj are V
  ///    derivatives with respect to bra (j) side conjugated WP
  ///    parameters conj(a) and conj(b) dv_ak, dv_bk are V derivatives
  ///    with respect to ket (k) side WP parameters a and b.  The
  ///    factor 1./sqrt(norm_e1*norm_e2) is differentiated BUT must be
  ///    inculded in the prefactor (not included automatically). 
  ////   When \a external=1, the contribution is attributed to the external
  ///    force and F_extra is updated (not added to internal
  ///    force).
  ///    Otherwise E_der is updated.
  void eterm_deriv(int ic1,int s1, int c1,int k1,int ic2,int s2, int c2,int j2,cdouble pref,
                   const OverlapDeriv &o,cdouble v,cdouble dv_aj_conj,
                   cdouble dv_ak,cVector_3 dv_bj_conj, cVector_3 dv_bk, int external=0);

  struct packet_index_info{
    int index;
    int spin;
  };
  void eterm_deriv(packet_index_info packet_1, packet_index_info packet_2, cdouble pref,
                   const OverlapDeriv &o,cdouble v,cdouble dv_aj_conj,
                   cdouble dv_ak,cVector_3 dv_bj_conj, cVector_3 dv_bk, int external=0);
  ///\en Calculates energy derivative contribution from interaction with
  ///the potential box. Parameters meaning is the same as in eterm_deriv.
  void box_eterm_deriv(
      const unsigned int ic1,
      const unsigned int s1,
      const unsigned int c1,
      const unsigned int j1,
      const unsigned int ic2,
      const unsigned int s2,
      const unsigned int c2,
      const unsigned int k2,
      const cdouble& pref);

  ///\en OpenMP version of eterm_deriv
  void eterm_deriv_omp(vector<double> &E_der1,int ic1,int s1,int c1,int k1,
                       vector<double> &E_der2,int ic2,int s2,int c2,int j2,cdouble pref,
                       const OverlapDeriv &o,cdouble v,cdouble dv_aj_conj,
                       cdouble dv_ak,cVector_3 dv_bj_conj, cVector_3 dv_bk);

  ///\en Adds the derivatives of Y for the term v*Y[s](c2,c1)
  void y_deriv(cdouble v,int s, int c2, int c1, int external =0 );


  ///\en Calculates L and M matrices and the Norm matrix depending on force flag and approximation
  void calc_LMN(int flag);

  ///\en Calculates L and M matrices and the Norm matrix for spin s depending on approximation
  void calc_LMN_s(int s);

  ///\en Transform the calculated norm matrix to physical variables for given spin s.
  int norm_matrix2phys(int s);

  int set_var_constraint_by_id(int s, int var_id, int fix_flag);

  void prepare_hartree_nm(int s, int c, int indw, int indn, sqmatrix<double> &matr, vector<double> &rhs, size_t &nconstr, bool apply_lin_constr=true);


public:
  ///\en Returns the number of imposed degeneracy constraints
  size_t get_constr_count() const ;

  int prepare_constraints(int s, int c);

  ///\en Gets current number of splits for spin variable \a spin (may be 0 or 1) and electron id.
  int get_split_number(int spin, int electron_id) const {
    if(spin<0 || spin >1)
      return 0;
    if(electron_id<0 || electron_id >ne[spin])
      return 0;
    return nspl[spin][electron_id];
  }

  ///\en Converts internal forces into physical coordinates
  void forces2phys();

  void forces2phsy(packet_index_info packet_index, WavePacket const &packet,
      double* eforce, double *erforce, double *ervforce);

  ///\en Evaluates expectation value of a two-center term over many-body system wavefefunction.
  ///    The term is defined by operator(): cdouble term(const WavePacket &bra, const WavePAcket &ket), where
  ///    the arguments \a bra and \a ket are supplied to the operator in direct form (not conjugated).
  ///    SHOULD BE CALLED AFTER \ref interaction() to ensure the norms and overlaps are calculated!!!
  ///    Performs for spin s, if s=-1 performs summation over spins, s=-2 cross spin (not implemented)
  template <class term2_t>
  typename term2_t::result_t exp_value_2(term2_t &term, int s=-1);


  template <class term4_t>
  void exp_value_4(term4_t &term, int s=-1, bool self_term = false);

  ///\en Calculates projection of the total many-body wave function onto
  ///    the given harmonic state. The harmonic well is defined by 
  ///    BoxHamiltonian, see \ref set_box(). Prameters \a num and \a axis define
  ///    the state number (valid values are from 0 to \a max_order of the BoxHamiltonian) and
  ///    the axis (0=x, 1=y, 2=z). If \a axis value is -1, then the sum of all states contributing
  ///    to the energy h_plank*omega*num from all axes are included. 
  double project_harm_state(int num, int axis = -1);
# if 1
  ///\en Calculates mean value of an eigenfunction of 3D harmonic well at temperatureaT: <Q|f*exp(-Hharm/T)|Q> 
  ///    over the current total many-body wavefunction Q. The result uses harmonic well partition 
  ///    function dZ(Q) = det( <wp_i | exp(-Hharm/T)| wp_j >) and is
  ///    obtained via projection of individual wavepackets onto harmonic eigenstates.
  ///    Class \a eigenfunc_t defines operator(int axis, double h_omega, int num) giving the expectation value of a function \a f,
  ///    commuting with the energy of harmonic well for a well component in coordinate \a axis with quantum \a h_omega,
  ///    over an eigenthete of number \a num of the 1D harmonic well.
  template<class eigenfunc_t>
  typename eigenfunc_t::result_t get_harm_mean_value(const eigenfunc_t &f, double T, double *dz=NULL){

    typedef typename eigenfunc_t::result_t result_t;
    typedef typename eigenfunc_t::cresult_t cresult_t;
    typedef typename eigenfunc_t::value_t value_t;
    typedef typename eigenfunc_t::cvalue_t cvalue_t;
    result_t result(0.);
    int max_order = box.get_max_proj_order();
    if(max_order<0) // in case the box or projections are undefined
      return result;

    chmatrix Z[2];// density matrices for UHF
    hmatrix<cresult_t> F[2];
    vector<double> zd[2]; // diagonal densities for HARTREE
    vector<result_t> fd[2];
    double z_det_log[2];

    Vector_3 h_omega = box.get_eigen_energies();

    for(int s =0;s<2; s++){
      int nes = ne[s];
      z_det_log[s] = 0.;
      if(nes == 0)
        continue;

      recmatrix<cdouble> proj[3];
      //vector<cdouble> left[2][3], right[2][3];
      for(int i=0;i<3;i++){ // filling axial projections
        proj[i].init(nes,max_order+1);
        proj[i].Set(cdouble(0.,0.));

        for(int c1=0,ic1=0; c1<nes; ic1+=nspl[s][c1], c1++){
          double norm_pref = norm_mode==NORMALIZE ? 1./sqrt(wf_norm[s][c1])  : 1.;

          for(int j1=0;j1<nspl[s][c1];j1++){
            double cj_re=split_c[s][ic1+j1][0];
            double cj_im=split_c[s][ic1+j1][1];
            cdouble cj=norm_pref*cdouble(cj_re,cj_im);
            WavePacket wj=wp[s][ic1+j1];
            for(int pj=0;pj<=max_order;pj++)
              proj[i](c1,pj) = cj*conj(box.eigen_proj(pj, wj ,i));
            //proj[i](c1,pj)  =  i==0 || pj==0? cdouble(100.*((double)rand()/RAND_MAX-0.5),100.*((double)rand()/RAND_MAX-0.5)) : 0.;
            //proj[i](c1,pj)  =   pj<=1 ? cdouble(100.*((double)rand()/RAND_MAX-0.5),100.*((double)rand()/RAND_MAX-0.5)) : 0.;
          }
          // checking projection norm
          double pnorm = 0.;
          for(int pj=0;pj<=max_order;pj++)
            pnorm+=norm(proj[i](c1,pj));
          pnorm = sqrt(pnorm);
          // normalizing
          for(int pj=0;pj<=max_order;pj++)
            proj[i](c1,pj)/=pnorm;

        }
      }
      // setting overlaps
      if(approx!=HARTREE){
        Z[s].init(nes);
        F[s].init(nes);
      }
      else{
        zd[s].resize(nes);
        fd[s].resize(nes);
      }

      int ik=0;
      for(int k=0;k<nes;k++){

        int il=0;
        for(int l=0/*k+1*/;l<nes;il+=nspl[s][l],l++){ // incrementing block2 wp address
          if(l<k)
            continue;
          if(approx==HARTREE && l>k)
            break;

          cdouble Zkl = 0.;
          cresult_t Fkl(0.,0.);
          // multiplexing: spanning all projections
          for(int ix =0; ix<=max_order; ix++){  // contributing nums from x
            double Ex = h_omega[0]*ix + 1./2.;
            cdouble vx = proj[0](k,ix)*conj(proj[0](l,ix))*exp(-Ex/T);
            for(int iy =0; iy<=max_order; iy++){ // contributing nums from y, z is uniquely defined now
              double Ey = h_omega[1]*iy + 1./2.;
              cdouble vy = proj[1](k,iy)*conj(proj[1](l,iy))*exp(-Ey/T);
              for(int iz=0;iz<=max_order;iz++){
                double Ez = h_omega[2]*iz + 1./2.;
                cdouble vz = proj[2](k,iz)*conj(proj[2](l,iz))*exp(-Ez/T);
                cdouble dzkl = vx*vy*vz;
                Zkl += dzkl;
                cvalue_t dfkl = f(h_omega, iVector_3(ix,iy,iz));
                dfkl *= dzkl;
                Fkl = Fkl+dfkl;
              }
            }
          }
          if(approx!=HARTREE){
            Z[s].set(k,l,Zkl);
            F[s].set(k,l,Fkl);
          }
          else{
            zd[s][k] = real(Zkl);
            fd[s][k] = real(Fkl);
          }
        }
        ik+=nspl[s][k]; // incrementing block1 wp address
      }
      //3. inverting the density matrix
      int info=0;
      if(approx!=HARTREE){
        /*FILE *f1=fopen(fmt("matrO_%d.d",s),"wt");
        fileout(f1,Y[s],"%15g");
        fclose(f1);8*/

        ZPPTRF("L",&nes,Z[s].arr,&info);
        // analyze return code here
        if(info<0)
          return LOGERR(info,fmt_iv("AWPMD.get_harm_mean_value: call to ZPTRF failed (exitcode %d)!",info),LINFO);


        for(int i=0; i<nes; i++)
          z_det_log[s] += log(real(Z[s](i,i)));

        ZPPTRI("L",&nes,Z[s].arr,&info);
        if(info<0)
          return LOGERR(info,fmt_iv("AWPMD.get_harm_mean_value: call to ZPTRI failed (exitcode %d)!",info),LINFO);
        /*f1=fopen(fmt("matrY_%d.d",s),"wt");
        fileout(f1,Y[s],"%15g");
        fclose(f1);*/
      }
      else{ // HARTREE
        for(int i=0; i<nes; i++){
          z_det_log[s] += log(zd[s][i]);
          zd[s][i]= 1./zd[s][i];  //inverting
        }
      }
      // calculating mean value
      for(int c1=0;c1<nes;c1++){
        for(int c2=0;c2<nes;c2++){
          if(c2<c1) // taken as M factor into account
            continue;
          cdouble yy;
          if(approx==HARTREE)
            yy=cdouble(zd[s][c1],0.);
          else
            yy=Z[s](c2,c1);

          // only diagonal terms for Hartree, and overlap nonzero
          if((approx!=HARTREE || c2==c1)){
            int M12=(c1==c2 ? 1: 2);
            cresult_t t;
            if(approx==HARTREE)
              t = fd[s][c1];
            else
              t = F[s](c1,c2);
            t*= M12*yy; // result_t must define operator*=
            result+= real(t); // result_t must define operator+=

          } // !HARTREE or spins or overlap
        }// c2
      }// c1
    } //s

    if(dz)
      *dz = z_det_log[0]+z_det_log[1];
    return result;
  }
# endif


  double density1(const Vector_3 &x, int spin =-1) {
    density_term1_t dens(x, x);
    cdouble res = 0.;
    if (ne[0] && spin<=0)
      res += exp_value_2(dens, 0) / ne[0];
    if(ne[1] && (spin<0 || spin ==1))
      res += exp_value_2(dens, 1) / ne[1];
    return real(res);
  }

  double density2(const Vector_3 &x1, const Vector_3 &x2, int spin = -1) {
    density_term1_t dens_11(x1, x1), dens_22(x2,x2);
    cdouble res = 0.;
    cdouble rho10 = 0., rho11 = 0., rho20 = 0., rho21 = 0.;
    if (ne[0] > 0 && spin <=0) {
      rho10 = exp_value_2(dens_11, 0);
      rho20 = exp_value_2(dens_22, 0);
    }
    if (ne[1] > 0 && (spin<0 || spin == 1)) {
      rho11 = exp_value_2(dens_11, 1);
      rho21 = exp_value_2(dens_22, 1);
    }

    if (ne[0]>1 && (spin<0 || spin == 0))
      res += rho10*rho20 / (ne[0] -1)/ ne[0] ;
    if (ne[1]>1 && (spin<0 || spin == 1))
      res += rho11*rho21 / (ne[1] - 1)/ ne[1] ;
    if (ne[0]>0 && ne[1]>0 && spin<0)
      res += (rho10*rho21 + rho11*rho20) / ne[0] / ne[1] /2;

    if (approx == HARTREE) {
      density_term2_t dens2_12(x1, x2);
      if (ne[0] > 1 && (spin<0 || spin == 0))
        res -= exp_value_2(dens2_12, 0) / (ne[0] - 1) / ne[0];
      if (ne[1] > 1 && (spin<0 || spin == 1))
        res -= exp_value_2(dens2_12, 1) / (ne[1] - 1) / ne[1];
    }
    else{
      density_term1_t dens_12(x1, x2), dens_21(x2, x1);
      if (ne[0]>1 && (spin<0 || spin == 0))
        res -= exp_value_2(dens_12, 0)*exp_value_2(dens_21, 0) / (ne[0] - 1)/ne[0];
      if (ne[1]>1 && (spin<0 || spin == 1))
        res -= exp_value_2(dens_12, 1)*exp_value_2(dens_21, 1) / (ne[1] - 1)/ne[1];

    }
    return real(res);
  }

};





template <class term2_t>
typename term2_t::result_t AWPMD_split::exp_value_2(term2_t &term, int s){
  typename term2_t::result_t res(0.);
  // 1. calculating overlaps if needed (internal check inside)
  calc_norms(0);
  //2. calculating overlap matrix if needed (internal check inside)
  int info=calc_overlaps();
  if (info < 0) {
    LOGERR(info, fmt_iv("AWPMD_split.exp_value_2: overlap matrix inversion failed!"), LINFO);
    return res;
  }

  int si = s >= 0 ? s : 0;
  int sf = s >= 0 ? s + 1 : 2;

  for(int s1=si; s1<sf; s1++){
    //  single particle contribution
    int ic1=0;

#if 0  //def _OPENMP
    int nthreads = omp_get_max_threads();
    vector<double> *E_der1 = new vector<double>[nthreads];
    int nvars1 = nvar[s1];
    for(int ith=0; ith<nthreads; ith++)
      E_der1[ith].assign(nvars1,0.);
#endif

    for(int c1=0;c1<ne[s1];c1++){

      int ic2=0;
      for(int c2=0;c2<ne[s1];ic2+=nspl[s1][c2],c2++){
        if( !Oflg[s1](c1,c2) )
          continue; // non-overlapping WPs

        if(c2<c1) // taken as M factor into account
          continue;

        double sq_norm12= norm_mode==NORMALIZE ? sqrt(wf_norm[s1][c1]*wf_norm[s1][c2]) : 1.;
        double pref_norm=1./sq_norm12; // for external energy

        cdouble yy;
        if(approx==HARTREE)
          yy=1.;
        else
          yy=Y[s1](c2,c1);

        // only diagonal terms for Hartree, and overlap nonzero
        if((approx!=HARTREE || c2==c1) && check_overlap(s1,c1,c2)){
          // WP blocks:
          for(int j1=0;j1<nspl[s1][c1];j1++){
            cdouble cj1(split_c[s1][ic1+j1][0],split_c[s1][ic1+j1][1]);
            WavePacket wj1=wp[s1][ic1+j1];

            //OverlapDeriv o12;
            //if(flag&(0x8|0x4)) //electron forces needed
            //o12.set1(wj1);

            for(int k2=(c1==c2 ? j1: 0); k2<nspl[s1][c2];k2++){
              int M12=(c1==c2 && j1==k2 ? 1: 2);
              double M12pe, M12pf;
              _mytie(M12pe,M12pf)=check_part1(s1,ic1+j1,s1,ic2+k2)*M12;
              cdouble ck2(split_c[s1][ic2+k2][0],split_c[s1][ic2+k2][1]);

              WavePacket wk2=wp[s1][ic2+k2];
              if(pbc)
                move_to_image(wj1,wk2);

              //WavePacket wjk12=conj(wj1)*wk2;
              //cdouble I012 = wjk12.integral();

              //if(norm(I012)<1e-22) // zero overlap !
              //continue;

              cdouble part_jk12=conj(cj1)*ck2;
              //cVector_3 djk12=wjk12.b/(2.*wjk12.a);

              // 2 center term contribution
              if(M12pf){
                typename term2_t::result_t t = term(wj1,wk2);
                cdouble pref = M12*part_jk12*yy*pref_norm;
                t*= pref; // result_t must define operator*=
                res+= t; // result_t must define operator+=

                term.add2(pref, wj1, wk2);
              }

            }// k2
          }// j1
          //if(flag&(0x8|0x4) && approx!=HARTREE){ //electron forces needed
          //  // adding Y derivative term for all variables (te and tei terms)
          //  y_deriv(Tc1c2,s1,c2,c1);
          //  y_deriv(Tc1c2e,s1,c2,c1,1); // external energy term
          //}
        } // !HARTREE or spins or overlap
        // END single particle contribution
      }// c2
      ic1+=nspl[s1][c1]; // incrementing block1 wp address
    }// c1

#if 0 // def _OPENMP
    for(int ith=0; ith<nthreads; ith++)
      for(int i=0; i<nvars1; i++)
        E_der[s1][i] += E_der1[ith][i];

    delete[] E_der1;
#endif
  } // s1
  return res;
}

# if 0
///\en Four-argument density term
template <class term4_t>
void AWPMD_split::exp_value_4(term4_t &term, int s, bool self_term){
  // 1. calculating overlaps if needed (internal check inside)
  calc_norms(0);
  //2. calculating overlap matrix if needed (internal check inside)
  int info = calc_overlaps();
  if (info < 0) {
    LOGERR(info, fmt("AWPMD_split.exp_value_4: overlap matrix inversion failed!"), LINFO);
    return;
  }

  // BEGIN  main energy loop
  int si = s >= 0 ? s : 0;
  int sf = s >= 0 ? s + 1 : 2;
  //typename term4_t::result_t res(0.);
  for (int s1 = si; s1<sf; s1++) {
    //  BEGIN single particle contribution
    int ic1 = 0;
    for (int c1 = 0; c1<ne[s1]; c1++) {
      int ic2 = 0;
      for (int c2 = 0; c2<ne[s1]; ic2 += nspl[s1][c2], c2++) {
        if (!check_overlap(s1,c1,c2))
          continue; // non-overlapping WPs

        if(c2<c1) // taken as M factor into account
          continue;

        double sq_norm12 = norm_mode == NORMALIZE ? sqrt(wf_norm[s1][c1] * wf_norm[s1][c2]) : 1.;
        //double pref_ee = -h2_me / 2;  //-h2_me/(2*sq_norm12); // ekin
        //double pref_ei = coul_pref / sq_norm12;
        double pref_norm = 1. / sq_norm12; // for external energy

        // pair by pair sum
        // second block
        // e-e interaction
        for (int s2 = s1; s2<2; s2++) {
          if (approx == HARTREE && c1 != c2) // only Vkmkm terms for Hartree
            continue;
          if (/*s1==s2 &&*/ c2<c1) // pair selection term
            continue;

          int ic3 = 0; // starting index of the wp for current electron
          for (int c3 = 0; c3<ne[s2]; ic3 += nspl[s2][c3], c3++) { // incrementing block2 wp address
            if (!self_term) {
              //if (approx == HARTREE && s1 == s2 && c3 == c1) // [Vkkmn contribution for same spin is 0], also  no Vkkkk terms for Hartree (=)
              //  continue;
              if (s1 == s2 && c3 <= c1) // and general Coulomb sum: i<j (<) //??
                continue;
            }
            else { // for self_term (ex dft extensions) we add self interaction as well
              if (s1 == s2 && c3 < c1) // and general Coulomb sum: i<=j (<) //??
                continue;
            }

            int ic4 = 0;
            for (int c4 = 0; c4<ne[s2]; ic4 += nspl[s2][c4], c4++) {
              if (approx == HARTREE && c4 != c3) // only Vkmkm terms for Hartree
                continue;
              //if(s1==s2 && c4==c2) // Vklmm contribution for same spin is 0 for antisymmetrized approximations, also Vmmmm is 0 for Hartree
              //  continue;
              if (!self_term) {
                if (s1 == s2 && c4 <= c2) // pair selection term: Vklmm
                  continue;
              }
              else { // for self term (ex dft extensions) we add self interaction as well
                if (s1 == s2 && c4 < c2) // pair selection term: Vklmm
                  continue;
              }

              if (/*s1==s2 &&*/ c2 == c1 && c4<c3) // pair selection term
                continue;

              double sq_norm34 = norm_mode == NORMALIZE ? sqrt(wf_norm[s2][c3] * wf_norm[s2][c4]) : 1.;
              double pref_ee = coul_pref / (sq_norm12*sq_norm34);


              cdouble yy;
              double K = 1.;
              if (approx == HARTREE) {
                yy = 1.;
              }
              else {
                if (s1 == s2) { // same spin antisymmetrized term
                  yy = Y[s1](c2, c1)*Y[s1](c4, c3) - Y[s1](c2, c3)*Y[s1](c4, c1);  // check the order of c
                }
                else { // different spin atisymmetrized term
                  //if(approx==UHF)
                  //K=2.;
                  yy = K*Y[s2](c4, c3)*Y[s1](c2, c1);
                }
              }

              for (int j1 = 0; j1<nspl[s1][c1]; j1++) {
                cdouble cj1(split_c[s1][ic1 + j1][0], split_c[s1][ic1 + j1][1]);
                WavePacket wj1 = wp[s1][ic1 + j1];

                for (int k2 = (approx == HARTREE ? j1 : 0); k2<nspl[s1][c2]; k2++) {
                  int M12 = (c1 == c2 && j1 == k2 ? 1 : 2);
                  cdouble ck2(split_c[s1][ic2 + k2][0], split_c[s1][ic2 + k2][1]);
                  WavePacket wk2 = wp[s1][ic2 + k2];
                  if (pbc)
                    move_to_image(wj1, wk2);

                  WavePacket wjk12 = conj(wj1)*wk2;
                  cdouble I012 = wjk12.integral();
                  cdouble part_jk12 = conj(cj1)*ck2;
                  cVector_3 djk12 = wjk12.b / (2.*wjk12.a);


                  for (int j3 = 0; j3<nspl[s2][c3]; j3++) {

                    cdouble cj3(split_c[s2][ic3 + j3][0], split_c[s2][ic3 + j3][1]);
                    const WavePacket& wj3 = wp[s2][ic3 + j3];

                    // 3-2
                    WavePacket wjk32;
                    cdouble I032 = 0., part_jk32;
                    cVector_3 djk32;
                    if (s1 == s2 || approx == UHF) {
                      WavePacket wk2 = wp[s2][ic2 + k2];
                      if (pbc)
                        move_to_image(wj3, wk2);
                      wjk32 = conj(wj3)*wk2;
                      I032 = wjk32.integral();
                      part_jk32 = conj(cj3)*ck2;
                      djk32 = wjk32.b / (2 * wjk32.a);
                    }


                    for (int k4 = (approx == HARTREE ? j3 : 0); k4<nspl[s2][c4]; k4++) {
                      int M34 = (c3 == c4 && j3 == k4 ? 1 : 2);
                      double M0;
                      if (approx == HARTREE)
                        M0 = M12*M34;
                      else {
                        M0 = (c1 == c2 && c3 == c4 ? 1 : 2); // will have exchange term for different pairs instead of M12*M34 factor
                      }
                      if (self_term) {
                        if (s1 == s2 && c1 == c3 && c2 == c4) // this is self energy
                          M0 *= 0.5;
                      }
                      double Me, Mf;
                      _mytie(Me, Mf) = check_part1(s1, ic1 + j1, ic2 + k2, s2, ic3 + j3, ic4 + k4)*M0;
                      if (!Mf)
                        continue;

                      cdouble ck4(split_c[s2][ic4 + k4][0], split_c[s2][ic4 + k4][1]);

                      // 3-4
                      WavePacket wk4 = wp[s2][ic4 + k4];
# if 1
                      if (pbc)
                        move_to_image(wj3, wk4);
                      WavePacket wjk34 = conj(wj3)*wk4;
                      cdouble I034 = wjk34.integral();
                      if (norm(I034)>ovl_tolerance && norm(I012)>ovl_tolerance) {
                        cdouble part_jk34 = conj(cj3)*ck4;
                        cdouble Kj1j3k2k4 = Me*yy*I012*I034*part_jk12*part_jk34;  // Vklmn
                        term.add4(Kj1j3k2k4, wjk12, wjk34);
                      }// ovl_tolerance
# endif
# if 1
                      // 1-4
                      if (approx != HARTREE && (approx == UHF && s1 == s2)) {
                        wk4 = wp[s2][ic4 + k4];
                        if (pbc)
                          move_to_image(wj1, wk4);
                        WavePacket wjk14 = conj(wj1)*wk4;
                        cdouble I014 = wjk14.integral();
                        if (norm(I032)>ovl_tolerance && norm(I014)>ovl_tolerance) {

                          cdouble part_jk14 = conj(cj1)*ck4;
                          cdouble Kj1j3k4k2 = Me*yy*I032*I014*part_jk32*part_jk14;  // Vklmn
                          term.add4(-Kj1j3k4k2, wjk32, wjk14);

                        }// ovl_tolerance
                      } // s1==s2
# endif
                    }// k4
                  }// j3
                } // k2
              } // j1

            }
          } // c4

        } // s2
      }// c2
      ic1 += nspl[s1][c1]; // incrementing block1 wp address
    }// c1
  } // s1
  // END main energy loop
}

# endif

///\en converts the contents of Vector_2 container form complex to polar representation
template <class inp_it>
void convert2polar(inp_it beg, inp_it end){
  for(;beg!=end;++beg){
    complex<double> cc((*beg)[0],(*beg)[1]);
    double rho=sqrt(norm(cc));
    double alpha=arg(cc);
    *beg=Vector_2(rho,alpha);
  }
}


///\en Auxiliary class for summing complex values stored in std::vector<cdouble> of arbitrary size.
///    Definies initialization and operations +=, *=.
struct vec_sum_t{
  cdouble ini;
  vector<cdouble> v;
  vec_sum_t(double ini_re=0., double ini_im=0.):ini(ini_re,ini_im){}

  vec_sum_t &operator+=(const vec_sum_t &other){
    v.resize(other.v.size(),ini);
    for(size_t i=0;i<v.size();i++)
      v[i]+=other.v[i];
    return *this;
  }
  vec_sum_t &operator*=(const cdouble &c){
    for(size_t i=0;i<v.size();i++)
      v[i]*=c;
    return *this;
  }
};

inline vec_sum_t real(const vec_sum_t &vec){
  vec_sum_t res;
  res.v.resize(vec.v.size());
  for(size_t i=0;i<vec.v.size();i++)
    res.v[i]= cdouble(real(vec.v[i]),0.);
  return res;
}

inline vec_sum_t conj(const vec_sum_t &vec){
  vec_sum_t res;
  res.v.resize(vec.v.size());
  for(size_t i=0;i<vec.v.size();i++)
    res.v[i]= conj(vec.v[i]);
  return res;
}


///\en Term to ne supplied to AWPMD_split::exp_value_2 to calculate all projections contributing to state with
///    certain energy (given by sum of axial numbers). The sum is specified by \ref set_sate_sum()
class box_projection_term_sum{

protected:
  AWPMD_split *awpmds;
  int max_order;
  bool fill_all;
public:
  typedef vec_sum_t result_t;

  ///\en Constructs the term: should be called with \a awpmds with valid norms and overlaps
  box_projection_term_sum(AWPMD_split *awpmds_, int num = -1, bool fill_all_ = true):awpmds(awpmds_),max_order(num), fill_all(fill_all_){
    if(max_order<0) // take from box
      max_order =  awpmds->box.get_max_proj_order();

  }



  result_t operator()(const WavePacket &bra, const WavePacket &ket){
    result_t result;
    if(max_order<0) // in case the box or projections are undefined
      return result;
    vector<cdouble> left[3], right[3];
    for(int i=0;i<3;i++){ // filling axial projections
      left[i].resize(max_order+1);
      right[i].resize(max_order+1);
      for(int j=0;j<=max_order;j++){
        left[i][j]  = awpmds->box.eigen_proj(j, bra ,i);
        right[i][j] = awpmds->box.eigen_proj(j, ket ,i);
      }
    }

    int istart = max_order, iend = max_order;
    if(fill_all){ // filling all results to max_order to reuse projections
      istart = 0;
      result.v.resize(max_order+1,cdouble(0,0));
    }
    else // need only one element
      result.v.resize(1,cdouble(0,0));

    for(int num=istart;num<=iend;num++){
      // multiplexing: spanning all ix+iy+iz = num,  (num+1)*(num+2)/2 terms altogether
      for(int ix =0; ix<=num; ix++){  // contributing nums from x
        for(int iy =0; iy<=num-ix; iy++){ // contributing nums from y, z is uniquely defined now
          int iz = num-iy-ix;
          // the order reflects projection having wp at left: <wp|harm_n>
          result.v[num-istart] += (left[0][ix]*left[1][iy]*left[2][iz])*conj(right[0][ix]*right[1][iy]*right[2][iz]);
        }
      }
    }
    return result;
  }
  void add2(cdouble k, const WavePacket &bra, const WavePacket &ket) {}
};







# endif
