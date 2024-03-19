// See LICENSE_CELLO file for license and copyright information
/// @file     EnzoInitialShearStream.hpp
/// @author   Matthew Abruzzo (matthewabruzzo@gmail.com)
/// @date     Fri May 19 2023
/// @brief    [\ref Enzo] Initialization routine for a shear stream

#ifndef ENZO_ENZO_INITIAL_SHEAR_STREAM_HPP
#define ENZO_ENZO_INITIAL_SHEAR_STREAM_HPP

class EnzoInitialShearStream : public Initial {
  /// @class    EnzoInitialShearStream
  /// @ingroup  Enzo
  /// @brief    [\ref Enzo] Initialize a shearing stream. Both gas phases have
  /// bulk flows (anti-)parallel to the x-axis

public: // interface

  /// Constructor
  EnzoInitialShearStream(int cycle, double time, double vshear, double chi, double lambda_pert, double rho_hot, double vel_pert, double radius, double smoothing_thickness);

  /// CHARM++ PUP::able declaration
  PUPable_decl(EnzoInitialShearStream);

  /// CHARM++ migration constructor
  EnzoInitialShearStream(CkMigrateMessage *m)
  : Initial (m), vshear_(0.), chi_(0.), lambda_pert_(0.), rho_hot_(0.), vel_pert_(0.), radius_(0.), smoothing_thickness_(0.)
  { }

  /// CHARM++ Pack / Unpack function
  void pup (PUP::er &p);
  /*
  /// CHARM++ Pack / Unpack function
  void pup (PUP::er &p)
  {
    // NOTE: update whenever attributes change

    TRACEPUP;
    Initial::pup(p);
  }
  */

public: // virtual methods

  /// Initialize the block
  virtual void enforce_block
  ( Block * block, const Hierarchy * hierarchy ) throw();

private: // attributes

  double vshear_;
  double chi_;
  double lambda_pert_;
  double rho_hot_;
  double vel_pert_;
  double radius_;
  double smoothing_thickness_;

};


#endif //ENZO_ENZO_INITIAL_SHEAR_STREAM_HPP
