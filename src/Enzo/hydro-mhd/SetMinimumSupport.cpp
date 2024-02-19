// See LICENSE_ENZO file for license and copyright information

/// @file      enzo_SetMinimumSupport.cpp
/// @author    Greg Bryan
/// @date      November, 1998
/// @brief     Set the energy to provide minimal pressure support

#include "cello.hpp"

#include "enzo.hpp"
 
//----------------------------------------------------------------------
 
int EnzoBlock::SetMinimumSupport(enzo_float &MinimumSupportEnergyCoefficient,
                                 enzo_float minimum_pressure_support_parameter,
				 bool comoving_coordinates)
{
  const int in = cello::index_static();

  Field field = data()->field();
  if (field.num_permanent() > 0) {  // TODO: revisit if-clause. This could be
                                    // improved. (plus we probably want to
                                    // report an error when false)
 
    /* Compute cosmology factors. */
 
    enzo_float cosmo_a = 1.0, cosmo_dadt = 0.0;

    EnzoPhysicsCosmology * cosmology = enzo::cosmology();

    ASSERT ("EnzoBlock::SetMinimumSupport()",
	    "comoving_coordinates enabled but missing EnzoPhysicsCosmology",
	    ! (comoving_coordinates && (cosmology == NULL)) );

    if (comoving_coordinates) {
      cosmology ->compute_expansion_factor(&cosmo_a, &cosmo_dadt,time());
    }

    enzo_float CosmoFactor = 1.0/cosmo_a;
 
    /* Determine the size of the grids. */

    const int GridRank = cello::rank();
 
    int dim, size = 1, i;
    for (dim = 0; dim < GridRank; dim++)
      size *= GridDimension[dim];

    enzo_float * density         = (enzo_float*) field.values("density");
    enzo_float * total_energy    = (enzo_float *)field.values("total_energy");
    enzo_float * internal_energy = (enzo_float *)field.values("internal_energy");
    enzo_float * velocity_x      = (enzo_float*) field.values("velocity_x");
    enzo_float * velocity_y      = (enzo_float*) field.values("velocity_y");
    enzo_float * velocity_z      = (enzo_float*) field.values("velocity_z");

    /* Set minimum GE. */

    // ToDo: figure out how to properly configure this variable. It used to be
    // a static global variable of EnzoBlock, but it could never be configured
    const enzo_float GravitationalConstant = 1.0;

    const enzo_float gamma = enzo::fluid_props()->gamma();
    MinimumSupportEnergyCoefficient =
      GravitationalConstant/(4.0*cello::pi) / (cello::pi * (gamma*(gamma-1.0))) *
      CosmoFactor * minimum_pressure_support_parameter *
      CellWidth[0] * CellWidth[0];

    /* PPM: set GE. */
    
    const EnzoDualEnergyConfig& de_config
        = enzo::fluid_props()->dual_energy_config();
    bool uses_dual_energy_formalism;
    if (de_config.bryan95_formulation()){
      uses_dual_energy_formalism = true;
    } else if (de_config.is_disabled()){
      uses_dual_energy_formalism = false;
    } else { // de_config.modern_formulation() == true
      // it's unlikely that there would be issues, but raise error to be safe
      ERROR("EnzoBlock::SetMinimumSupport",
            "the method is untested with this formulation of the dual "
            "energy formalism");
      uses_dual_energy_formalism = true;
    }

    if (uses_dual_energy_formalism) {
      for (i = 0; i < size; i++)
	internal_energy[i] = MAX(internal_energy[i],
				 MinimumSupportEnergyCoefficient*density[i]);
      if (GridRank != 3) return ENZO_FAIL;
      for (i = 0; i < size; i++)
	total_energy[i] = 
	  MAX((enzo_float)
	      (  internal_energy[i] +
		 0.5*(velocity_x[i]*velocity_x[i] +
		      velocity_y[i]*velocity_y[i] +
		      velocity_z[i]*velocity_z[i])),
	      total_energy[i]);
								
    }
    else {
      fprintf(stderr, "not implemented.\n");
      return ENZO_FAIL;
    }
 
  } // end: if (field.num_permanent() > 0)
 
  return ENZO_SUCCESS;
}