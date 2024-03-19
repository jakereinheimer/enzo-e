// See LICENSE_CELLO file for license and copyright information
/// @file     EnzoInitialShearStream.cpp
/// @author   Matthew Abruzzo (matthewabruzzo@gmail.com)
/// @date     Fri May 19 2023
/// @brief    [\ref Enzo] Implementation of EnzoInitialShearStream

#include "enzo.hpp"
#include "charm_enzo.hpp"
#include "cello.hpp"

#include <cmath>      // std::tanh, std::exp
#include <random>     // random noise


//----------------------------------------------------------------------

EnzoInitialShearStream::EnzoInitialShearStream(int cycle, double time, double vshear, double chi,double lambda_pert, double rho_hot, double vel_pert, double radius, double smoothing_thickness)
  : Initial(cycle, time),vshear_(vshear),chi_(chi),lambda_pert_(lambda_pert),rho_hot_(rho_hot),vel_pert_(vel_pert),radius_(radius),smoothing_thickness_(smoothing_thickness)
{}

void EnzoInitialShearStream::pup (PUP::er &p)
{
  // NOTE: update whenever attributes change

  TRACEPUP;

  Initial::pup(p);

  p | vshear_;
  p | chi_;
  p | lambda_pert_;
  p | rho_hot_;
  p | vel_pert_;
  p | radius_;
  p | smoothing_thickness_;
}

//----------------------------------------------------------------------

static void handle_bfields_(Field& field, enzo_float target_bfield_x) {
  // this currently just handles the case where the bfield is uniform and it is
  //   - 0 along all axes (compatible with pure-hydro)
  //   - the x-axis is the only non-zero component
  bool has_bfield = (field.is_field("bfield_x") ||
                     field.is_field("bfield_y") ||
                     field.is_field("bfield_z"));
  bool has_interface_bfield = (field.is_field("bfieldi_x") ||
                               field.is_field("bfieldi_y") ||
                               field.is_field("bfieldi_z"));

  if ( !(has_bfield || has_interface_bfield) ){
    ASSERT("handle_bfields_",
           "Issue encountered while initializing shear-stream. User has "
           "specified a non-zero bfield, but hasn't initialized Enzo-E in a "
           "way that can handle bfields",
           target_bfield_x == 0);
    return;
  }

  // at this moment in time, no MHD solver can handle the case where has_bfield
  // is true & has_interface_bfield is false. Thus we assume are both true
  // right now (an error will be reported down below when we attempt to load a
  // non-existant field).
  // - In the future, we may want to revisit this if we add solvers (e.g.
  //   divergence cleaning) that don't need interface bfields

  using EFlt3DView = CelloView<enzo_float,3>;

  // define lambda function that assigns to all values in a CelloView
  auto assign_val = [](EFlt3DView arr, enzo_float val) {
    for (int iz = 0; iz < arr.shape(2); iz++){
      for (int iy = 0; iy < arr.shape(1); iy++){
        for (int ix = 0; ix < arr.shape(0); ix++){

          

          arr(iz,iy,ix) = val;
        }
      }
    }
  };

  // assign cell-centered values:
  assign_val(field.view<enzo_float>("bfield_x"), target_bfield_x);
  assign_val(field.view<enzo_float>("bfield_y"), 0.0);
  assign_val(field.view<enzo_float>("bfield_z"), 0.0);

  // assign face-centered values:
  assign_val(field.view<enzo_float>("bfieldi_x"), target_bfield_x);
  assign_val(field.view<enzo_float>("bfieldi_y"), 0.0);
  assign_val(field.view<enzo_float>("bfieldi_z"), 0.0);

  // now, let's update the total energy (recall: it's specific total energy)
  EFlt3DView density = field.view<enzo_float>("density");
  EFlt3DView total_energy = field.view<enzo_float>("total_energy");

  for (int iz = 0; iz < density.shape(2); iz++){
    for (int iy = 0; iy < density.shape(1); iy++){
      for (int ix = 0; ix < density.shape(0); ix++){
        total_energy(iz,iy,ix) +=
          0.5 * (target_bfield_x * target_bfield_x) / density(iz,iy,ix);
      }
    }
  }
}

//----------------------------------------------------------------------

void EnzoInitialShearStream::enforce_block
(Block * block, const Hierarchy * hierarchy) throw()
{

  const double gamma = enzo::fluid_props()->gamma();

  // parameters: to be read from the parameter file

  enzo_float vshear = vshear_; //The velocity difference inside vs outside the cylinder
  enzo_float chi = chi_; //density difference between these phases
  enzo_float lambda_pert = lambda_pert_; //how often the perturbations happen, 0.0 means just one pert and any higher number is the number of times there is one centered at the left edge
  enzo_float rho_hot = rho_hot_; //base density, the density of the outside
  enzo_float vel_pert = vel_pert_; //The max velocity of the perturbation
  enzo_float radius = radius_; //The radius of the cylinder
  enzo_float smoothing_thickness = smoothing_thickness_; //How dramatic the density contrast is at the boundary
  enzo_float hot_temp = pow(10,4); //The max temperature of the outside gas
  enzo_float cool_temp = 1.0; //The cool temp of the cylinder


  const double smoothing_thickness_lambda = 1.0; //How large the max pert value is, higher value means smaller
  const double pert_width = 1.0; //How wide the perturbation is, higher means wider
  const double eint_hot = 10.0/(rho_hot * (gamma - 1.0)); // specific thermal
  const double target_bfield_x = pow(10,-5); //this is 10^-5 teslas or 10^-3 microgauss
  const double max_random_noise = 0.0; //maximum value of random noise added to the perturbation


  Field field = block->data()->field();
  using EFlt3DView = CelloView<enzo_float,3>;
  EFlt3DView density = field.view<enzo_float>("density");
  EFlt3DView velocity_x = field.view<enzo_float>("velocity_x");
  EFlt3DView velocity_y = field.view<enzo_float>("velocity_y");
  EFlt3DView velocity_z = field.view<enzo_float>("velocity_z");
  EFlt3DView total_energy = field.view<enzo_float>("total_energy");
  EFlt3DView temperature = field.view<enzo_float>("temperature");
  EFlt3DView metal_density = field.view<enzo_float>("metal_density");

  EFlt3DArray internal_energy; // specific internal (i.e. thermal) energy field
  const bool dual_energy
    = enzo::fluid_props()->dual_energy_config().any_enabled();
  if (dual_energy) {
    internal_energy = field.view<enzo_float>("total_energy");
  }

  // fetch shape of grid (including ghost zones). There aren't typos
  const int mz = density.shape(0);
  const int my = density.shape(1);
  const int mx = density.shape(2);

  // calculate x,y,z values at all cell centers (including ghost zones)
  int gx,gy,gz;
  field.ghost_depth(field.field_id("density"),&gx,&gy,&gz);
  std::vector<double> x_vals = std::vector<double>(mx);
  std::vector<double> y_vals = std::vector<double>(my);
  std::vector<double> z_vals = std::vector<double>(mz);
  block->data()->field_cells (x_vals.data(), y_vals.data(), z_vals.data(),
                              gx, gy, gz);

  // Now, let's fill in the grid (this loop currently includes ghost zones)
  for (int iz = 0; iz<mz; iz++){
    for (int iy = 0; iy<my; iy++){
      for (int ix = 0; ix<mx; ix++){

        const double r_cyl = std::sqrt((y_vals[iy] * y_vals[iy]) +
                                       (z_vals[iz] * z_vals[iz]));

        const double local_tanh = std::tanh((r_cyl-radius)/smoothing_thickness);

        // compute local_overdensity, relative to the wind density.
        // -> currently, it's always less than chi (nominal density contrast)
        double local_overdensity = 0.5 * (chi + 1 + (chi - 1) * -local_tanh);

        //double local_temp = 0.5*(hot_temp-cool_temp) * local_tanh + 0.5*(hot_temp-cool_temp) + 1;

        double local_density = rho_hot * local_overdensity;

        // compute the velocity component along transverse axes (y & z axes)
        double local_vy = 0.0;
        double local_vz = 0.0;

        //perturbations:
        if (lambda_pert > 0.0) { // if stated over 0.0, it will create a sin wave that implements perturbations

          // In the original implementation, there was a comment saying that we
          // only want to apply the perturbations to the areas between the 2
          // densities (this is a faithful transcription of the logic).
          // -> NOTE: in current implementation local_overdensity is ALWAYS
          //          smaller than chi (the nominal density contrast).
          if ((local_overdensity > rho_hot) && (local_overdensity < chi)) { //This catches in the region between the two densitiy regions

            // the perturbation has a gaussian profile
            double profile = std::exp(-1 * pert_width * pow((r_cyl-radius),2)) / smoothing_thickness_lambda;

            // compute perturbation along longitudinal direction (along x-axis)
            double longitude_factor = std::sin(2*cello::pi*x_vals[ix]/lambda_pert);
              /*
            } else {
              ERROR("EnzoInitialCloud::enforce_block", "uncomment this later");
              //longitude_factor =
              //  std::exp(-1 * ((x_vals[ix]-pert_loc)*(x_vals[ix]-pert_loc)) /
              //           pert_width);
            }
            */

            // get the angle at the current location
            // - (this is like the azimuthal in cylindrical coordinates, except
            //    we have replaced z-axis with x-axis)
            // - Should we be using std::atan2(y_vals[iy], z_vals[iz])?
            double theta = std::atan2(y_vals[iy],z_vals[iz]);
            local_vz = vel_pert * profile * longitude_factor * std::cos(theta);
            local_vy = vel_pert * profile * longitude_factor * std::sin(theta);
            }

        } else { //if the user just specifies 0.0 or anything else, it just does the pertibation at 1 next to the left border (x=-2)
          if ((local_overdensity > rho_hot) && (local_overdensity < chi)) {

            // the perturbation has a gaussian profile
            double profile = std::exp(-1 * ((r_cyl-radius) * (r_cyl-radius)) /
                                      smoothing_thickness_lambda);

            // compute perturbation along longitudinal direction (along x-axis)
            double longitude_factor = std::exp(-1*pert_width*pow((x_vals[ix]+(x_vals[0]+1)),2))/smoothing_thickness_lambda;

            // get the angle at the current location
            // - (this is like the azimuthal in cylindrical coordinates, except
            //    we have replaced z-axis with x-axis)
            // - Should we be using std::atan2(y_vals[iy], z_vals[iz])?
            double theta = std::atan(y_vals[iy]/z_vals[iz]);
            local_vz = vel_pert * profile * longitude_factor * std::cos(theta);
            local_vy = vel_pert * profile * longitude_factor * std::sin(theta);
          }
        }
        
        //Random noise

        std::mt19937 rng(std::random_device{}());

        std::uniform_real_distribution<double> dist((-1)*max_random_noise, max_random_noise);

        local_vy+=dist(rng);
        local_vz+=dist(rng);

        // the following assumes uniform thermal pressure!
        double local_eint = eint_hot / local_density;

        enzo_float metal_fraction = 0.0125;

        density(iz,iy,ix) = static_cast<enzo_float>(local_density);
        metal_density(iz,iy,ix)=density(iz,iy,ix)*metal_fraction;
        //temperature(iz,iy,ix) = static_cast<enzo_float>(local_temp);
        velocity_x(iz,iy,ix) = static_cast<enzo_float>(vshear * -local_tanh);
        velocity_y(iz,iy,ix) = static_cast<enzo_float>(local_vy);
        velocity_z(iz,iy,ix) = static_cast<enzo_float>(local_vz);

        if (dual_energy) {
          internal_energy(iz,iy,ix) = static_cast<enzo_float>(local_eint);
        }

        total_energy(iz,iy,ix) =
          (static_cast<enzo_float>(local_eint) +
           0.5* ( (velocity_x(iz,iy,ix) * velocity_x(iz,iy,ix)) +
                  (velocity_y(iz,iy,ix) * velocity_y(iz,iy,ix)) +
                  (velocity_z(iz,iy,ix) * velocity_z(iz,iy,ix)) ) );
      }
    }
  }



  // handle the bfields (if there are any)
  handle_bfields_(field, static_cast<enzo_float>(target_bfield_x));

  // this is where we would initialize bfields
  block->initial_done();
}


