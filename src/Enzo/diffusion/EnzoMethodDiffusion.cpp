// See LICENSE_CELLO file for license and copyright information

/// @file     enzo_EnzoMethodDiffusion.cpp
/// @author   Jake Reinheimer
/// @date     
/// @brief    Implements the EnzoMethodDiffusion class

#include <map>

#include "Enzo/assorted/assorted.hpp"

#include "Cello/cello.hpp"
#include "Enzo/enzo.hpp"

#include "EnzoMethodDiffusion.hpp"

#include "enzo_typedefs.hpp"

//typedef CelloView<enzo_float,3> EFlt3DArray;

//----------------------------------------------------------------------

EnzoMethodDiffusion::EnzoMethodDiffusion (ParameterGroup p)
  : Method(),
  visc_(p.value_logical("visc",false)),
  cond_(p.value_logical("cond",false)),
  iso_(p.value_logical("iso",false)),
  aniso_(p.value_logical("aniso",false))
{
  
    if (visc_ == true){
        cello::define_field ("nu");

        cello::simulation()->refresh_set_name(ir_post_,name());
        Refresh * refresh = cello::refresh(ir_post_);
        refresh->add_field("nu");

    } else if (cond_ == true) {
        cello::define_field("kappa");

        cello::simulation()->refresh_set_name(ir_post_,name());
        Refresh * refresh = cello::refresh(ir_post_);
        refresh->add_field("kappa");

        std::map<int , EFlt3DArray> cndflx;

    } else {
        cello::define_field ("nu");
        cello::define_field("kappa");

        cello::simulation()->refresh_set_name(ir_post_,name());
        Refresh * refresh = cello::refresh(ir_post_);
        refresh->add_field("nu");
        refresh->add_field("kappa");
    }
    
    
}

//----------------------------------------------------------------------

void EnzoMethodDiffusion::pup (PUP::er &p)
{

  // NOTE: change this function whenever attributes change

  TRACEPUP;

  Method::pup(p);

  p | visc_;
  p | cond_;
  p | iso_;
  p | aniso_;
}

//----------------------------------------------------------------------

void EnzoMethodDiffusion::compute ( Block * block) throw()
{

  if (block->is_leaf()) {

    Field field = block->data()->field();

    EFlt3DArray rho = field.view<enzo_float>("density");

    //compute_ (block,rho);
  }

  block->compute_done();
}

//----------------------------------------------------------------------

double EnzoMethodDiffusion::timestep ( Block * block ) throw()
{
  // initialize_(block);

  Data * data = block->data();
  Field field = data->field();

  const int id_temp = field.field_id ("temperature");

  int mx,my,mz;
  field.dimensions (id_temp,&mx,&my,&mz);
  const int rank = ((mz == 1) ? ((my == 1) ? 1 : 2) : 3);

  double hx,hy,hz;
  block->cell_width(&hx,&hy,&hz);

  double h_min = std::numeric_limits<double>::max();
  if (rank >= 1) h_min = std::min(h_min,hx);
  if (rank >= 2) h_min = std::min(h_min,hy);
  if (rank >= 3) h_min = std::min(h_min,hz);

  return 0.5*courant_*h_min*h_min/0.001;
}

//======================================================================

void EnzoMethodDiffusion::compute_ (Block * block,EFlt3DArray U) throw()
{

  Data * data = block->data();
  Field field = data->field();

  //use density field id to get dimensions of sim
  int mx,my,mz;
  int gx,gy,gz;

  field.dimensions  (field.field_id ("density"),&mx,&my,&mz);
  field.ghost_depth (field.field_id ("density"),&gx,&gy,&gz);

  double hx,hy,hz;
  block->cell_width(&hx,&hy,&hz);

  double dt = timestep(block);

  //create new array with new rho information for after flux calcs
  EFlt3DArray rho_new=EFlt3DArray(U);

  if (visc_){


  } else if (cond_) {

    //fill cndflx map with the right size arrays for flux
    cndflx[1]=EFlt3DArray(mz,my,mx+1);
    cndflx[2]=EFlt3DArray(mz,my+1,mx);
    cndflx[3]=EFlt3DArray(mz+1,my,mx);


    if (iso_){
    //call thermalflux_iso method, described in EnzoMethodConduction. Fills flux arrays of cndflx
    ThermalFlux_iso(block,cndflx);

    update_cons(U,rho_new,cndflx[1],cndflx[2],cndflx[3],dt,hx,hy,hz);

    U=rho_new;

    } else if (aniso_){
      
    }


  } else {



  }

    




}

//====================================================================
//update conservative variable

typedef double enzo_float;
typedef CelloView<enzo_float,3> EFlt3DArray;

void  EnzoMethodDiffusion::update_cons(EFlt3DArray &u, EFlt3DArray &out,
                  EFlt3DArray &xflux, EFlt3DArray &yflux,
                  EFlt3DArray &zflux, enzo_float dt, enzo_float dx,
                  enzo_float dy, enzo_float dz){
  enzo_float dtdx = dt/dx;
  enzo_float dtdy = dt/dy;
  enzo_float dtdz = dt/dz;

  for (int iz=1; iz<u.shape(0)-1; iz++) {
    for (int iy=1; iy<u.shape(1)-1; iy++) {
      for (int ix=1; ix<u.shape(2)-1; ix++) {
        out(iz,iy,ix) = (u(iz,iy,ix) -
                         dtdx*(xflux(iz,iy,ix) - xflux(iz,iy,ix-1)) -
                         dtdy*(yflux(iz,iy,ix) - yflux(iz,iy-1,ix)) -
                         dtdz*(zflux(iz,iy,ix) - zflux(iz-1,iy,ix)));
      }
    }
  }
}
