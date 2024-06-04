#include <algorithm> //std::min and std::max

#include "Cello/cello.hpp"
#include "Enzo/enzo.hpp"

#include "EnzoMethodDiffusion.hpp"

//typedef CelloView<enzo_float,3> EFlt3DArray;

void EnzoMethodDiffusion::viscosity(Block * block, std::map<int , EFlt3DArray>  *cndflx){

  //get field data
  Data * data = block->data();
  Field field = data->field();

  //use density field id to get dimensions of sim
  int mx,my,mz;
  int gx,gy,gz;

  field.dimensions  (field.field_id ("density"),&mx,&my,&mz);
  field.ghost_depth (field.field_id ("density"),&gx,&gy,&gz);

  //use dimensions to prepare bounds without ghost zones
  int is = mx; int ie = mx-gx;
  int js = my; int je = my-gy;
  int ks = mz; int ke = mz-gz;

  //retrieve the data we need for calc
  EFlt3DArray rho = field.view<enzo_float>("density");
  EFlt3DArray press = field.view<enzo_float>("press");
  EFlt3DArray kappa = field.view<enzo_float>("kappa");

  //get block location
  //int ix, iy, iz, nx, ny, nz;
  //block->index_global(&ix, &iy, &iz, &nx, &ny, &nz);

  //get block size
  double hx,hy,hz;
  block->cell_width(&hx,&hy,&hz);

  //get inverse square area for finite differences method
  //double dxi = 1.0/(hx*hx);
  //double dyi = 1.0/(hy*hy);
  //double dzi = 1.0/(hz*hz);

  // prepare other vars we need for calc
  enzo_float kappaf, dTdx, dTdy, dTdz;

  //more:

  return;

}
