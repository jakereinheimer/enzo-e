#include <algorithm> //std::min and std::max
#include <map>

#include "Cello/cello.hpp"
#include "Enzo/enzo.hpp"

#include "EnzoMethodDiffusion.hpp"

static enzo_float limiter2(const enzo_float A, const enzo_float B);
static enzo_float limiter4(const enzo_float A, const enzo_float B, const enzo_float C, const enzo_float D);
static enzo_float vanleer (const enzo_float A, const enzo_float B);
static enzo_float minmod(const enzo_float A, const enzo_float B);
static enzo_float mc(const enzo_float A, const enzo_float B);

//typedef CelloView<enzo_float,3> EFlt3DArray;

void EnzoMethodDiffusion::ThermalFlux_iso(Block * block, std::map<int,EFlt3DArray> cndflx){

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


  //x direction
  for (int k=ks; k<=ke; ++k) {
    for (int j=js; j<=je; ++j) {
      for (int i=is; i<=ie+1; ++i) {

        kappaf = 0.5*(kappa(k,j,i)+kappa(k,j,i-1));

        dTdx = (press(k,j,i)/rho(k,j,i)
                -press(k,j,i-1)/rho(k,j,i-1))
                /hx;

        cndflx[1](k,j,i) = kappaf*dTdx;
      } //end x loop
    } //end y loop
  } //end z loop

  //y direction
   for (int k=ks; k<=ke; ++k) {
    for (int j=js; j<=je+1; ++j) {
      for (int i=is; i<=ie; ++i) {

          kappaf = 0.5*(kappa(k,j,i)+kappa(k,j-1,i));

          dTdy = (press(k,j,i)/rho(k,j,i)
                  -press(k,j-1,i)/rho(k,j-1,i))
                  /hy;

          cndflx[2](k,j,i) = kappaf*dTdy;
        } //end x loop
      } //end y loop
    } //end z loop

  //z direction
  for (int k=ks; k<=ke+1; ++k) {
    for (int j=js; j<=je; ++j) {
      for (int i=is; i<=ie; ++i) {

          kappaf = 0.5*(kappa(k,j,i)+kappa(k-1,j,i));

          dTdz = (press(k,j,i)/rho(k,j,i)
                  -press(k-1,j,i)/rho(k-1,j,i))
                  /hz;
                  
          cndflx[3](k,j,i) = kappaf*dTdz;
        } //end x loop
      } //end y loop
    } //end z loop

  return;

}


/*----------------------------------------------------------------------------*/
/* limiter2 and limiter4: call slope limiters to preserve monotonicity
 */

static enzo_float limiter2(const enzo_float A, const enzo_float B)
{
    /* slope limiter */
    return mc(A,B);
}

static enzo_float limiter4(const enzo_float A, const enzo_float B, const enzo_float C, const enzo_float D)
{
    return limiter2(limiter2(A,B),limiter2(C,D));
}

/*----------------------------------------------------------------------------*/
/* vanleer: van Leer slope limiter
 */

static enzo_float vanleer(const enzo_float A, const enzo_float B)
{
    if (A*B > 0) {
        return 2.0*A*B/(A+B);
    } else {
        return 0.0;
    }
}

/*----------------------------------------------------------------------------*/
/* minmod: minmod slope limiter
 */

static enzo_float minmod(const enzo_float A, const enzo_float B)
{
    if (A*B > 0) {
        if (A > 0) {
            return std::min(A,B);
        } else {
            return std::max(A,B);
        }
    } else {
        return 0.0;
    }
}

/*----------------------------------------------------------------------------*/
/* mc: monotonized central slope limiter
 */

static enzo_float mc(const enzo_float A, const enzo_float B)
{
    return minmod(2.0*minmod(A,B), (A + B)/2.0);
}