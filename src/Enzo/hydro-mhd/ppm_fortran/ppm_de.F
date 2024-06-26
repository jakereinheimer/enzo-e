c     See LICENSE_ENZO file for license and copyright information

#include "fortran.h"

c=======================================================================
c/////////////////////////  SUBROUTINE PPM_DE  \\\\\\\\\\\\\\\\\\\\\\\\\
c
      subroutine ppm_de(d, e, u, v, w, ge,
     &     gravity, gr_xacc, gr_yacc, gr_zacc,
     &     gamma, dt, nhy, dx, dy, dz,
     &     rank, in, jn, kn, start, pend,
     &     iflatten, ipresfree,
     &     iconsrec, iposrec,
     &     idiff, isteepen, idual, eta1, eta2,
     &     nsubgrids, lface, rface,
     &     fistart, fiend, fjstart, fjend,
     &     array, dindex, eindex,
     &     uindex, vindex, windex, geindex, tmp,
     &     ncolor, colorpt, coloff, colindex,
     &     pmin, dmin, error,
     &     ie_error_x,ie_error_y,ie_error_z,num_ie_error)
c     
c  PERFORMS PPM (DIRECT EULERIAN) UPDATE FOR ONE TIMESTEP (WRAPPER)
c
c  written by: Greg Bryan
c       (some of this is based on code originally written by Jim Stone)
c  date:       June, 1993
c  modified1:  November, 1994; by GB
c              Modified to work as an entirely callable routine within
c              an AMR package.  There is now only one exterior parameter:
c              MAX_ANY_SINGLE_DIRECTION.
c  modified2:  Summer(?), 2002; by Robert Harkness
c              Extended sweeps into ghost zones to update boundary conditions
c              within the directionally split scheme.
c  modified3:  July, 2003; by Alexei Kritsuk
c              Corrected permutation pattern for 2D runs. 
c
c  PURPOSE:  Advances the fluid equations by one timestep using the
c    PPM Direct Eulerian algorithm in 3-D.  This routine is a callable
c    wrapper for the 1d sweeps.  It assumes (in classic fortran style)
c    that the arguments are defined to dimension three; however, it does
c    the right thing if called with one or more collapsed dimensions
c    (i.e. with the width set to 1).
c
c  EXTERNALS:
c    x,y,zeuler_sweep - routines to compute the Eulerian step in
c                       one dimension
c
c  INPUTS:
c     d       - density field (includes boundary zones)
c     dx,y,z  - zone width arrays for each dimension
c     e       - total specific energy field
c     pend     - array (of dimension 3) specifying the end of the active
c               region for reach dimension (zero based)
c     eta1    - (dual) selection parameter for gas energy (typically ~0.1)
c     eta2    - (dual) selection parameter for total energy (typically ~0.001)
c     ge      - gas energy (used when idual = 1)
c     gr_x,y,zacc - gravitational acceleration fields
c     gravity - flag indicating whether or not to use gravity field (1 = yes)
c     i,j,kn  - dimensions of field arrays
c     idiff   - diffusion flag (0 = off)
c     idual   - dual energy formalism flag (0 = off)
c     ipresfree - pressure free flag (0 = off, 1 = on, i.e. p=0)
c     nhy     - cycle number (for better operator splitting)
c     rank    - dimension of problem (unused until modification3)
c     start   - array (of dimension 3) specifying the start of the active
c               region fo reach dimension (zero based)
c     tmp     - temporary work space (30 * largest_slice)
c     u       - x-velocity field
c     v       - y-velocity field
c     w       - z-velocity field
c
c  LOCALS:
c
c-----------------------------------------------------------------------
      implicit NONE
#define FORTRAN
#include "fortran_types.h"
c-----------------------------------------------------------------------
c
c  Arguments
c
      integer gravity, idiff, idual, iflatten, isteepen, nhy, rank,
     &        ipresfree, pend(3), in, jn, kn, nsubgrids, start(3),
     &        ncolor, coloff(ncolor)
      integer  iconsrec, iposrec
      
      integer fistart(nsubgrids*3), fiend(nsubgrids*3),
     &        fjstart(nsubgrids*3), fjend(nsubgrids*3), 
     &        lface(nsubgrids*3), rface(nsubgrids*3)
      integer dindex(nsubgrids*6), eindex(nsubgrids*6),
     &        uindex(nsubgrids*6), vindex(nsubgrids*6),
     &        windex(nsubgrids*6),geindex(nsubgrids*6),
     &     colindex(nsubgrids*6,ncolor)
      integer error
      integer ie_error_x(*),ie_error_y(*),ie_error_z(*)
      integer num_ie_error
      ENZO_REAL d(in,jn,kn), e(in,jn,kn), u(in,jn,kn),
     &     v(in,jn,kn), w(in,jn,kn),ge(in,jn,jn),
     &     gr_xacc(in,jn,kn), gr_yacc(in,jn,kn), gr_zacc(in,jn,kn),
     &        dx(in),dy(jn),dz(kn)
      ENZO_REAL dt, eta1, eta2, gamma, pmin, dmin
      ENZO_REAL array(1), colorpt(1)
c
c  Parameters
c
c
c  Locals
c
      integer i, ie, is, ixyz, j, je, js, k, ke, ks,
     &        n, nxz, nyz, nzz, ms
      integer i1, i2, j1, j2, k1, k2,ii
      integer ntmp

      ENZO_REAL tmp(1)
c
c\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\////////////////////////////////
c=======================================================================
c
c  Error check
c     
      if (max(in,jn,kn) .gt. MAX_ANY_SINGLE_DIRECTION) then
         write(6,*) 'PPM_DE: A grid dimension is too long.'
         write(6,*) '   (increase MAX_ANY_SINGLE_DIRECTION.)'
         write(6,*) 'in=',in,'jn=',jn,'kn=',kn,'max=',
     *        MAX_ANY_SINGLE_DIRECTION
         stop 
      endif
c
c  Convert arguments to usable form
c
      is = start(1) + 1
      js = start(2) + 1
      ks = start(3) + 1
      ie = pend(1) + 1
      je = pend(2) + 1
      ke = pend(3) + 1
      ms = max(in*jn, jn*kn, kn*in)
      nxz = ie-is+1
      nyz = je-js+1
      nzz = ke-ks+1

      i1 = 1
      i2 = in
      j1 = 1
      j2 = jn
      k1 = 1
      k2 = kn
c
c  Loop over directions, using a Strang-type splitting
c
      ixyz = mod(nhy,rank)

      ntmp = ms*(32+ncolor*4)
      do ii=1,ntmp
         tmp(ii)=0.0
      end do
      do n=ixyz,ixyz+rank-1
c
c  Update in x-direction
c
         if (mod(n,rank) .eq. 0 .and. nxz .gt. 1) then
c
c*$* ASSERT CONCURRENT CALL
c$DOACROSS LOCAL(k)
           do k=k1, k2
              do ii=1,ntmp
                 tmp(ii)=0.0
              end do
              call xeuler_sweep(k, d, e, u, v, w, ge, in, jn, kn,
     &             gravity, gr_xacc, idual, eta1, eta2,
     &             is, ie, js, je, ks, ke,
     &             gamma, pmin, dmin, dt, dx, dy, dz,
     &             idiff, iflatten, isteepen,
     &             iconsrec, iposrec,
     &             ipresfree,  nsubgrids, lface, rface,
     &             fistart, fiend, fjstart, fjend,
     &             dindex, eindex, geindex,
     &             uindex, vindex, windex, array,
     &             ncolor, colorpt, coloff, colindex,
     &             tmp(1+ms*0), tmp(1+ms*1), tmp(1+ms*2), tmp(1+ms*3),
     &             tmp(1+ms*4), tmp(1+ms*5), tmp(1+ms*6), tmp(1+ms*7),
     &             tmp(1+ms*8), tmp(1+ms*9), tmp(1+ms*10),tmp(1+ms*11),
     &             tmp(1+ms*12),tmp(1+ms*13),tmp(1+ms*14),tmp(1+ms*15),
     &             tmp(1+ms*16),tmp(1+ms*17),tmp(1+ms*18),tmp(1+ms*19),
     &             tmp(1+ms*20),tmp(1+ms*21),tmp(1+ms*22),tmp(1+ms*23),
     &             tmp(1+ms*24),tmp(1+ms*25),tmp(1+ms*26),tmp(1+ms*27),
     &             tmp(1+ms*28),tmp(1+ms*29),tmp(1+ms*30),
     &             tmp(1+ms*(31+0*ncolor)),tmp(1+ms*(31+1*ncolor)), 
     &             tmp(1+ms*(31+2*ncolor)),tmp(1+ms*(31+3*ncolor)),
     &             error,
     &             ie_error_x,ie_error_y,ie_error_z,num_ie_error
     &             )
           enddo
        endif
c     
c     Update in y-direction
c
        if (mod(n,rank) .eq. 1 .and. nyz .gt. 1) then
c           
c*$* ASSERT CONCURRENT CALL
c$DOACROSS LOCAL(i)
           do i=i1, i2
              do ii=1,ntmp
                 tmp(ii)=0.0
              end do
              call yeuler_sweep(i, d, e, u, v, w, ge, in, jn, kn,
     &             gravity, gr_yacc, idual, eta1, eta2,
     &             is, ie, js, je, ks, ke,
     &             gamma, pmin, dmin, dt, dx, dy, dz,
     &             idiff, iflatten, isteepen,
     &             iconsrec, iposrec,
     &             ipresfree,
     &             nsubgrids, lface, rface,
     &             fistart, fiend, fjstart, fjend,
     &             dindex, eindex, geindex,
     &             uindex, vindex, windex, array,
     &             ncolor, colorpt, coloff, colindex,
     &             tmp(1+ms*0), tmp(1+ms*1), tmp(1+ms*2), tmp(1+ms*3),
     &             tmp(1+ms*4), tmp(1+ms*5), tmp(1+ms*6), tmp(1+ms*7),
     &             tmp(1+ms*8), tmp(1+ms*9), tmp(1+ms*10),tmp(1+ms*11),
     &             tmp(1+ms*12),tmp(1+ms*13),tmp(1+ms*14),tmp(1+ms*15),
     &             tmp(1+ms*16),tmp(1+ms*17),tmp(1+ms*18),tmp(1+ms*19),
     &             tmp(1+ms*20),tmp(1+ms*21),tmp(1+ms*22),tmp(1+ms*23),
     &             tmp(1+ms*24),tmp(1+ms*25),tmp(1+ms*26),tmp(1+ms*27),
     &             tmp(1+ms*28),tmp(1+ms*29),tmp(1+ms*30),
     &             tmp(1+ms*(31+0*ncolor)),tmp(1+ms*(31+1*ncolor)), 
     &             tmp(1+ms*(31+2*ncolor)),tmp(1+ms*(31+3*ncolor)),
     &             error,
     &             ie_error_x,ie_error_y,ie_error_z,num_ie_error
     &             )
           enddo
c     
        endif
c
c  Update in z-direction
c
        if (mod(n,rank) .eq. 2 .and. nzz .gt. 1) then
c
c*$* ASSERT CONCURRENT CALL
c$DOACROSS LOCAL(j)
           do j=j1, j2
              do ii=1,ntmp
                 tmp(ii)=0.0
              end do
              call zeuler_sweep(j, d, e, u, v, w, ge, in, jn, kn,
     &             gravity, gr_zacc, idual, eta1, eta2,
     &             is, ie, js, je, ks, ke,
     &             gamma, pmin, dmin, dt, dx, dy, dz,
     &             idiff, iflatten, isteepen,
     &             iconsrec, iposrec,
     &             ipresfree,
     &             nsubgrids, lface, rface,
     &             fistart, fiend, fjstart, fjend,
     &             dindex, eindex, geindex,
     &             uindex, vindex, windex, array,
     &             ncolor, colorpt, coloff, colindex,
     &             tmp(1+ms*0), tmp(1+ms*1), tmp(1+ms*2), tmp(1+ms*3),
     &             tmp(1+ms*4), tmp(1+ms*5), tmp(1+ms*6), tmp(1+ms*7),
     &             tmp(1+ms*8), tmp(1+ms*9), tmp(1+ms*10),tmp(1+ms*11),
     &             tmp(1+ms*12),tmp(1+ms*13),tmp(1+ms*14),tmp(1+ms*15),
     &             tmp(1+ms*16),tmp(1+ms*17),tmp(1+ms*18),tmp(1+ms*19),
     &             tmp(1+ms*20),tmp(1+ms*21),tmp(1+ms*22),tmp(1+ms*23),
     &             tmp(1+ms*24),tmp(1+ms*25),tmp(1+ms*26),tmp(1+ms*27),
     &             tmp(1+ms*28),tmp(1+ms*29),tmp(1+ms*30),
     &             tmp(1+ms*(31+0*ncolor)),tmp(1+ms*(31+1*ncolor)), 
     &             tmp(1+ms*(31+2*ncolor)),tmp(1+ms*(31+3*ncolor)),
     &             error,
     &             ie_error_x,ie_error_y,ie_error_z,num_ie_error
     &             )
           enddo
        endif
c
      enddo
c
      return
      end

