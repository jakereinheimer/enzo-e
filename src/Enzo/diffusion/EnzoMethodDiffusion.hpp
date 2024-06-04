#pragma once

#include <map>

#include "Cello/cello.hpp"
#include "Enzo/enzo.hpp"

//typedef CelloView<enzo_float,3> EFlt3DArray;

class EnzoMethodDiffusion : public Method{
    public:
        EnzoMethodDiffusion(ParameterGroup p);

        ~EnzoMethodDiffusion()
        {

        }

        /// Charm++ PUP::able migration constructor
        EnzoMethodDiffusion (CkMigrateMessage *m)
            : Method (m)
        { }

        PUPable_decl(EnzoMethodDiffusion);

        void pup (PUP::er &p) override;

        void ThermalFlux_iso(Block * block, std::map<int , EFlt3DArray>  cndflx);

        void viscosity(Block * block, std::map<int , EFlt3DArray>  *cndflx);

        void compute( Block * block) throw() override;

        std::string name () throw () override
        { return "diffusion"; }

        /// Compute maximum timestep for this method
        double timestep ( Block * block) throw() override;

        protected: // methods

        void compute_ (Block * block, EFlt3DArray U ) throw();

        void update_cons(EFlt3DArray &u, EFlt3DArray &out,
                  EFlt3DArray &xflux, EFlt3DArray &yflux,
                  EFlt3DArray &zflux, enzo_float dt, enzo_float dx,
                  enzo_float dy, enzo_float dz);

        protected: // attributes

        EFlt3DArray x1flux,x2flux,x3flux;

        enzo_float kappa_iso, kappa_aniso;
        std::map<int,EFlt3DArray> cndflx;
        bool visc_, cond_,iso_,aniso_;
        EFlt3DArray kappa;


};