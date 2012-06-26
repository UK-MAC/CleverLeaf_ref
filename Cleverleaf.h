#ifndef CLEVERLEAF_H 
#define CLEVERLEAF_H 

#include "LagrangianEulerianPatchStrategy.h"
#include "LagrangianEulerianIntegrator.h"

#include "SAMRAI/tbox/Pointer.h"
#include "SAMRAI/appu/VisItDataWriter.h"
#include "SAMRAI/mesh/StandardTagAndInitStrategy.h"
#include "SAMRAI/hier/PatchHierarchy.h"
#include "SAMRAI/hier/PatchLevel.h"
#include "SAMRAI/hier/IntVector.h"
#include "SAMRAI/pdat/CellVariable.h"
#include "SAMRAI/pdat/NodeVariable.h"
#include "SAMRAI/pdat/EdgeVariable.h"
#include "SAMRAI/geom/CartesianGridGeometry.h"

using namespace SAMRAI;

using namespace std;

#include <vector>

class Cleverleaf:
    public LagrangianEulerianPatchStrategy
{
    public:
        Cleverleaf(
                tbox::Pointer<hier::PatchHierarchy>,
                const tbox::Dimension& dim,
                tbox::Pointer<geom::CartesianGridGeometry> grid_geometry);

        void registerVisItDataWriter(tbox::Pointer<appu::VisItDataWriter>);
        void registerModelVariables(LagrangianEulerianIntegrator* integrator);

        void initializeDataOnPatch(
                hier::Patch&,
                double init_data_time,
                bool initial_time);

        void accelerate(
            hier::Patch& patch,
            double dt);

        void ideal_gas_knl(
                hier::Patch& patch,
                bool predict);

        void viscosity_knl(
                hier::Patch& patch);

        double calc_dt_knl(
                hier::Patch& patch);

        void pdv_knl(
            hier::Patch& patch,
            double dt,
            bool predict);

        void flux_calc_knl(
                hier::Patch& patch,
                double dt);
 
    private:

        tbox::Pointer<hier::PatchHierarchy> d_hierarchy;
        tbox::Pointer<appu::VisItDataWriter> d_visit_writer;

        tbox::Pointer<geom::CartesianGridGeometry> d_grid_geometry;

        const tbox::Dimension d_dim;

        hier::IntVector d_nghosts;

        /*
         * Variables
         */
        tbox::Pointer<pdat::NodeVariable<double> > d_velocity;
        tbox::Pointer<pdat::EdgeVariable<double> > d_massflux;
        tbox::Pointer<pdat::EdgeVariable<double> > d_volflux;
        tbox::Pointer<pdat::CellVariable<double> > d_pressure;
        tbox::Pointer<pdat::CellVariable<double> > d_viscosity;
        tbox::Pointer<pdat::CellVariable<double> > d_soundspeed;
        tbox::Pointer<pdat::CellVariable<double> > d_density;
        tbox::Pointer<pdat::CellVariable<double> > d_energy;
        tbox::Pointer<pdat::CellVariable<double> > d_volume;

        tbox::Pointer<pdat::CellVariable<double> > d_celldeltas;
        tbox::Pointer<pdat::CellVariable<double> > d_cellcoords;

        tbox::Pointer<pdat::NodeVariable<double> > d_vertexdeltas;
        tbox::Pointer<pdat::NodeVariable<double> > d_vertexcoords;

        /*
         * Variable contexts
         */
        tbox::Pointer<hier::VariableContext> d_plot_context;
};
#endif
