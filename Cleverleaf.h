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

        double computeStableDtOnPatch(
            hier::Patch& patch,
            const bool initial_time,
            const double dt_time);

        void accelerate(
            hier::Patch& patch,
            double dt);

        void ideal_gas_knl(
                int xmin,
                int xmax,
                int ymin,
                int ymax,
                double* density,
                double* energy,
                double* pressure,
                double* soundspeed);

        void viscosity_knl(
                int xmin,
                int xmax,
                int ymin,
                int ymax,
                double* density,
                double* pressure,
                double* viscosity,
                double* xvel0,
                double* yvel0,
                double* celldx,
                double* celldy);

        double calc_dt_knl(
            int xmin,
            int xmax,
            int ymin,
            int ymax,
            double* celldx,
            double* celldy,
            double* soundspeed,
            double* viscosity,
            double* pressure,
            double* xvel0,
            double* yvel0,
            double* density,
            double* energy,
            double* xarea,
            double* yarea,
            double* volume,
            double* cellx,
            double* celly);


        void pdv_knl(
            hier::Patch& patch,
            double dt,
            bool predict);

    private:

        tbox::Pointer<hier::PatchHierarchy> d_hierarchy;
        tbox::Pointer<appu::VisItDataWriter> d_visit_writer;

        tbox::Pointer<geom::CartesianGridGeometry> d_grid_geometry;

        const tbox::Dimension d_dim;

        hier::IntVector d_nghosts;

        /*
         * Variables
         */
        tbox::Pointer<pdat::CellVariable<double> > d_velocity;
        tbox::Pointer<pdat::CellVariable<double> > d_massflux;
        tbox::Pointer<pdat::CellVariable<double> > d_volflux;
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
