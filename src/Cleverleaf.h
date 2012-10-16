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
#include "SAMRAI/appu/CartesianBoundaryDefines.h"

using namespace SAMRAI;

using namespace std;

#include <vector>


/**
 * @class Cleverleaf
 *
 * Cleverleaf extends LagrangianEulerianPatchStrategy to provide necessary physics.
 *
 * Cleverleaf implements the abstract methods in the LagrangianEulerianPatchStrategy
 * class with the concrete versions of the methods needed to run the required physics
 * on a patch.
 */
class Cleverleaf:
    public LagrangianEulerianPatchStrategy
{
    public:

        /**
         * Constructor for Cleverleaf class.
         *
         * @param hierarchy Pointer to the PatchHierarchy we are using.
         * @param dim The dimension of the problem.
         * @param grid_geometry The geometry of the problem.
         */
        Cleverleaf(
                tbox::Pointer<hier::PatchHierarchy> hierarchy,
                const tbox::Dimension& dim,
                tbox::Pointer<geom::CartesianGridGeometry> grid_geometry);

        /**
         * Register a VisitDataWriter with the class.
         *
         * @param writer The VisitDataWriter to register.
         */
        void registerVisItDataWriter(tbox::Pointer<appu::VisItDataWriter> writer);

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

        void advec_cell(hier::Patch& patch,
                int sweep_number,
                ADVEC_DIR direction);

        void advec_mom(hier::Patch& patch,
                int sweep_number,
                ADVEC_DIR direction,
                ADVEC_DIR which_vel);
 
        void setPhysicalBoundaryConditions(
                    hier::Patch& patch,
                    const double fill_time,
                    const hier::IntVector& ghost_width_to_fill);

        void field_summary(
                hier::Patch& patch,
                double* vol,
                double* mass,
                double* press,
                double* ie,
                double* ke);

        /**
         * 
         * Reflect data at physical boundaries.
         *
         * This method applies the reflective boundary condition up to a give depth for the data passed in. Use this method for cell-centered quantities such as pressure. The reflected values are simply copies.
         *
         * @param data The data to array to reflect boundaries of.
         * @param boundary Which boundayr to reflect.
         * @param depth The depth of cells to reflect the boundary over.
         * @param ifirst The first non-ghost cell of each dimension.
         * @param ilast The last non-ghost cell of each dimension.
         * @param xmin The xmin index of the entire patch (including ghosts)
         * @param xmax The xmax index of the entire patch (including ghosts)
         * @param ymin The ymin index of the entire patch (including ghosts)
         * @param ymax The ymax index of the entire patch (including ghosts)
         * @param nx The number of cells in the x dimension
         */
        void reflectPhysicalBoundary(
                double* data,
                BdryLoc::Type boundary,
                int depth,
                hier::Index ifirst,
                hier::Index ilast,
                int xmin,
                int xmax,
                int ymin,
                int ymax,
                int nx);

        /**
         * 
         * Reflect node/edge data at physical boundaries.
         *
         * This method applies the reflective boundary condition up to a give depth for the data passed in. Use this method for node and edge centered quantities such as velocities or fluxes. The reflected values are negated in the approriate direction copies.
         *
         * @param data The data to array to reflect boundaries of.
         * @param boundary Which boundayr to reflect.
         * @param depth The depth of cells to reflect the boundary over.
         * @param ifirst The first non-ghost cell of each dimension.
         * @param ilast The last non-ghost cell of each dimension.
         * @param xmin The xmin index of the entire patch (including ghosts)
         * @param xmax The xmax index of the entire patch (including ghosts)
         * @param ymin The ymin index of the entire patch (including ghosts)
         * @param ymax The ymax index of the entire patch (including ghosts)
         * @param nx The number of cells in the x dimension
         */
        void reflectXNodeBoundary(
                double* data,
                BdryLoc::Type boundary,
                int depth,
                hier::Index ifirst,
                hier::Index ilast,
                int xmin,
                int xmax,
                int ymin,
                int ymax,
                int nx);

        /**
         * 
         * Reflect node/edge data at physical boundaries.
         *
         * This method applies the reflective boundary condition up to a give depth for the data passed in. Use this method for node and edge centered quantities such as velocities or fluxes. The reflected values are negated in the approriate direction copies.
         *
         * @param data The data to array to reflect boundaries of.
         * @param boundary Which boundayr to reflect.
         * @param depth The depth of cells to reflect the boundary over.
         * @param ifirst The first non-ghost cell of each dimension.
         * @param ilast The last non-ghost cell of each dimension.
         * @param xmin The xmin index of the entire patch (including ghosts)
         * @param xmax The xmax index of the entire patch (including ghosts)
         * @param ymin The ymin index of the entire patch (including ghosts)
         * @param ymax The ymax index of the entire patch (including ghosts)
         * @param nx The number of cells in the x dimension
         */
        void reflectYNodeBoundary(
                double* data,
                BdryLoc::Type boundary,
                int depth,
                hier::Index ifirst,
                hier::Index ilast,
                int xmin,
                int xmax,
                int ymin,
                int ymax,
                int nx);

        void reflectXEdgeBoundary(
                double* data,
                BdryLoc::Type boundary,
                int depth,
                hier::Index ifirst,
                hier::Index ilast,
                int xmin,
                int xmax,
                int ymin,
                int ymax,
                int nx);

        void reflectYEdgeBoundary(
                double* data,
                BdryLoc::Type boundary,
                int depth,
                hier::Index ifirst,
                hier::Index ilast,
                int xmin,
                int xmax,
                int ymin,
                int ymax,
                int nx);

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
