/*
 * Copyright 2013 David Beckingsale.
 * 
 * This file is part of CleverLeaf.
 * 
 * CleverLeaf is free software: you can redistribute it and/or modify it under
 * the terms of the GNU General Public License as published by the Free
 * Software Foundation, either version 3 of the License, or (at your option)
 * any later version.
 * 
 * CleverLeaf is distributed in the hope that it will be useful, but WITHOUT
 * ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
 * FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for
 * more details.
 * 
 * You should have received a copy of the GNU General Public License along with
 * CleverLeaf. If not, see http://www.gnu.org/licenses/.
 */ 
#ifndef CLEVERLEAF_CLEVERLEAF_H_
#define CLEVERLEAF_CLEVERLEAF_H_

#include <vector>

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

#include "LagrangianEulerianPatchStrategy.h"
#include "LagrangianEulerianLevelIntegrator.h"

using namespace SAMRAI;
using namespace std;

/**
 * @class Cleverleaf
 *
 * Cleverleaf extends LagrangianEulerianPatchStrategy to provide necessary
 * physics.
 *
 * Cleverleaf implements the abstract methods in the
 * LagrangianEulerianPatchStrategy class with the concrete versions of the
 * methods needed to run the required physics on a patch.
 */
class Cleverleaf:
  public LagrangianEulerianPatchStrategy
{
  public:
    /**
     * Create a new CleverLeaf object
     *
     * @param input_database The InputDatabase containing setup parameters.
     * @param hierarchy The PatchHierarchy to use.
     * @param dim The dimension of the problem.
     * @param grid_geometry The GridGeometry to use.
     */
    Cleverleaf(
        boost::shared_ptr<tbox::Database> input_database,
        boost::shared_ptr<hier::PatchHierarchy> hierarchy,
        const tbox::Dimension& dim,
        boost::shared_ptr<geom::CartesianGridGeometry> grid_geometry);

    /**
     * Register a VisitDataWriter with the class.
     *
     * Variables are then added to the VisitDataWriter in order to be output
     * during visulisation dumps.
     *
     * @param writer The VisitDataWriter to register.
     */
    void registerVisItDataWriter(
        boost::shared_ptr<appu::VisItDataWriter> writer);

    void registerModelVariables(LagrangianEulerianLevelIntegrator* integrator);

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
        double* total_volume,
        double* total_mass,
        double* total_pressure,
        double* total_internal_energy,
        double* total_kinetic_energy);

    virtual void tagGradientDetectorCells(
        hier::Patch& patch,
        const double regrid_time,
        const bool initial_error,
        const int tag_index);

    void debug_knl(hier::Patch& patch);
  private:
    boost::shared_ptr<hier::PatchHierarchy> d_hierarchy;
    boost::shared_ptr<appu::VisItDataWriter> d_visit_writer;

    boost::shared_ptr<geom::CartesianGridGeometry> d_grid_geometry;

    const tbox::Dimension d_dim;

    hier::IntVector d_nghosts;

    boost::shared_ptr<tbox::Database> input_db;
    const std::string state_prefix;

    bool d_tag_all;
    double d_tag_q_threshold;
    double d_tag_density_gradient;
    double d_tag_energy_gradient;

    boost::shared_ptr<pdat::NodeVariable<double> > d_velocity;
    boost::shared_ptr<pdat::EdgeVariable<double> > d_massflux;
    boost::shared_ptr<pdat::EdgeVariable<double> > d_volflux;
    boost::shared_ptr<pdat::CellVariable<double> > d_pressure;
    boost::shared_ptr<pdat::CellVariable<double> > d_viscosity;
    boost::shared_ptr<pdat::CellVariable<double> > d_soundspeed;
    boost::shared_ptr<pdat::CellVariable<double> > d_density;
    boost::shared_ptr<pdat::CellVariable<double> > d_energy;
    boost::shared_ptr<pdat::CellVariable<double> > d_volume;

    boost::shared_ptr<pdat::CellVariable<double> > d_celldeltas;
    boost::shared_ptr<pdat::CellVariable<double> > d_cellcoords;

    boost::shared_ptr<pdat::NodeVariable<double> > d_vertexdeltas;
    boost::shared_ptr<pdat::NodeVariable<double> > d_vertexcoords;

    boost::shared_ptr<pdat::CellVariable<int> > d_level_indicator;

    boost::shared_ptr<hier::VariableContext> d_plot_context;

    int* d_exchange_fields;

    static const int FIELD_DENSITY0 = 0;
    static const int FIELD_DENSITY1 = 1;
    static const int FIELD_ENERGY0 = 2;
    static const int FIELD_ENERGY1 = 3;
    static const int FIELD_PRESSURE = 4;
    static const int FIELD_VISCOSITY = 5;
    static const int FIELD_SOUNDSPEED = 6;
    static const int FIELD_XVEL0 = 7;
    static const int FIELD_XVEL1 = 8;
    static const int FIELD_YVEL0 = 9;
    static const int FIELD_YVEL1 = 10;
    static const int FIELD_VOL_FLUX_X = 11;
    static const int FIELD_VOL_FLUX_Y = 12;
    static const int FIELD_MASS_FLUX_X = 13;
    static const int FIELD_MASS_FLUX_Y = 14;

    const static int g_rectangle = 1;
    const static int g_circle = 2;
    const static int g_point = 4;

    const static double g_small = 1.0e-16;
    const static double g_big = 1.0e+21;
};
#endif
