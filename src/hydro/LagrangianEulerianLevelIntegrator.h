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

#ifndef CLEVERLEAF_LAGRANGIANEULERIANLEVELINTEGRATOR_H_
#define CLEVERLEAF_LAGRANGIANEULERIANLEVELINTEGRATOR_H_

#include "LagrangianEulerianPatchStrategy.h"

#include <string>
#include <list>

#include "SAMRAI/algs/TimeRefinementLevelStrategy.h"
#include "SAMRAI/mesh/GriddingAlgorithmStrategy.h"
#include "SAMRAI/hier/PatchLevel.h"
#include "SAMRAI/hier/PatchHierarchy.h"
#include "SAMRAI/mesh/StandardTagAndInitStrategy.h"
#include "SAMRAI/tbox/Database.h"
#include "SAMRAI/xfer/RefineAlgorithm.h"
#include "SAMRAI/xfer/CoarsenAlgorithm.h"
#include "SAMRAI/hier/BaseGridGeometry.h"
#include "SAMRAI/tbox/Timer.h"
#include "SAMRAI/tbox/StartupShutdownManager.h"

using namespace SAMRAI;

/**
 * @class LagrangianEulerianLevelIntegrator
 *
 * Controls the steps necessary to advance the solution across one level of the
 * patch hierarchy.
 */
class LagrangianEulerianLevelIntegrator:
  public mesh::StandardTagAndInitStrategy
{
  public:
    /**
     * @name Types 
     * Constants representing the type of a given variable.
     * @{
     */
    /** Normal variable. */
    const static int NORMAL = 1;
    /** Field variable, copied back to tl0 at end of timestep. */
    const static int FIELD = 2;
    /** Revert variable, copied back to tl0 at half-step. */
    const static int REVERT = 4;
    /** Level indicator **/
    const static int INDICATOR = 8;
    /**
     * @}
     */

    /**
     * @name Exchanges 
     * Constants for the various exchange points, used when registering
     * variables.
     * @{
     */
    /** Not exchanged */
    const static int NO_EXCH = 256;
    /** Exchanged at startup */
    const static int PRIME_CELLS_EXCH = 1;
    /** Exchanged before the Lagrangian step */
    const static int PRE_LAGRANGE_EXCH = 2;
    /** Exchanged after viscosity kernel */
    const static int POST_VISCOSITY_EXCH = 4;
    /** Exchanged at half step */
    const static int HALF_STEP_EXCH = 8;
    /** Exchanged before cell advection sweep 1 */
    const static int PRE_SWEEP_1_CELL_EXCH = 16;
    /** Exchanged before momentum advection sweep 1 */
    const static int PRE_SWEEP_1_MOM_EXCH = 32;
    /** Exchanged before momentum advection sweep 2 */
    const static int PRE_SWEEP_2_MOM_EXCH = 64;
    /** Used to fill a new level */
    const static int FIELD_EXCH = 128;
    /**
     * @}
     */

    /**
     * Create a new LagrangianEulerianLevelIntegrator.
     *
     * The patch_strategy object is the key parameter, as here we provide the
     * class implementing the LagrangianEulerianPatchStrategy methods in such a
     * way as to perform the desired physics on a given patch.
     *
     * @param input_db Input database containing setup parameters.
     * @param patch_strategy Patch strategy object to use.
     */
    LagrangianEulerianLevelIntegrator(
        const boost::shared_ptr<tbox::Database>& input_db,
        LagrangianEulerianPatchStrategy* patch_strategy);

    /**
     * Default empty destructor.
     */
    ~LagrangianEulerianLevelIntegrator();

    /**
     * Register a variable with the integrator, allowing it to be correctly
     * transferred at halo exchanges as well as coarsen/refine times.
     *
     * @param var The variable to register.
     * @param var_type The type of variable.
     * @param var_exchanges The exchanges this variable is involved in.
     * @param ghosts The number of ghosts this variable has.
     * @param transfer_geom The type of geometry being used.
     */
    void registerVariable(
        const boost::shared_ptr<hier::Variable>& var,
        const int var_type,
        const int var_exchanges,
        hier::IntVector nghosts,
        const boost::shared_ptr<hier::BaseGridGeometry>& transfer_geom);

    /**
     * Get the context to be used for visualisation.
     *
     * @returns The data context used for visualisation.
     */
    boost::shared_ptr<hier::VariableContext> getPlotContext();

    /**
     * Initialize the levelIntegrator.
     *
     * Sets up the variables being used for the simulation.
     *
     * @param gridding_alg
     */
    void initializeLevelIntegrator(
        const boost::shared_ptr<mesh::GriddingAlgorithmStrategy>& gridding_alg);

    void standardLevelSynchronization(
        const boost::shared_ptr<hier::PatchHierarchy>& hierarchy,
        const int coarsest_level,
        const int finest_level,
        const double sync_time,
        const double old_time);

    void synchronizeNewLevels(
        const boost::shared_ptr<hier::PatchHierarchy>& hierarchy,
        const int coarsest_level,
        const int finest_level,
        const double sync_time,
        const bool initial_time);

    void initializeLevelData(
        const boost::shared_ptr<hier::PatchHierarchy>& hierarchy,
        const int level_number,
        const double init_data_time,
        const bool can_be_refined,
        const bool initial_time,
        const boost::shared_ptr<hier::PatchLevel>& old_level =
        boost::shared_ptr<hier::PatchLevel>(),
        const bool allocate_data=true);

    void resetHierarchyConfiguration(
        const boost::shared_ptr<hier::PatchHierarchy>& hierarchy,
        const int coarsest_level,
        const int finest_level);

    void applyGradientDetector(
        const boost::shared_ptr<hier::PatchHierarchy>& hierarchy,
        const int level_number,
        const double error_data_time,
        const int tag_index,
        const bool initial_time,
        const bool uses_richardson_extrapolation_too);

    /**
     * Run the Lagrangian predictor step on a given level.
     *
     * @param level The level to work on.
     * @param dt The dt to advance by.
     */
    void lagrangianPredictor(
        const boost::shared_ptr<hier::PatchLevel>& level,
        const double dt);

    /**
     * Perform the half-step halo exchange.
     *
     * @param level The level to work on.
     * @param hierarchy The current PatchHierarchy.
     * @param current_time The current simulation time.
     */
    void halfStepHaloExchange(
        const boost::shared_ptr<hier::PatchLevel>& level,
        const boost::shared_ptr<hier::PatchHierarchy>& hierarchy,
        const double current_time);

    /**
     * Run the Lagrangian corrector step on a given level.
     *
     * @param level The level to work on.
     * @param dt The dt to advance by.
     */
    void lagrangianCorrector(
        const boost::shared_ptr<hier::PatchLevel>& level,
        const double dt);

    /**
     * Perform the pre-cell advection halo exchange.
     *
     * @param level The level to work on.
     * @param hierarchy The current PatchHierarchy.
     * @param current_time The current simulation time.
     */
    void preCellHaloExchange(
        const boost::shared_ptr<hier::PatchLevel>& level,
        const boost::shared_ptr<hier::PatchHierarchy>& hierarchy,
        const double current_time);

    /**
     * Run the first cell advection sweep on a given level.
     *
     * @param level The level to work on.
     */
    void advecCellSweep1(const boost::shared_ptr<hier::PatchLevel>& level);

    /**
     * Perform the pre-momentum advection halo exchange.
     *
     * @param level The level to work on.
     * @param hierarchy The current PatchHierarchy.
     * @param current_time The current simulation time.
     */
    void preMomSweep1HaloExchange(
        const boost::shared_ptr<hier::PatchLevel>& level,
        const boost::shared_ptr<hier::PatchHierarchy>& hierarchy,
        const double current_time);

    /**
     * Run the first momentum advection sweep on a given level.
     *
     * @param level The level to work on.
     */
    void advecMomSweep1(const boost::shared_ptr<hier::PatchLevel>& level);

    /**
     * Perform the second pre-momentum halo exchange.
     *
     * @param level The level to work on.
     * @param hierarchy The current PatchHierarchy.
     * @param current_time The current simulation time.
     */
    void preMomSweep2HaloExchange(
        const boost::shared_ptr<hier::PatchLevel>& level,
        const boost::shared_ptr<hier::PatchHierarchy>& hierarchy,
        const double current_time);

    /**
     * Run the second cell advection sweep on a given level.
     *
     * @param level The level to work on.
     */
    void advecCellSweep2(const boost::shared_ptr<hier::PatchLevel>& level);

    /**
     * Run the second momentum advection sweep on a given level.
     *
     * @param level The level to work on.
     */
    void advecMomSweep2(const boost::shared_ptr<hier::PatchLevel>& level);

    /**
     * Run the timestep equation of state on a given level.
     *
     * @param level The level to work on.
     */
    void timestepEoS(const boost::shared_ptr<hier::PatchLevel>& level);

    /**
     * Perform the pre-Lagrange halo exchange.
     *
     * @param level The level to work on.
     * @param hierarchy The current PatchHierarchy.
     * @param current_time The current simulation time.
     */
    void preLagrangeHaloExchange(
        const boost::shared_ptr<hier::PatchLevel>& level,
        const boost::shared_ptr<hier::PatchHierarchy>& hierarchy,
        const double current_time);

    /**
     * Perform the prime boundary halo exchange.
     *
     * @param level The level to work on.
     * @param hierarchy The current PatchHierarchy.
     * @param current_time The current simulation time.
     */
    void primeBoundaryHaloExchange(
        const boost::shared_ptr<hier::PatchLevel>& level,
        const boost::shared_ptr<hier::PatchHierarchy>& hierarchy,
        const double current_time);

    /**
     * Run the viscosity step on a given level.
     *
     * @param level The level to work on.
     */
    void viscosity(const boost::shared_ptr<hier::PatchLevel>& level);

    /**
     * Perform the post-viscosity halo exchange.
     *
     * @param level The level to work on.
     * @param hierarchy The current PatchHierarchy.
     * @param current_time The current simulation time.
     */
    void postViscosityHaloExchange(
        const boost::shared_ptr<hier::PatchLevel>& level,
        const boost::shared_ptr<hier::PatchHierarchy>& hierarchy,
        const double current_time);

    /**
     * Calculate the dt for a given level.
     *
     * @param level The level to work on.
     *
     * @returns The safe dt for the level.
     */
    double calcDt(const boost::shared_ptr<hier::PatchLevel>& level);

    /**
     * Set the data time for a given level.
     *
     * @param level The level to work on.
     * @param current_time The time to stamp the data to.
     */
    void stampDataTime(
        const boost::shared_ptr<hier::PatchLevel>& level,
        const double current_time);

    /**
     * Run the debug kernel on a given level.
     *
     * @param level The level to work on.
     */
    void debugLevel(const boost::shared_ptr<hier::PatchLevel>& level);

    /**
     * Swap the direction of the advection sweeps.
     */
    void swapAdvecDir();

    /**
     * Copy new field variable values back to the current time level.
     *
     * @param level The level to work on..
     */
    void resetField(const boost::shared_ptr<hier::PatchLevel>& level);

    /**
     * Copy variables in revert_vars back to the current time level.
     *
     * @param level The level to work on.
     */
    void revert(const boost::shared_ptr<hier::PatchLevel>& level);

    /**
     * Calculate the field summary quantities for a given level.
     *
     * @param level The level to work on.
     * @param level_volume The total volume of the level.
     * @param level_mass The total mass of the level.
     * @param level_pressure The total pressure of the level.
     * @param level_internal_energy The total internal energy of the level.
     * @param level_kinetic_energy The total kinetic eneryg of the level.
     */
    void getFieldSummary(
        const boost::shared_ptr<hier::PatchLevel>& level,
        double* level_volume,
        double* level_mass,
        double* level_pressure,
        double* level_internal_energy,
        double* level_kinetic_energy,
        int* level_effective_cells);
  private:
    LagrangianEulerianPatchStrategy* d_patch_strategy;

    const tbox::Dimension d_dim;

    /**
     * @name Variable contexts
     *
     * Variable contexts for the simulation
     * @{
     */
    /** "New" timelevel 1 variables. */
    boost::shared_ptr<hier::VariableContext> d_new;
    /** "Old" timelevel 0 variables. */
    boost::shared_ptr<hier::VariableContext> d_current;
    /** Context used to write ViSiT dumps. */
    boost::shared_ptr<hier::VariableContext> d_plot_context;
    /**
     * @}
     */

    hier::ComponentSelector d_var_cur_data;
    hier::ComponentSelector d_var_new_data;

    std::list<boost::shared_ptr<hier::Variable> > d_field_vars;
    std::list<boost::shared_ptr<hier::Variable> > d_revert_vars;

    /**
     * @name Communication exchanges
     * 
     * Algorithms for various communication points in the application.
     *
     * @{
     */
    boost::shared_ptr<xfer::RefineAlgorithm> d_bdry_fill_half_step;
    boost::shared_ptr<xfer::RefineAlgorithm> d_bdry_fill_prime_halos;
    boost::shared_ptr<xfer::RefineAlgorithm> d_bdry_fill_pre_lagrange;
    boost::shared_ptr<xfer::RefineAlgorithm> d_bdry_fill_post_viscosity;
    boost::shared_ptr<xfer::RefineAlgorithm> d_bdry_fill_pre_sweep1_cell;
    boost::shared_ptr<xfer::RefineAlgorithm> d_bdry_fill_pre_sweep1_mom;
    boost::shared_ptr<xfer::RefineAlgorithm> d_bdry_fill_pre_sweep2_mom;

    boost::shared_ptr<xfer::RefineAlgorithm> d_fill_new_level;
    boost::shared_ptr<xfer::CoarsenAlgorithm> d_coarsen_field_data;

    boost::shared_ptr<xfer::CoarsenAlgorithm> d_coarsen_level_indicator;

    std::vector<boost::shared_ptr<xfer::RefineSchedule> >
      d_half_step_schedules;
    std::vector<boost::shared_ptr<xfer::RefineSchedule> >
      d_prime_halos_schedules;
    std::vector<boost::shared_ptr<xfer::RefineSchedule> >
      d_pre_lagrange_schedules;
    std::vector<boost::shared_ptr<xfer::RefineSchedule> >
      d_post_viscosity_schedules;
    std::vector<boost::shared_ptr<xfer::RefineSchedule> >
      d_pre_sweep1_cell_schedules;
    std::vector<boost::shared_ptr<xfer::RefineSchedule> >
      d_pre_sweep1_mom_schedules;
    std::vector<boost::shared_ptr<xfer::RefineSchedule> >
      d_pre_sweep2_mom_schedules;
    /**
     * @}
     */

    bool advect_x;

    boost::shared_ptr<hier::PatchHierarchy> d_current_hierarchy;

    static boost::shared_ptr<tbox::Timer> t_synchronize_levels_create;
    static boost::shared_ptr<tbox::Timer> t_synchronize_levels_fill;

    static boost::shared_ptr<tbox::Timer> t_fill_new_levels_create;
    static boost::shared_ptr<tbox::Timer> t_fill_new_levels_fill;

    static boost::shared_ptr<tbox::Timer> t_half_step_exchange_create;
    static boost::shared_ptr<tbox::Timer> t_half_step_exchange_fill;

    static boost::shared_ptr<tbox::Timer> t_prime_halos_exchange_create;
    static boost::shared_ptr<tbox::Timer> t_prime_halos_exchange_fill;

    static boost::shared_ptr<tbox::Timer> t_pre_lagrange_exchange_create;
    static boost::shared_ptr<tbox::Timer> t_pre_lagrange_exchange_fill;

    static boost::shared_ptr<tbox::Timer> t_post_viscosity_exchange_create;
    static boost::shared_ptr<tbox::Timer> t_post_viscosity_exchange_fill;

    static boost::shared_ptr<tbox::Timer> t_pre_sweep1_cell_exchange_create;
    static boost::shared_ptr<tbox::Timer> t_pre_sweep1_cell_exchange_fill;

    static boost::shared_ptr<tbox::Timer> t_pre_sweep1_mom_exchange_create;
    static boost::shared_ptr<tbox::Timer> t_pre_sweep1_mom_exchange_fill;

    static boost::shared_ptr<tbox::Timer> t_pre_sweep2_mom_exchange_create;
    static boost::shared_ptr<tbox::Timer> t_pre_sweep2_mom_exchange_fill;

    static boost::shared_ptr<tbox::Timer> t_tag_gradient_detector_cells;

    static boost::shared_ptr<tbox::Timer> t_kernel_pdv;
    static boost::shared_ptr<tbox::Timer> t_kernel_ideal_gas;
    static boost::shared_ptr<tbox::Timer> t_kernel_revert;
    static boost::shared_ptr<tbox::Timer> t_kernel_reset;
    static boost::shared_ptr<tbox::Timer> t_kernel_accelerate;
    static boost::shared_ptr<tbox::Timer> t_kernel_flux_calc;
    static boost::shared_ptr<tbox::Timer> t_kernel_advec_cell;
    static boost::shared_ptr<tbox::Timer> t_kernel_advec_mom;
    static boost::shared_ptr<tbox::Timer> t_kernel_viscosity;
    static boost::shared_ptr<tbox::Timer> t_kernel_calc_dt;
    static boost::shared_ptr<tbox::Timer> t_kernel_initialize_data;

    static void initializeCallback();
    static void finalizeCallback();

    static tbox::StartupShutdownManager::Handler s_initialize_handler;
};
#endif
