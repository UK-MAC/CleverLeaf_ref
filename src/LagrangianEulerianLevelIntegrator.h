#ifndef included_LagrangianEulerianLevelIntegrator
#define included_LagrangianEulerianLevelIntegrator

#include "LagrangianEulerianPatchStrategy.h"

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

#include <string>
#include <list>

using namespace SAMRAI;

/**
 * @class LagrangianEulerianLevelIntegrator
 *
 * Controls the steps necessary to advance the solution across the geometry.
 *
 * These methods are largely defined in SAMRAI::algs::TimeRefinementLevelStrategy, and
 * provide a mechanism for SAMRAI to use the LagrangianEulerianLevelIntegrator to
 * perform a complete timestep on the given solution state.
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
         * Constants for the various exchange points, used when registering variables.
         *
         * @{
         */
        /** Not exchanged */
        const static int NO_EXCH = 0;
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
         * class implementing the LagrangianEulerianPatchStrategy methods in such
         * a way as to perform the desired physics on a give patch.
         *
         * @param object_name Name for this object, currently unused.
         * @param input_db Input database containing necessary setup info.
         * @param patch_strategy Patch strategy object to use.
         */
        LagrangianEulerianLevelIntegrator(
                const std::string& object_name,
                const boost::shared_ptr<tbox::Database>& input_db,
                LagrangianEulerianPatchStrategy* patch_strategy);

        /**
         * Default empty destructor.
         */
        ~LagrangianEulerianLevelIntegrator();

        /**
         * Register a variable with the integrator, allowing
         * it to be correctly transferred at halo exchanges as well as
         * coarsen/refine times
         *
         * @param var The variable to register
         * @param var_type The type of variable
         * @param var_exchanges The exchanges this variable is involved in
         * @param ghosts The number of ghosts this variable has
         * @param transfer_geom
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

        /*
         * TimeRefinementLevelStrategy methods
         */
        
        /**
         * Initialize the levelIntegrator.
         *
         * Sets up the model variables.
         *
         * @param gridding_alg
         */
        void initializeLevelIntegrator(
                const boost::shared_ptr<mesh::GriddingAlgorithmStrategy>& gridding_alg);

        double getMaxFinerLevelDt(
                    const int finer_level_number,
                    const double coarse_dt,
                    const hier::IntVector& ratio);

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

        bool usingRefinedTimestepping() const;

        /*
         * StandardTagAndInitialize methods
         */
        void initializeLevelData (
                const boost::shared_ptr<hier::PatchHierarchy>& hierarchy,
                const int level_number,
                const double init_data_time,
                const bool can_be_refined,
                const bool initial_time,
                const boost::shared_ptr<hier::PatchLevel>& old_level =
                    boost::shared_ptr<hier::PatchLevel>(),
                const bool allocate_data=true);

        void resetHierarchyConfiguration (
                const boost::shared_ptr<hier::PatchHierarchy>& hierarchy,
                const int coarsest_level,
                const int finest_level);

        void applyGradientDetector (
                const boost::shared_ptr<hier::PatchHierarchy>& hierarchy,
                const int level_number,
                const double error_data_time,
                const int tag_index,
                const bool initial_time,
                const bool uses_richardson_extrapolation_too);

        void lagrangianPredictor(
                const boost::shared_ptr<hier::PatchLevel>& level,
                const double dt);

        void halfStepHaloExchange(
                const boost::shared_ptr<hier::PatchLevel>& level,
                const boost::shared_ptr<hier::PatchHierarchy>& hierarchy,
                const double current_time);

        void lagrangianCorrector(
                const boost::shared_ptr<hier::PatchLevel>& level,
                const double dt);

        void preCellHaloExchange(
                const boost::shared_ptr<hier::PatchLevel>& level,
                const boost::shared_ptr<hier::PatchHierarchy>& hierarchy,
                const double current_time);

        void advecCellSweep1(
                const boost::shared_ptr<hier::PatchLevel>& level);

        void preMomSweep1HaloExchange(
                const boost::shared_ptr<hier::PatchLevel>& level,
                const boost::shared_ptr<hier::PatchHierarchy>& hierarchy,
                const double current_time);

        void advecMomSweep1(
                const boost::shared_ptr<hier::PatchLevel>& level);

        void preMomSweep2HaloExchange(
                const boost::shared_ptr<hier::PatchLevel>& level,
                const boost::shared_ptr<hier::PatchHierarchy>& hierarchy,
                const double current_time);

        void advecCellSweep2(
                const boost::shared_ptr<hier::PatchLevel>& level);

        void advecMomSweep2(
                const boost::shared_ptr<hier::PatchLevel>& level);

        void timestepEoS(
                const boost::shared_ptr<hier::PatchLevel>& level);

        void preLagrangeHaloExchange(
                const boost::shared_ptr<hier::PatchLevel>& level,
                const boost::shared_ptr<hier::PatchHierarchy>& hierarchy,
                const double current_time);

        void primeBoundaryHaloExchange(
                const boost::shared_ptr<hier::PatchLevel>& level,
                const boost::shared_ptr<hier::PatchHierarchy>& hierarchy,
                const double current_time);

        void viscosity(
                const boost::shared_ptr<hier::PatchLevel>& level);

        void postViscosityHaloExchange(
                const boost::shared_ptr<hier::PatchLevel>& level,
                const boost::shared_ptr<hier::PatchHierarchy>& hierarchy,
                const double current_time);

        double calcDt(
                const boost::shared_ptr<hier::PatchLevel>& level);

        void stampDataTime(
                const boost::shared_ptr<hier::PatchLevel>& level,
                const double current_time);

        void debugLevel(
                const boost::shared_ptr<hier::PatchLevel>& level);
        
        void swapAdvecDir();

        /**
         * Copy new field variable values back to timelevel 0.
         *
         * @param level The level we are working on.
         */
        void resetField(
                const boost::shared_ptr<hier::PatchLevel>& level);

        /**
         * Copy variables in revert_vars back to timelevel 0.
         *
         * @param level The level we are working on.
         * @param hierarchy The patch hierarchy we are working on.
         * @param current_time The current simulation time.
         */
        void revert(
                const boost::shared_ptr<hier::PatchLevel>& level);

        void fillBoundaries();

        /*
         * Print out the field summary
         */
        void getFieldSummary(
                const boost::shared_ptr<hier::PatchLevel>& level,
                double* level_volume,
                double* level_mass,
                double* level_pressure,
                double* level_internal_energy,
                double* level_kinetic_energy);

    protected:
        /**
         * PatchStrategy contains user-specified methods needed for
         * operating on a patch in the AMR hierarchy.
         */
        LagrangianEulerianPatchStrategy* d_patch_strategy;

        /**
         * The name of this object.
         */
        std::string d_object_name;

        /**
         * The dimension of the problem.
         */
        const tbox::Dimension d_dim;

        /**
         * @name Variable contexts
         *
         * Variable contexts for the simulation
         * @{
         */
        /** "New" timelevel 1 variables. */
        boost::shared_ptr<hier::VariableContext> d_new;
        /** Scratch variables. */
        boost::shared_ptr<hier::VariableContext> d_scratch;
        /** Scratch space for "new" variables. */
        boost::shared_ptr<hier::VariableContext> d_scratch_new;
        /** "Old" timelevel 0 variables. */
        boost::shared_ptr<hier::VariableContext> d_current;
        /** Context used to write ViSiT dumps. */
        boost::shared_ptr<hier::VariableContext> d_plot_context;
        /**
         * @}
         */

        hier::ComponentSelector d_var_scratch_data;
        hier::ComponentSelector d_var_scratch_new_data;
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

        std::vector<boost::shared_ptr<xfer::RefineSchedule> > d_half_step_schedules;
        std::vector<boost::shared_ptr<xfer::RefineSchedule> > d_prime_halos_schedules;
        std::vector<boost::shared_ptr<xfer::RefineSchedule> > d_pre_lagrange_schedules;
        std::vector<boost::shared_ptr<xfer::RefineSchedule> > d_post_viscosity_schedules;
        std::vector<boost::shared_ptr<xfer::RefineSchedule> > d_pre_sweep1_cell_schedules;
        std::vector<boost::shared_ptr<xfer::RefineSchedule> > d_pre_sweep1_mom_schedules;
        std::vector<boost::shared_ptr<xfer::RefineSchedule> > d_pre_sweep2_mom_schedules;

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
