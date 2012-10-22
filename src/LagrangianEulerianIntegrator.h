#ifndef included_LagrangianEulerianIntegrator
#define included_LagrangianEulerianIntegrator

#include "LagrangianEulerianPatchStrategy.h"

#include "SAMRAI/algs/TimeRefinementLevelStrategy.h"
#include "SAMRAI/mesh/GriddingAlgorithmStrategy.h"
#include "SAMRAI/hier/PatchLevel.h"
#include "SAMRAI/hier/PatchHierarchy.h"
#include "SAMRAI/mesh/StandardTagAndInitStrategy.h"
#include "SAMRAI/tbox/Serializable.h"
#include "SAMRAI/tbox/Pointer.h"
#include "SAMRAI/tbox/Database.h"
#include "SAMRAI/tbox/List.h"
#include "SAMRAI/xfer/RefineAlgorithm.h"

#include <string>

using namespace SAMRAI;

/**
 * @class LagrangianEulerianIntegrator
 *
 * Controls the steps necessary to advance the solution across the geometry.
 *
 * These methods are largely defined in SAMRAI::algs::TimeRefinementLevelStrategy, and
 * provide a mechanism for SAMRAI to use the LagrangianEulerianIntegrator to
 * perform a complete timestep on the given solution state.
 */
class LagrangianEulerianIntegrator:
    public algs::TimeRefinementLevelStrategy,
    public mesh::StandardTagAndInitStrategy,
    public tbox::Serializable
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
        /**
         * @}
         */

        /**
         * Create a new LagrangianEulerianIntegrator.
         *
         * The patch_strategy object is the key parameter, as here we provide the
         * class implementing the LagrangianEulerianPatchStrategy methods in such
         * a way as to perform the desired physics on a give patch.
         *
         * @param object_name Name for this object, currently unused.
         * @param input_db Input database containing necessary setup info.
         * @param patch_strategy Patch strategy object to use.
         */
        LagrangianEulerianIntegrator(
                const std::string& object_name,
                tbox::Pointer<tbox::Database> input_db,
                LagrangianEulerianPatchStrategy* patch_strategy
                );

        /**
         * Default empty destructor.
         */
        ~LagrangianEulerianIntegrator();

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
                tbox::Pointer<hier::Variable> var,
                const int var_type,
                const int var_exchanges,
                hier::IntVector nghosts,
                const tbox::Pointer<hier::GridGeometry> transfer_geom);

        /**
         * Get the context to be used for visualisation.
         *
         * @returns The data context used for visualisation.
         */
        tbox::Pointer<hier::VariableContext> getPlotContext();

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
                tbox::Pointer<mesh::GriddingAlgorithmStrategy> gridding_alg);

        /**
         * Return the stable dt for the given level.
         *
         * @param level Level to compute dt for.
         * @param dt_time The current dt time.
         * @param initial_time True if we are at the initial time.
         *
         * @returns 
         */
        double getLevelDt(
                const tbox::Pointer<hier::PatchLevel> level,
                const double dt_time,
                const bool initial_time);

        /**
         * Advances the level from current_time to new_time.
         *
         * Takes the steps necessary to advance the level one timestep.
         * This routine utilizes the methods from the LagrangianEulerianPatchStrategy
         * to perform the necessary physics on each patch in the level.
         *
         * Communication and boundary computation is carried out at the appropriate
         * points. 
         *
         * @param level level to advance.
         * @param hierarchy the hierarchy of patches.
         * @param current_time the current simulation time.
         * @param new_time the time to advance to.
         * @param first_step true if this is the first step.
         * @param last_step true if this is the last step
         * @param regrid_advance
         *
         * @returns the next dt value.
         */
        double advanceLevel(
                const tbox::Pointer<hier::PatchLevel> level,
                const tbox::Pointer<hier::PatchHierarchy> hierarchy,
                const double current_time,
                const double new_time,
                const bool first_step,
                const bool last_step,
                const bool regrid_advance=false);

        void standardLevelSynchronization(
                const tbox::Pointer<hier::PatchHierarchy> hierarchy,
                const int coarsest_level,
                const int finest_level,
                const double sync_time,
                const double old_time);

        void synchronizeNewLevels(
                const tbox::Pointer<hier::PatchHierarchy> hierarchy,
                const int coarsest_level,
                const int finest_level,
                const double sync_time,
                const bool initial_time);

        void resetTimeDependentData(
                const tbox::Pointer<hier::PatchLevel> level,
                const double new_time,
                const bool can_be_refined);

        void resetDataToPreadvanceState(
                const tbox::Pointer<hier::PatchLevel> level);

        bool usingRefinedTimestepping() const;

        /*
         * StandardTagAndInitialize methods
         */
        void initializeLevelData (
                const tbox::Pointer<hier::PatchHierarchy> hierarchy,
                const int level_number,
                const double init_data_time,
                const bool can_be_refined,
                const bool initial_time,
                const tbox::Pointer<hier::PatchLevel> old_level=tbox::Pointer<hier::PatchLevel>(
                    NULL),
                const bool allocate_data=true);

        void resetHierarchyConfiguration (
                const tbox::Pointer<hier::PatchHierarchy> hierarchy,
                const int coarsest_level,
                const int finest_level);

        void applyGradientDetector (
                const tbox::Pointer<hier::PatchHierarchy> hierarchy,
                const int level_number,
                const double error_data_time,
                const int tag_index,
                const bool initial_time,
                const bool uses_richardson_extrapolation_too);

        /**
         * Serializable methods.
         */

        void putToDatabase(tbox::Pointer<tbox::Database> database);

        /**
         * Copy new field variable values back to timelevel 0.
         *
         * @param level The level we are working on.
         */
        void resetField(
                const tbox::Pointer<hier::PatchLevel> level);

        /**
         * Copy variables in revert_vars back to timelevel 0.
         *
         * @param level The level we are working on.
         * @param hierarchy The patch hierarchy we are working on.
         * @param current_time The current simulation time.
         */
        void revert(
                const tbox::Pointer<hier::PatchLevel> level);

        /**
         * Perform cell-centered and momentum advection.
         *
         * @param level The level we are working on.
         * @param hierarchy The patch hierarchy we are working on.
         * @param current_time The current simulation time.
         */
        void advection(
                const tbox::Pointer<hier::PatchLevel> level,
                const tbox::Pointer<hier::PatchHierarchy> hierarchy,
                double current_time);

        void fillBoundaries();

        /*
         * Print out the field summary
         */
        void printFieldSummary(
                double time,
                int step);

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
        tbox::Pointer<hier::VariableContext> d_new;
        /** Scratch variables. */
        tbox::Pointer<hier::VariableContext> d_scratch;
        /** "Old" timelevel 0 variables. */
        tbox::Pointer<hier::VariableContext> d_current;
        /** Context used to write ViSiT dumps. */
        tbox::Pointer<hier::VariableContext> d_plot_context;
        /**
         * @}
         */

        hier::ComponentSelector d_var_scratch_data;
        hier::ComponentSelector d_var_cur_data;
        hier::ComponentSelector d_var_new_data;

        tbox::List<tbox::Pointer<hier::Variable> > d_field_vars;
        tbox::List<tbox::Pointer<hier::Variable> > d_revert_vars;

        /**
         * @name Communication exchanges
         * 
         * Algorithms for various communication points in the application.
         *
         * @{
         */
        tbox::Pointer<xfer::RefineAlgorithm> d_bdry_fill_half_step;
        tbox::Pointer<xfer::RefineAlgorithm> d_bdry_fill_prime_halos;
        tbox::Pointer<xfer::RefineAlgorithm> d_bdry_fill_pre_lagrange;
        tbox::Pointer<xfer::RefineAlgorithm> d_bdry_fill_post_viscosity;
        tbox::Pointer<xfer::RefineAlgorithm> d_bdry_fill_pre_sweep1_cell;
        tbox::Pointer<xfer::RefineAlgorithm> d_bdry_fill_pre_sweep1_mom;
        tbox::Pointer<xfer::RefineAlgorithm> d_bdry_fill_pre_sweep2_mom;

        tbox::Pointer<xfer::RefineAlgorithm> d_fill_new_level;
        /**
         * @}
         */

        bool advect_x;

        tbox::Pointer<hier::PatchHierarchy> d_current_hierarchy;
};

#endif
