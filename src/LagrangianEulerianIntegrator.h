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
         * \defgroup Types Constants representing the type of a given variable.
         * @{
         */
        const static int NORMAL = 1;
        const static int FIELD = 2;
        const static int REVERT = 4;
        /**
         * @}
         */

        /**
         * \defgroup Exchanges Constants for the various exchange points, used when registering variables.
         *
         * @{
         */
        const static int NO_EXCH = 0;
        const static int PRIME_CELLS_EXCH = 1;
        const static int PRE_LAGRANGE_EXCH = 2;
        const static int POST_VISCOSITY_EXCH = 4;
        const static int HALF_STEP_EXCH = 8;
        const static int PRE_SWEEP_1_CELL_EXCH = 16;
        const static int PRE_SWEEP_1_MOM_EXCH = 32;
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

        /*
         * Serializable methods.
         */

        void putToDatabase(tbox::Pointer<tbox::Database> database);

        /*
         * Copy new field variable values back to timelevel 0.
         */
        void resetField(
                const tbox::Pointer<hier::PatchLevel> level);

        void revert(
                const tbox::Pointer<hier::PatchLevel> level);

        void advection(
                const tbox::Pointer<hier::PatchLevel> level,
                const tbox::Pointer<hier::PatchHierarchy> hierarchy,
                double current_time);

        void fillBoundaries();

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
         * Variable contexts.
         */

        tbox::Pointer<hier::VariableContext> d_new;
        tbox::Pointer<hier::VariableContext> d_scratch;
        tbox::Pointer<hier::VariableContext> d_current;
        tbox::Pointer<hier::VariableContext> d_plot_context;

        hier::ComponentSelector d_temp_var_scratch_data;
        hier::ComponentSelector d_temp_var_cur_data;
        hier::ComponentSelector d_temp_var_new_data;

        tbox::List<tbox::Pointer<hier::Variable> > d_field_vars;
        tbox::List<tbox::Pointer<hier::Variable> > d_revert_vars;

        tbox::Pointer<xfer::RefineAlgorithm> d_bdry_fill_pressure;
        tbox::Pointer<xfer::RefineAlgorithm> d_bdry_fill_prime_halos;
        tbox::Pointer<xfer::RefineAlgorithm> d_bdry_fill_pre_lagrange;

        bool advect_x;
};

#endif
