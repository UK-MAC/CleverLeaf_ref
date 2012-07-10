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

class LagrangianEulerianIntegrator:
    public algs::TimeRefinementLevelStrategy,
    public mesh::StandardTagAndInitStrategy,
    public tbox::Serializable
{
    public:

        enum VAR_TYPE {
            FIELD = 0,
            FLUX = 1,
            NORMAL = 2 };


        LagrangianEulerianIntegrator(
                const std::string& object_name,
                tbox::Pointer<tbox::Database> input_db,
                LagrangianEulerianPatchStrategy* patch_strategy
                );

        ~LagrangianEulerianIntegrator();

        /**
         * Register a variable with the integrator, allowing
         * it to be correctly transferred at halo exchanges as well as
         * coarsen/refine times
         *
         * @param var The variable to register
         * @param ghosts The number of ghosts this variable has
         * @param transfer_geom
         */
        void registerVariable(
                tbox::Pointer<hier::Variable> var,
                const VAR_TYPE var_type,
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
        void initializeLevelIntegrator(
                tbox::Pointer<mesh::GriddingAlgorithmStrategy> gridding_alg);

        double getLevelDt(
                const tbox::Pointer<hier::PatchLevel> level,
                const double dt_time,
                const bool initial_time);

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
                const tbox::Pointer<hier::PatchHierarchy> hierarchy);

    protected:
        /*
         * PatchStrategy contains user-specified methods needed for
         * operating on a patch in the AMR hierarchy.
         */
        LagrangianEulerianPatchStrategy* d_patch_strategy;

        /*
         * The name of this object.
         */
        std::string d_object_name;

        /*
         * The dimension of the problem.
         */
        const tbox::Dimension d_dim;

        /*
         * Variable contexts.
         */
        tbox::Pointer<hier::VariableContext> d_current;
        tbox::Pointer<hier::VariableContext> d_new;
        tbox::Pointer<hier::VariableContext> d_scratch;
        tbox::Pointer<hier::VariableContext> d_plot_context;

        hier::ComponentSelector d_temp_var_scratch_data;
        hier::ComponentSelector d_temp_var_cur_data;
        hier::ComponentSelector d_temp_var_new_data;

        tbox::List<tbox::Pointer<hier::Variable> > d_field_vars;
        tbox::List<tbox::Pointer<hier::Variable> > d_revert_vars;

        tbox::Pointer<xfer::RefineAlgorithm> d_bdry_fill_density;

        bool advect_x;
};

#endif
