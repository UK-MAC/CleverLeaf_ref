#ifndef included_LagrangianEulerianPatchStrategy
#define included_LagrangianEulerianPatchStrategy

#include "SAMRAI/hier/Patch.h"
#include "SAMRAI/hier/VariableContext.h"
#include "SAMRAI/tbox/Dimension.h"
#include "SAMRAI/tbox/Pointer.h"

using namespace SAMRAI;

class LagrangianEulerianIntegrator;

/**
 * @class LagrangianEulerianPatchStrategy
 *
 * A class to describe the abstract methods needed by LagrangianEulerianIntegrator.
 *
 * The abstract methods here provide all the operations needed to integrate a
 * hierarchy of patches through the Lagrangian-Eulerian scheme seen in
 * Cloverleaf. This includes the Lagrangian step methods and the advection
 * methods for remapping the grid.
 */
class LagrangianEulerianPatchStrategy
{
    public:
        /**
         * Enum for the direction of advection.
         */
        enum ADVEC_DIR {
            X = 1,
            Y = 2 };

        /**
         * Simple constructor to set the dimension of the problem.
         *
         * @param dim The dimension of the problem
         */
        LagrangianEulerianPatchStrategy(const tbox::Dimension& dim);

        /**
         * Register model variables with the LagrangianEulerianIntegrator.
         *
         * @param integrator The integrator being used.
         */
        virtual void registerModelVariables(
                LagrangianEulerianIntegrator* integrator) = 0;

        /**
         * Set up data on the patch.
         *
         * @param patch The patch to initialize data on.
         * @param init_data_time The initial data time.
         * @param initial_time True if it is the initial problem time.
         */
        virtual void initializeDataOnPatch(
                hier::Patch& patch,
                double init_data_time,
                bool initial_time) = 0;

        /**
         * Calculate acceleration due to density and pressure. 
         *
         * @param patch The patch to work on.
         * @param dt Current dt value.
         */
        virtual void accelerate(
                hier::Patch& patch,
                double dt) = 0;

        /**
         * Calculate the safe timestep for the given patch.
         *
         * @param patch Patch to calculate timestep for.
         *
         * @returns The safe dt value.
         */
        virtual double calc_dt_knl(
                hier::Patch& patch) = 0;

        /**
         * Calculate equation of state to compute pressure and soundspeed.
         *
         * @param patch Patch to work on.
         * @param predict If true, update using tl 1 copies of the density and energy arrays.
         */
        virtual void ideal_gas_knl(
                hier::Patch& patch,
                bool predict) = 0;

        /**
         * Calculate viscosity on the give patch.
         *
         * @param patch Patch to work on.
         */
        virtual void viscosity_knl(
                hier::Patch& patch) = 0;

        /**
         * Calculate the updated energy and density values.
         *
         * @param patch Patch to work on.
         * @param dt Current dt value.
         * @param predict If true, update timelevel 1 copies of density and energy.
         */
        virtual void pdv_knl(
                hier::Patch& patch,
                double dt,
                bool predict) = 0;

        /**
         * Calculate the volume flux in the X and Y directions.
         *
         * @param patch Patch to work on.
         * @param dt Current dt value.
         */
        virtual void flux_calc_knl(
                hier::Patch& patch,
                double dt) = 0;

        /**
         * Compute the cell-centered advection.
         *
         * @param patch Patch to work on.
         * @param sweep_number Which sweep.
         * @param direction Which direction to advect in.
         */
        virtual void advec_cell(hier::Patch& patch,
                int sweep_number,
                ADVEC_DIR direction) = 0;

        /**
         * Calculate the momentum advection.
         *
         * @param patch Patch to work on.
         * @param sweep_number Which sweep.
         * @param direction Which direction to advect in.
         * @param which_vel Which velocity to advect (X or Y).
         */
        virtual void advec_mom(hier::Patch& patch,
                int sweep_number,
                ADVEC_DIR direction,
                ADVEC_DIR which_vel) = 0;

        virtual void tagGradientDetectorCells(
                hier::Patch& patch,
                const double regrid_time,
                const bool initial_error,
                const int tag_index);

        /**
         * Get the data context corresponding to the current time.
         *
         * @returns Current data context.
         */
        tbox::Pointer<hier::VariableContext> getCurrentDataContext();

        /**
         * Get the data context corresponding to the new time.
         *
         * @returns New data context.
         */
        tbox::Pointer<hier::VariableContext> getNewDataContext();

        /**
         * Set the data context for the current time.
         *
         * This context contains the timelevel 0 data.
         *
         * @param context Context for current time.
         */
        void setCurrentDataContext(
                tbox::Pointer<hier::VariableContext> context);

        /**
         * Set the data context for the new time.
         *
         * This context contains the timelevel 1 data.
         *
         * @param context Context for new time.
         */
        void setNewDataContext(
                tbox::Pointer<hier::VariableContext> context);

        /**
         * Get the dimension of the problem.
         *
         * This is ALWAYS 2 at the moment.
         *
         * @returns The dimension of the problem.
         */
        const tbox::Dimension& getDim() const;

    private:
        const tbox::Dimension d_dim;

        tbox::Pointer<hier::VariableContext> d_new_data_context;

        tbox::Pointer<hier::VariableContext> d_current_data_context;
};
#endif
