#ifndef included_LagrangianEulerianPatchStrategy
#define included_LagrangianEulerianPatchStrategy

#include "SAMRAI/hier/Patch.h"
#include "SAMRAI/hier/VariableContext.h"
#include "SAMRAI/tbox/Dimension.h"
#include "SAMRAI/xfer/RefinePatchStrategy.h"

using namespace SAMRAI;

class LagrangianEulerianLevelIntegrator;

/**
 * @class LagrangianEulerianPatchStrategy
 *
 * A class to describe the abstract methods needed by LagrangianEulerianLevelIntegrator.
 *
 * The abstract methods here provide all the operations needed to integrate a
 * hierarchy of patches through the Lagrangian-Eulerian scheme seen in
 * Cloverleaf. This includes the Lagrangian step methods and the advection
 * methods for remapping the grid.
 */
class LagrangianEulerianPatchStrategy:
    public xfer::RefinePatchStrategy
{
    public:
        /**
         * Enum for the direction of advection.
         */
        enum ADVEC_DIR {
            X = 1,
            Y = 2
        };

        /**
         * Simple constructor to set the dimension of the problem.
         *
         * @param dim The dimension of the problem
         */
        LagrangianEulerianPatchStrategy(const tbox::Dimension& dim);

        /**
         * Register model variables with the LagrangianEulerianLevelIntegrator.
         *
         * @param integrator The integrator object being used.
         */
        virtual void registerModelVariables(
                LagrangianEulerianLevelIntegrator* integrator) = 0;

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

        virtual void field_summary(
                hier::Patch& patch,
                double* vol,
                double* mass,
                double* press,
                double* ie,
                double* ke) = 0;

        virtual void tagGradientDetectorCells(
                hier::Patch& patch,
                const double regrid_time,
                const bool initial_error,
                const int tag_index) = 0;

        virtual void debug_knl(hier::Patch& patch) = 0;

        /**
         * Get the data context corresponding to the current time.
         *
         * @returns Current data context.
         */
        boost::shared_ptr<hier::VariableContext> getCurrentDataContext();

        /**
         * Get the data context corresponding to the new time.
         *
         * @returns New data context.
         */
        boost::shared_ptr<hier::VariableContext> getNewDataContext();

        /**
         * Get the data context corresponding to the scratch storage.
         *
         * @returns Scratch data context.
         */
        boost::shared_ptr<hier::VariableContext> getScratchDataContext();

        boost::shared_ptr<hier::VariableContext> getScratchNewDataContext();

        /**
         * Set the data context for the current time.
         *
         * This context contains the timelevel 0 data.
         *
         * @param context Context for current time.
         */
        void setCurrentDataContext(
                boost::shared_ptr<hier::VariableContext> context);

        /**
         * Set the data context for the new time.
         *
         * This context contains the timelevel 1 data.
         *
         * @param context Context for new time.
         */
        void setNewDataContext(
                boost::shared_ptr<hier::VariableContext> context);

        /**
         * Set the data context for the scratch space.
         *
         * @param context Context for scratch space.
         */
        void setScratchDataContext(
                boost::shared_ptr<hier::VariableContext> context);

        void setScratchNewDataContext(
                boost::shared_ptr<hier::VariableContext> context);
        /**
         * Get the dimension of the problem.
         *
         * This is ALWAYS 2 at the moment.
         *
         * @returns The dimension of the problem.
         */
        const tbox::Dimension& getDim() const;

        /*
         * RefinePatchStrategy methods.
         */

        /**
         * Set user-defined boundary conditions at the physical domain boundary.
         */
        virtual void
            setPhysicalBoundaryConditions(
                    hier::Patch& patch,
                    const double fill_time,
                    const hier::IntVector& ghost_width_to_fill) = 0;

        /**
         * Return maximum stencil width needed for user-defined
         * data interpolation operations.  Default is to return
         * zero, assuming no user-defined operations provided.
         *
         * Note that this function is not pure virtual. It is given a
         * dummy implementation here so that users may ignore it when
         * inheriting from this class.
         */
        virtual hier::IntVector
            getRefineOpStencilWidth( const tbox::Dimension &dim ) const;

        /**
         * Pre- and post-processing routines for implementing user-defined
         * spatial interpolation routines applied to variables.  The
         * interpolation routines are used in the hyperbolic AMR algorithm
         * for filling patch ghost cells before advancing data on a level
         * and after regridding a level to fill portions of the new level
         * from some coarser level.  These routines are called automatically
         * from within patch boundary filling schedules; thus, some concrete
         * function matching these signatures must be provided in the user's
         * patch routines.  However, the routines only need to perform some
         * operations when "USER_DEFINED_REFINE" is given as the interpolation
         * method for some variable when the patch routines register variables
         * with the hyperbolic level integration algorithm, typically.  If the
         * user does not provide operations that refine such variables in either
         * of these routines, then they will not be refined.
         *
         * The order in which these operations are used in each patch
         * boundary filling schedule is:
         *
         * - \b (1) {Call user's preprocessRefine() routine.}
         * - \b (2) {Refine all variables with standard interpolation operators.}
         * - \b (3) {Call user's postprocessRefine() routine.}
         *
         * Note that these functions are not pure virtual. They are given
         * dummy implementations here so that users may ignore them when
         * inheriting from this class.
         */
        virtual void
            preprocessRefine(
                    hier::Patch& fine,
                    const hier::Patch& coarse,
                    const hier::Box& fine_box,
                    const hier::IntVector& ratio);

        ///
        virtual void
            postprocessRefine(
                    hier::Patch& fine,
                    const hier::Patch& coarse,
                    const hier::Box& fine_box,
                    const hier::IntVector& ratio);

        void setExchangeFlag(const int exchange);

    private:
        const tbox::Dimension d_dim;

        boost::shared_ptr<hier::VariableContext> d_new_data_context;
        boost::shared_ptr<hier::VariableContext> d_current_data_context;
        boost::shared_ptr<hier::VariableContext> d_scratch_data_context;
        boost::shared_ptr<hier::VariableContext> d_scratch_new_data_context;

    protected:
        int d_which_exchange;

};

#endif
