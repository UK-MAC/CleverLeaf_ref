//////////////////////////////////////////////////////////////////////////////
// Crown Copyright 2014 AWE, Copyright 2014 David Beckingsale.
//
// This file is part of CleverLeaf.
//
// CleverLeaf is free software: you can redistribute it and/or modify it under
// the terms of the GNU General Public License as published by the Free Software
// Foundation, either version 3 of the License, or (at your option) any later
// version.
//
// CleverLeaf is distributed in the hope that it will be useful, but WITHOUT ANY
// WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR
// A PARTICULAR PURPOSE. See the GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License along with
// CleverLeaf. If not, see http://www.gnu.org/licenses/.
//////////////////////////////////////////////////////////////////////////////
#ifndef CLEVERLEAF_LAGRANGIANEULERIANPATCHSTRATEGY_H_
#define CLEVERLEAF_LAGRANGIANEULERIANPATCHSTRATEGY_H_

#include "SAMRAI/hier/Patch.h"
#include "SAMRAI/hier/VariableContext.h"
#include "SAMRAI/tbox/Dimension.h"
#include "SAMRAI/xfer/RefinePatchStrategy.h"

using namespace SAMRAI;

class LagrangianEulerianLevelIntegrator;

/**
 * @class LagrangianEulerianPatchStrategy
 *
 * A class to describe the abstract methods needed by
 * LagrangianEulerianLevelIntegrator.
 *
 * The abstract methods here provide all the operations needed to integrate a
 * patch through the Lagrangian-Eulerian scheme seen in Cloverleaf. This
 * includes the Lagrangian-step methods to advance the solution, and the
 * advection methods for remapping the grid.
 */
class LagrangianEulerianPatchStrategy: public xfer::RefinePatchStrategy
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
     * @param dim The dimension of the problem.
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
     * @param initial_time True, if it is the initial problem time.
     */
    virtual void initializeDataOnPatch(
        hier::Patch& patch,
        double init_data_time,
        bool initial_time) = 0;

    /**
     * Calculate acceleration due to density and pressure. 
     *
     * @param patch The patch to work on.
     * @param dt The current dt value.
     */
    virtual void accelerate(hier::Patch& patch, double dt) = 0;

    /**
     * Calculate the safe timestep for the given patch.
     *
     * @param patch The patch to work on.
     *
     * @returns The safe dt value.
     */
    virtual double calc_dt_knl(hier::Patch& patch) = 0;

    /**
     * Apply the equation of state to compute pressure and soundspeed.
     *
     * @param patch The patch to work on.
     * @param predict If true, update new copies of the density and energy
     *                arrays.
     */
    virtual void ideal_gas_knl(hier::Patch& patch, bool predict) = 0;

    /**
     * Calculate viscosity on the given patch.
     *
     * @param patch The patch to work on.
     */
    virtual void viscosity_knl(hier::Patch& patch) = 0;

    /**
     * Calculate the updated energy and density values.
     *
     * @param patch The patch to work on.
     * @param dt The current dt value.
     * @param predict If true, update new copies of density and energy.
     */
    virtual void pdv_knl(hier::Patch& patch, double dt, bool predict) = 0;

    /**
     * Calculate the volume flux in the X and Y directions.
     *
     * @param patch The patch to work on.
     * @param dt The current dt value.
     */
    virtual void flux_calc_knl(hier::Patch& patch, double dt) = 0;

    /**
     * Compute the cell-centered advection.
     *
     * @param patch The patch to work on.
     * @param sweep_number The number of the sweep (1 or 2).
     * @param direction The direction to advect in.
     */
    virtual void advec_cell(
        hier::Patch& patch,
        int sweep_number,
        ADVEC_DIR direction) = 0;

    /**
     * Calculate the momentum advection.
     *
     * @param patch The patch to work on.
     * @param sweep_number The number of the sweep (1 or 2).
     * @param direction The direction to advect in.
     * @param which_vel The velocity to advect (X or Y).
     */
    virtual void advec_mom(
        hier::Patch& patch,
        int sweep_number,
        ADVEC_DIR direction,
        ADVEC_DIR which_vel) = 0;

    /**
     * Calculate field summary values for a patch.
     *
     * This method calculates the values for each of the summary quantities on
     * the patch.
     *
     * @param patch The patch to work on.
     * @param total_volume The total volume of the patch.
     * @param total_mass The total mass of the patch.
     * @param total_pressure The total pressure of the patch.
     * @param total_internal_energy The total internal energy of the patch.
     * @param total_kinetic_energy The total kinetic energy of the patch.
     */
    virtual void field_summary(
        hier::Patch& patch,
        double* total_volume,
        double* total_mass,
        double* total_pressure,
        double* total_internal_energy,
        double* total_kinetic_energy,
        int* total_effective_cells) = 0;

    /**
     * Method to call the debug method on a patch.
     *
     * This method can be called at an arbitrary point in the hydrodynamcis
     * cycle and then used to inspect the contect of the variables on the given
     * patch.
     *
     * @param patch The patch to work on.
     */
    virtual void debug_knl(hier::Patch& patch) = 0;

    /**
     * Get the data context corresponding to the current time.
     *
     * @returns Current data context.
     */
    std::shared_ptr<hier::VariableContext> getCurrentDataContext();

    /**
     * Get the data context corresponding to the new time.
     *
     * @returns New data context.
     */
    std::shared_ptr<hier::VariableContext> getNewDataContext();

    /**
     * Set the data context for the current time.
     *
     * This context contains the current data.
     *
     * @param context Context for current time.
     */
    void setCurrentDataContext(
        std::shared_ptr<hier::VariableContext> context);

    /**
     * Set the data context for the new time.
     *
     * This context contains the advanced data.
     *
     * @param context Context for new time.
     */
    void setNewDataContext(
        std::shared_ptr<hier::VariableContext> context);

    /**
     * Set the exchange flag variable.
     *
     * The exchange flag variable is used to control which of the variables
     * will have boundary conditions applied when the
     * setPhysicalBoundaryConditions method is called.
     *
     * @param exchange The exchange flag value.
     */
    void setExchangeFlag(const int exchange);

    /**
     * Get the dimension of the problem.
     *
     * For CleverLeaf, this is always 2.
     *
     * @returns The dimension of the problem.
     */
    const tbox::Dimension& getDim() const;

    virtual void tagGradientDetectorCells(
        hier::Patch& patch,
        const double regrid_time,
        const bool initial_error,
        const int tag_index) = 0;

    virtual void setPhysicalBoundaryConditions(
          hier::Patch& patch,
          const double fill_time,
          const hier::IntVector& ghost_width_to_fill) = 0;

    virtual hier::IntVector getRefineOpStencilWidth(
        const tbox::Dimension &dim ) const;

    virtual void preprocessRefine(
          hier::Patch& fine,
          const hier::Patch& coarse,
          const hier::Box& fine_box,
          const hier::IntVector& ratio);

    virtual void postprocessRefine(
          hier::Patch& fine,
          const hier::Patch& coarse,
          const hier::Box& fine_box,
          const hier::IntVector& ratio);

    virtual void fillLevelIndicator(
          hier::Patch& patch,
          const int level_number) = 0;
  protected:
    int d_which_exchange;
  private:
    const tbox::Dimension d_dim;

    std::shared_ptr<hier::VariableContext> d_new_data_context;
    std::shared_ptr<hier::VariableContext> d_current_data_context;
};
#endif
