#include "LagrangianEulerianIntegrator.h"

#include "SAMRAI/tbox/MathUtilities.h"
#include "SAMRAI/hier/VariableDatabase.h"
#include "SAMRAI/hier/Patch.h"

LagrangianEulerianIntegrator::LagrangianEulerianIntegrator(
        const std::string& object_name,
        tbox::Pointer<tbox::Database> input_db,
        LagrangianEulerianPatchStrategy* patch_strategy):
    d_dim(patch_strategy->getDim())
{
    d_object_name = object_name;

    d_patch_strategy = patch_strategy;

    /*
     * Default parameter values.
     */

    /*
     * Communication algorithms
     */

    /*
     * Variable contexts
     *
     * d_current corresponds to the current timelevel
     * d_new corresponds to timelevel 1
     */
    d_current = hier::VariableDatabase::getDatabase()->getContext("CURRENT");
    d_new = hier::VariableDatabase::getDatabase()->getContext("NEW");

    d_plot_context = d_current;
}

LagrangianEulerianIntegrator::~LagrangianEulerianIntegrator()
{
}


void LagrangianEulerianIntegrator::initializeLevelIntegrator(
        tbox::Pointer<mesh::GriddingAlgorithmStrategy> gridding_alg)
{
    /*
     * Register model variables here.
     */
    d_patch_strategy->registerModelVariables(this);
}

/*
 *************************************************************************
 *
 * Invoke dt calculation routines in patch strategy and take a min
 * over all patches on the level.  The result will be the max of the
 * next timestep on the level. If the boolean recompute_dt is true,
 * the max timestep on the level will be computed.  If it is false,
 * the method will simply access the latest dt stored in the time
 * refinement integrator.
 *
 *************************************************************************
 */
double LagrangianEulerianIntegrator::getLevelDt(
        const tbox::Pointer<hier::PatchLevel> level,
        const double dt_time,
        const bool initial_time)
{
    //   tbox::Pointer<hier::PatchLevel> patch_level(level);
    //
    //   const tbox::SAMRAI_MPI& mpi(patch_level->getBoxLevel()->getMPI());
    //
    //   double dt = tbox::MathUtilities<double>::getMax();
    //
    //   if (!d_use_ghosts_for_dt) {
    //
    //      d_patch_strategy->setDataContext(d_current);
    //
    //      for (hier::PatchLevel::Iterator p(patch_level); p; p++) {
    //         tbox::Pointer<hier::Patch> patch = *p;
    //
    //         patch->allocatePatchData(d_temp_var_scratch_data, dt_time);
    //
    //         double patch_dt;
    //         patch_dt = d_patch_strategy->
    //            computeStableDtOnPatch(*patch,
    //               initial_time,
    //               dt_time);
    //
    //         dt = tbox::MathUtilities<double>::Min(dt, patch_dt);
    //         //tbox::plog.precision(12);
    //         //tbox::plog << "Level " << patch_level->getLevelNumber()
    //         //           << " Patch " << p()
    //         //           << " box " << patch->getBox()
    //         //           << " has patch_dt " << patch_dt
    //         //           << " dt " << dt
    //         //           << std::endl;
    //
    //         patch->deallocatePatchData(d_temp_var_scratch_data);
    //      }
    //
    //      d_patch_strategy->clearDataContext();
    //
    //   } else {
    //
    //      //tbox::plog << "use ghosts for dt" << std::endl;
    //
    //      patch_level->allocatePatchData(d_saved_var_scratch_data, dt_time);
    //
    //      d_patch_strategy->setDataContext(d_scratch);
    //
    //      d_bdry_sched_advance[patch_level->getLevelNumber()]->fillData(dt_time);
    //
    //      for (hier::PatchLevel::Iterator ip(patch_level); ip; ip++) {
    //         tbox::Pointer<hier::Patch> patch = *ip;
    //
    //         patch->allocatePatchData(d_temp_var_scratch_data, dt_time);
    //
    //         double patch_dt;
    //         patch_dt = d_patch_strategy->
    //            computeStableDtOnPatch(*patch,
    //               initial_time,
    //               dt_time);
    //
    //         dt = tbox::MathUtilities<double>::Min(dt, patch_dt);
    //
    //         patch->deallocatePatchData(d_temp_var_scratch_data);
    //      }
    //
    //      d_patch_strategy->clearDataContext();
    //
    //      /*
    //       * Copy data from scratch to current and de-allocate scratch storage.
    //       * This may be excessive here, but seems necessary if the
    //       * computation of dt affects the state of the problem solution.
    //       * Also, this getLevelDt() routine is called at initialization only
    //       * in most cases.
    //       */
    //
    //      copyTimeDependentData(patch_level, d_scratch, d_current);
    //
    //      patch_level->deallocatePatchData(d_saved_var_scratch_data);
    //   }
    //
    //   /*
    //    * The level time increment is a global min over all patches.
    //    */
    //
    //   double global_dt = dt;
    //
    //   if (mpi.getSize() > 1) {
    //      mpi.AllReduce(&global_dt, 1, MPI_MIN);
    //   }
    //
    //   global_dt *= tbox::MathUtilities<double>::Min(d_cfl_init, d_cfl);
    //
    //   return global_dt;

    return 0.004;
}

double LagrangianEulerianIntegrator::advanceLevel(
        const tbox::Pointer<hier::PatchLevel> level,
        const tbox::Pointer<hier::PatchHierarchy> hierarchy,
        const double current_time,
        const double new_time,
        const bool first_step,
        const bool last_step,
        const bool regrid_advance)
{

    /*
     * This routine performs the required steps for the predictor and corrector,
     * and the advection. It mirrors the following routines occuring in Cloverleaf
     * as hydro.f90:
     *
     *     CALL PdV(.TRUE.)
     *     CALL accelerate()
     *     CALL PdV(.FALSE.)
     *     CALL flux_calc()
     *     CALL advection()
     *     CALL reset_field()
     *
     * All halo exchanges are handled by this routine, as well as ensuring correcting
     * stepping over all levels of the AMR hierarchy.
     */


    /*
     * TODO: Acceleration kernel.
     */


    /*
     * reset_field is used to copy density, energy and velocity
     * timelevel 1 values back to timelevel 0.
     */

    return 0.03;


}

void LagrangianEulerianIntegrator::standardLevelSynchronization(
        const tbox::Pointer<hier::PatchHierarchy> hierarchy,
        const int coarsest_level,
        const int finest_level,
        const double sync_time,
        const double old_time){}

void LagrangianEulerianIntegrator::synchronizeNewLevels(
        const tbox::Pointer<hier::PatchHierarchy> hierarchy,
        const int coarsest_level,
        const int finest_level,
        const double sync_time,
        const bool initial_time){}

void LagrangianEulerianIntegrator::resetTimeDependentData(
        const tbox::Pointer<hier::PatchLevel> level,
        const double new_time,
        const bool can_be_refined){}

void LagrangianEulerianIntegrator::resetDataToPreadvanceState(
        const tbox::Pointer<hier::PatchLevel> level){}

bool LagrangianEulerianIntegrator::usingRefinedTimestepping() const
{
    return false;
}


void LagrangianEulerianIntegrator::initializeLevelData (
        const tbox::Pointer<hier::PatchHierarchy> hierarchy,
        const int level_number,
        const double init_data_time,
        const bool can_be_refined,
        const bool initial_time,
        const tbox::Pointer<hier::PatchLevel> old_level,
        const bool allocate_data)
{

    tbox::Pointer<hier::PatchLevel> level(
            hierarchy->getPatchLevel(level_number));

    /* 
     * Allocate storage needed to initialize level and fill data
     * from coarser levels in AMR hierarchy, potentially. Since
     * time gets set when we allocate data, re-stamp it to current
     * time if we don't need to allocate.
     */
    if (allocate_data) {
        level->allocatePatchData(d_temp_var_cur_data, init_data_time); 
        level->allocatePatchData(d_temp_var_new_data, init_data_time); 
    } else {
        level->setTime(init_data_time, d_temp_var_cur_data); 
    }


    d_patch_strategy->setDataContext(d_current);

    for (hier::PatchLevel::Iterator p(level); p; p++) {
        tbox::Pointer<hier::Patch> patch(*p);

        patch->allocatePatchData(d_temp_var_scratch_data, init_data_time);

        d_patch_strategy->initializeDataOnPatch(*patch,
                init_data_time,
                initial_time);

        patch->deallocatePatchData(d_temp_var_scratch_data);
    }

    d_patch_strategy->clearDataContext();

}

void LagrangianEulerianIntegrator::resetHierarchyConfiguration (
        const tbox::Pointer<hier::PatchHierarchy> hierarchy,
        const int coarsest_level,
        const int finest_level){}

void LagrangianEulerianIntegrator::applyGradientDetector (
        const tbox::Pointer<hier::PatchHierarchy> hierarchy,
        const int level_number,
        const double error_data_time,
        const int tag_index,
        const bool initial_time,
        const bool uses_richardson_extrapolation_too){}

/*
 * Serializable methods.
 */

void LagrangianEulerianIntegrator::putToDatabase(
        tbox::Pointer<tbox::Database> database){}


void LagrangianEulerianIntegrator::registerVariable(
        tbox::Pointer<hier::Variable> var,
        const hier::IntVector ghosts,
        const tbox::Pointer<hier::GridGeometry> transfer_geom)
{
#ifdef DEBUG
    tbox::pout << "Registering variable: " << var->getName() << std::endl;
#endif

    const tbox::Dimension dim(ghosts.getDim());

    hier::VariableDatabase* variable_db = hier::VariableDatabase::getDatabase();

    const hier::IntVector& zero_ghosts(hier::IntVector::getZero(dim));

    int cur_id = variable_db->registerVariableAndContext(var,
            d_current,
            zero_ghosts);

    int new_id = variable_db->registerVariableAndContext(var,
            d_new,
            zero_ghosts);

    d_temp_var_cur_data.setFlag(cur_id);
    d_temp_var_new_data.setFlag(new_id);
}

tbox::Pointer<hier::VariableContext> LagrangianEulerianIntegrator::getPlotContext()
{
    return d_plot_context;
}
