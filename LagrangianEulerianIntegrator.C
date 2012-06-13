#include "LagrangianEulerianIntegrator.h"

void LagrangianEulerianIntegrator::initializeLevelIntegrator(
        tbox::Pointer<mesh::GriddingAlgorithmStrategy> gridding_alg);

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
   tbox::Pointer<hier::PatchLevel> patch_level(level);

   const tbox::SAMRAI_MPI& mpi(patch_level->getBoxLevel()->getMPI());

   double dt = tbox::MathUtilities<double>::getMax();

   if (!d_use_ghosts_for_dt) {

      d_patch_strategy->setDataContext(d_current);

      for (hier::PatchLevel::Iterator p(patch_level); p; p++) {
         tbox::Pointer<hier::Patch> patch = *p;

         patch->allocatePatchData(d_temp_var_scratch_data, dt_time);

         double patch_dt;
         patch_dt = d_patch_strategy->
            computeStableDtOnPatch(*patch,
               initial_time,
               dt_time);

         dt = tbox::MathUtilities<double>::Min(dt, patch_dt);
         //tbox::plog.precision(12);
         //tbox::plog << "Level " << patch_level->getLevelNumber()
         //           << " Patch " << p()
         //           << " box " << patch->getBox()
         //           << " has patch_dt " << patch_dt
         //           << " dt " << dt
         //           << std::endl;

         patch->deallocatePatchData(d_temp_var_scratch_data);
      }

      d_patch_strategy->clearDataContext();

   } else {

      //tbox::plog << "use ghosts for dt" << std::endl;

      patch_level->allocatePatchData(d_saved_var_scratch_data, dt_time);

      d_patch_strategy->setDataContext(d_scratch);

      d_bdry_sched_advance[patch_level->getLevelNumber()]->fillData(dt_time);

      for (hier::PatchLevel::Iterator ip(patch_level); ip; ip++) {
         tbox::Pointer<hier::Patch> patch = *ip;

         patch->allocatePatchData(d_temp_var_scratch_data, dt_time);

         double patch_dt;
         patch_dt = d_patch_strategy->
            computeStableDtOnPatch(*patch,
               initial_time,
               dt_time);

         dt = tbox::MathUtilities<double>::Min(dt, patch_dt);

         patch->deallocatePatchData(d_temp_var_scratch_data);
      }

      d_patch_strategy->clearDataContext();

      /*
       * Copy data from scratch to current and de-allocate scratch storage.
       * This may be excessive here, but seems necessary if the
       * computation of dt affects the state of the problem solution.
       * Also, this getLevelDt() routine is called at initialization only
       * in most cases.
       */

      copyTimeDependentData(patch_level, d_scratch, d_current);

      patch_level->deallocatePatchData(d_saved_var_scratch_data);
   }

   /*
    * The level time increment is a global min over all patches.
    */

   double global_dt = dt;

   if (mpi.getSize() > 1) {
      mpi.AllReduce(&global_dt, 1, MPI_MIN);
   }

   global_dt *= tbox::MathUtilities<double>::Min(d_cfl_init, d_cfl);

   return global_dt;

}

double LagrangianEulerianIntegrator::advanceLevel(
        const tbox::Pointer<hier::PatchLevel> level,
        const tbox::Pointer<hier::PatchHierarchy> hierarchy,
        const double current_time,
        const double new_time,
        const bool first_step,
        const bool last_step,
        const bool regrid_advance=false);

void LagrangianEulerianIntegrator::standardLevelSynchronization(
        const tbox::Pointer<hier::PatchHierarchy> hierarchy,
        const int coarsest_level,
        const int finest_level,
        const double sync_time,
        const double old_time);

void LagrangianEulerianIntegrator::synchronizeNewLevels(
        const tbox::Pointer<hier::PatchHierarchy> hierarchy,
        const int coarsest_level,
        const int finest_level,
        const double sync_time,
        const bool initial_time);

void LagrangianEulerianIntegrator::resetTimeDependentData(
        const tbox::Pointer<hier::PatchLevel> level,
        const double new_time,
        const bool can_be_refined);

void LagrangianEulerianIntegrator::resetDataToPreadvanceState(
        const tbox::Pointer<hier::PatchLevel> level);

bool LagrangianEulerianIntegrator::usingRefinedTimestepping()
{
    return false;
}


