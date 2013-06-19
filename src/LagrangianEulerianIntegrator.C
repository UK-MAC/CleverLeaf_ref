#include "LagrangianEulerianIntegrator.h"

LagrangianEulerianIntegrator::LagrangianEulerianIntegrator(
        const boost::shared_ptr<tbox::Database>& input_db,
        const boost::shared_ptr<hier::PatchHierarchy>& hierarchy,
        const boost::shared_ptr<LagrangianEulerianLevelIntegrator>& level_integrator,
        const boost::shared_ptr<mesh::GriddingAlgorithmStrategy>& gridding_algorithm):
    d_patch_hierarchy(hierarchy),
    d_level_integrator(level_integrator),
    d_gridding_algorithm(gridding_algorithm),
    d_integrator_step(0)
{

    d_regrid_interval = input_db->getIntegerWithDefault("regrid_interval", 1);
    d_start_time = input_db->getDoubleWithDefault("start_time", 0.0);
    d_end_time = input_db->getDoubleWithDefault("end_time", 1.0);
    d_end_step = input_db->getIntegerWithDefault("max_integrator_steps", 100000000);
    d_grow_dt = input_db->getDoubleWithDefault("grow_dt", 1.5);
    d_integrator_time = d_start_time;

    d_dt = 0.04;
}

double LagrangianEulerianIntegrator::initializeHierarchy()
{

    d_level_integrator->initializeLevelIntegrator(d_gridding_algorithm);

    d_gridding_algorithm->makeCoarsestLevel(d_start_time);

    /*
     * initialize level data here.
     */

    initializeLevelData(0);

    /*
     * After data on each level is initialized at simulation start time,
     * coarser levels are synchronized with finer levels that didn't exist
     * when the coarser level initial data was set.  This synchronization
     * process is defined by the integration algorithm.
     */
    if (d_patch_hierarchy->getFinestLevelNumber() > 0) {
        // "true" argument: const bool initial_time = true;
        d_level_integrator->synchronizeNewLevels(d_patch_hierarchy,
                0,
                d_patch_hierarchy->getFinestLevelNumber(),
                d_start_time,
                true);
    }

    /*
     * TODO: calculate min hierarchy dt
     */

    getMinHeirarchyDt(true);

    return d_dt;
}

void LagrangianEulerianIntegrator::initializeLevelData(const int level_number)
{
   const boost::shared_ptr<hier::PatchLevel> patch_level(
      d_patch_hierarchy->getPatchLevel(level_number));

   bool initial_time = true;

   if (d_patch_hierarchy->levelCanBeRefined(level_number)) {
      int tag_buffer = 4;

      double regrid_start_time =
          d_integrator_time - d_dt;

      d_gridding_algorithm->makeFinerLevel(
              tag_buffer,
              true,
              d_integrator_step,
              d_integrator_time,
              /*
               * TODO: check how this works!
               */
              regrid_start_time);

      /*
       * If new finer level is made, data on its patches is initialized.
       * Also, if new level can be refined and time integration is used
       * during regridding process, data on the current level is advanced
       * through a single time increment to provide boundary data for the
       * regridding process on finer levels.
       */

      if (d_patch_hierarchy->finerLevelExists(level_number)) {
         /*
          * RECURSIVE invocation of initialization on next finer level.
          */
         initializeLevelData(level_number + 1);
      }
   }
}

double LagrangianEulerianIntegrator::advanceHierarchy(const double dt)
{
   int finest_level_number = d_patch_hierarchy->getFinestLevelNumber();
   double dt_new = tbox::MathUtilities<double>::getMax();

   int level_num;

   for (level_num = 0; level_num <= finest_level_number; level_num++) {
       boost::shared_ptr<hier::PatchLevel> patch_level(d_patch_hierarchy->getPatchLevel(level_num));

       d_level_integrator->stampDataTime(patch_level, d_integrator_time);
   }

   for (level_num = 0; level_num <= finest_level_number; level_num++) {
       boost::shared_ptr<hier::PatchLevel> patch_level(d_patch_hierarchy->getPatchLevel(level_num));

       d_level_integrator->lagrangianPredictor(patch_level, dt);
   }

   for (level_num = 0; level_num <= finest_level_number; level_num++) {
      boost::shared_ptr<hier::PatchLevel> patch_level(d_patch_hierarchy->getPatchLevel(level_num));

     d_level_integrator->halfStepHaloExchange(patch_level,
            d_patch_hierarchy,
            d_integrator_time);
   }

   for (level_num = 0; level_num <= finest_level_number; level_num++) {
      boost::shared_ptr<hier::PatchLevel> patch_level(d_patch_hierarchy->getPatchLevel(level_num));

     d_level_integrator->lagrangianCorrector(patch_level, dt);
   }

   for (level_num = 0; level_num <= finest_level_number; level_num++) {
      boost::shared_ptr<hier::PatchLevel> patch_level(d_patch_hierarchy->getPatchLevel(level_num));

     d_level_integrator->preCellHaloExchange(patch_level,
            d_patch_hierarchy,
            d_integrator_time);
   }

   for (level_num = 0; level_num <= finest_level_number; level_num++) {
      boost::shared_ptr<hier::PatchLevel> patch_level(d_patch_hierarchy->getPatchLevel(level_num));

     d_level_integrator->advecCellSweep1(patch_level);
   }

   for (level_num = 0; level_num <= finest_level_number; level_num++) {
      boost::shared_ptr<hier::PatchLevel> patch_level(d_patch_hierarchy->getPatchLevel(level_num));

     d_level_integrator->preMomSweep1HaloExchange(patch_level,
            d_patch_hierarchy,
            d_integrator_time);
   }

   for (level_num = 0; level_num <= finest_level_number; level_num++) {
      boost::shared_ptr<hier::PatchLevel> patch_level(d_patch_hierarchy->getPatchLevel(level_num));

     d_level_integrator->advecMomSweep1(patch_level);
   }

   for (level_num = 0; level_num <= finest_level_number; level_num++) {
      boost::shared_ptr<hier::PatchLevel> patch_level(d_patch_hierarchy->getPatchLevel(level_num));

     d_level_integrator->advecCellSweep2(patch_level);
   }

   for (level_num = 0; level_num <= finest_level_number; level_num++) {

      boost::shared_ptr<hier::PatchLevel> patch_level(d_patch_hierarchy->getPatchLevel(level_num));

     d_level_integrator->preMomSweep2HaloExchange(patch_level,
            d_patch_hierarchy,
            d_integrator_time);
   }

   for (level_num = 0; level_num <= finest_level_number; level_num++) {
      boost::shared_ptr<hier::PatchLevel> patch_level(d_patch_hierarchy->getPatchLevel(level_num));

     d_level_integrator->advecMomSweep2(patch_level);
   }

   for (level_num = 0; level_num <= finest_level_number; level_num++) {
      boost::shared_ptr<hier::PatchLevel> patch_level(d_patch_hierarchy->getPatchLevel(level_num));

     d_level_integrator->resetField(patch_level);
   }

   getMinHeirarchyDt(false);

   d_integrator_time += dt;

   d_integrator_step++;

   int coarse_level_number = 0;

   if (finest_level_number > 0) {
      d_level_integrator->standardLevelSynchronization(
         d_patch_hierarchy,
         0,
         finest_level_number,
         d_integrator_time,
         d_integrator_time - dt);
   }

   /*
    * Are we ready to re-grid??
    */
   bool regrid_now = (d_integrator_step % d_regrid_interval == 0);

   if (regrid_now) {
      /*
       * Regrid all levels, from coarsest to finest.  If the error
       * estimation procedure uses time integration (e.g. Richardson
       * extrapolation) then we must supply the oldest time at which
       * data is stored.
       *
       * If the error coarsen ratio is two, data will be stored
       * from the previous timestep (at d_level_old_time).  If the
       * error coarsen ratio is three, data will be stored
       * from two previous timesteps (at d_level_old_old_time).
       *
       * If we are not using time integration, the oldest time
       * information should not be used, so it is set to NaNs
       * to throw an assertion if it is accessed.
       */

      tbox::Array<double> regrid_start_time;

      if (!d_gridding_algorithm->getTagAndInitializeStrategy()->
          usesTimeIntegration(d_integrator_step, d_integrator_time)) {

         int max_levels = d_patch_hierarchy->getMaxNumberOfLevels();
         regrid_start_time.resizeArray(max_levels);
         for (int i = 0; i < regrid_start_time.getSize(); i++) {
            regrid_start_time[i] = 0.;
         }

      } else {

//         if (d_gridding_algorithm->getTagAndInitializeStrategy()->getErrorCoarsenRatio() == 2) {
//            regrid_start_time = d_level_old_time;
//         } else if (d_gridding_algorithm->getTagAndInitializeStrategy()->getErrorCoarsenRatio() ==
//                    3) {
//            regrid_start_time = d_level_old_old_time;
//         } else {
//            TBOX_ERROR(
//               d_object_name << ": the supplied gridding "
//                             << "algorithm uses an error coarsen ratio of "
//                             << d_gridding_algorithm->
//               getTagAndInitializeStrategy()->getErrorCoarsenRatio()
//                             << " which is not supported in this class"
//                             << std::endl);
//         }

      }

      int max_levels = d_patch_hierarchy->getMaxNumberOfLevels();
      d_tag_buffer.resizeArray(max_levels);
      for (int i = 0; i < d_tag_buffer.getSize(); i++) {
          d_tag_buffer[i] = 2;
      }

      d_gridding_algorithm->
      regridAllFinerLevels(
         coarse_level_number,
         d_tag_buffer,
         d_integrator_step,
         d_integrator_time,
         regrid_start_time);

      /*
       * Synchronize data on new levels.
       */
      if (d_patch_hierarchy->getFinestLevelNumber() > 0) {

         const bool initial_time = false;
         d_level_integrator->synchronizeNewLevels(
            d_patch_hierarchy,
            coarse_level_number,
            d_patch_hierarchy->getFinestLevelNumber(),
            d_integrator_time,
            initial_time);
      }

   }

   return d_dt;
}

bool LagrangianEulerianIntegrator::stepsRemaining() const
{
    return (d_end_step - d_integrator_step) > 0;
}


void LagrangianEulerianIntegrator::getMinHeirarchyDt(const bool initial_time)
{
   int finest_level_number = d_patch_hierarchy->getFinestLevelNumber();
   int level_num;

   for (level_num = 0; level_num <= finest_level_number; level_num++) {
      boost::shared_ptr<hier::PatchLevel> patch_level(d_patch_hierarchy->getPatchLevel(level_num));
     d_level_integrator->timestepEoS(patch_level);
   }

   for (level_num = 0; level_num <= finest_level_number; level_num++) {
      boost::shared_ptr<hier::PatchLevel> patch_level(d_patch_hierarchy->getPatchLevel(level_num));
      d_level_integrator->primeBoundaryHaloExchange(patch_level, d_patch_hierarchy, d_integrator_time);
   }

   for (level_num = 0; level_num <= finest_level_number; level_num++) {
      boost::shared_ptr<hier::PatchLevel> patch_level(d_patch_hierarchy->getPatchLevel(level_num));

     d_level_integrator->viscosity(patch_level);
   }

   for (level_num = 0; level_num <= finest_level_number; level_num++) {
      boost::shared_ptr<hier::PatchLevel> patch_level(d_patch_hierarchy->getPatchLevel(level_num));

     d_level_integrator->postViscosityHaloExchange(patch_level, d_patch_hierarchy, d_integrator_time);
   }

   double level_dt;

   for (level_num = 0; level_num <= finest_level_number; level_num++) {
      boost::shared_ptr<hier::PatchLevel> patch_level(d_patch_hierarchy->getPatchLevel(level_num));

     level_dt = d_level_integrator->calcDt(patch_level);

     d_dt = tbox::MathUtilities<double>::Min(d_dt, level_dt);
   }
}
