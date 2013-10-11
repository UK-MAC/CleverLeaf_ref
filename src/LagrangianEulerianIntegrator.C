/*
 * Copyright 2013 David Beckingsale.
 * 
 * This file is part of CleverLeaf.
 * 
 * CleverLeaf is free software: you can redistribute it and/or modify it under 
 * the terms of the GNU General Public License as published by the 
 * Free Software Foundation, either version 3 of the License, or (at your option) 
 * any later version.
 * 
 * CleverLeaf is distributed in the hope that it will be useful, but 
 * WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or 
 * FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more 
 * details.
 * 
 * You should have received a copy of the GNU General Public License along with 
 * CleverLeaf. If not, see http://www.gnu.org/licenses/.
 */ 

#include "LagrangianEulerianIntegrator.h"

tbox::StartupShutdownManager::Handler LagrangianEulerianIntegrator::s_initialize_handler(
        LagrangianEulerianIntegrator::initializeCallback,
        0,
        0,
        LagrangianEulerianIntegrator::finalizeCallback,
        tbox::StartupShutdownManager::priorityTimers);

boost::shared_ptr<tbox::Timer> LagrangianEulerianIntegrator::t_initialize_hierarchy;
boost::shared_ptr<tbox::Timer> LagrangianEulerianIntegrator::t_advance_hierarchy;

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
    d_max_dt = input_db->getDoubleWithDefault("max_dt", tbox::MathUtilities<double>::getMax());
    d_fix_dt = input_db->getBoolWithDefault("fix_dt", false);

    d_integrator_time = d_start_time;

    d_dt = 0.04;
}

double LagrangianEulerianIntegrator::initializeHierarchy()
{

    t_initialize_hierarchy->start();

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

    getMinHeirarchyDt(true);

    t_initialize_hierarchy->stop();

    return d_dt;
}

void LagrangianEulerianIntegrator::initializeLevelData(const int level_number)
{
   const boost::shared_ptr<hier::PatchLevel> patch_level(
      d_patch_hierarchy->getPatchLevel(level_number));

   bool initial_time = true;

   if (d_patch_hierarchy->levelCanBeRefined(level_number)) {
      int tag_buffer = 2;

      double regrid_start_time =
          d_integrator_time - d_dt;

      d_gridding_algorithm->makeFinerLevel(
              tag_buffer,
              true,
              d_integrator_step,
              d_integrator_time,
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
    t_advance_hierarchy->start();

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

    d_level_integrator->swapAdvecDir();

    for (level_num = 0; level_num <= finest_level_number; level_num++) {
        boost::shared_ptr<hier::PatchLevel> patch_level(d_patch_hierarchy->getPatchLevel(level_num));

        d_level_integrator->resetField(patch_level);
    }

    getMinHeirarchyDt(false);

    d_integrator_time += dt;
    d_integrator_step++;

    if ((d_integrator_time + d_dt) > d_end_time) {
        d_dt = d_end_time - d_integrator_time;
    }

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

        int max_levels = d_patch_hierarchy->getMaxNumberOfLevels();
        d_tag_buffer.resize(max_levels);
        for (int i = 0; i < d_tag_buffer.size(); i++) {
            d_tag_buffer[i] = 2;
        }

        d_gridding_algorithm->
            regridAllFinerLevels(
                    coarse_level_number,
                    d_tag_buffer,
                    d_integrator_step,
                    d_integrator_time);

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

    t_advance_hierarchy->stop();

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
   double dt_old = d_dt;
   double dt = tbox::MathUtilities<double>::getMax();

   const tbox::SAMRAI_MPI& mpi(tbox::SAMRAI_MPI::getSAMRAIWorld());

   for (level_num = 0; level_num <= finest_level_number; level_num++) {
      boost::shared_ptr<hier::PatchLevel> patch_level(d_patch_hierarchy->getPatchLevel(level_num));
     d_level_integrator->timestepEoS(patch_level);
   }

   for (level_num = 0; level_num <= finest_level_number; level_num++) {
      boost::shared_ptr<hier::PatchLevel> patch_level(d_patch_hierarchy->getPatchLevel(level_num));
      d_level_integrator->preLagrangeHaloExchange(patch_level, d_patch_hierarchy, d_integrator_time);
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

     dt = tbox::MathUtilities<double>::Min(dt, level_dt);
   }

    if (mpi.getSize() > 1) {
        mpi.AllReduce(&dt, 1, MPI_MIN);
    }

   if(d_fix_dt) {
       d_dt = d_max_dt;
   } else {
       d_dt = tbox::MathUtilities<double>::Min(dt, (dt_old*d_grow_dt));
       d_dt = tbox::MathUtilities<double>::Min(d_dt, d_max_dt);
   }
}

void LagrangianEulerianIntegrator::initializeCallback()
{
    t_initialize_hierarchy = tbox::TimerManager::getManager()->
        getTimer("LagrangianEulerianIntegrator::initializeHierarchy()");
    t_advance_hierarchy = tbox::TimerManager::getManager()->
        getTimer("LagrangianEulerianIntegrator::advanceHierarchy()");
}

void LagrangianEulerianIntegrator::finalizeCallback()
{
    t_initialize_hierarchy.reset();
    t_advance_hierarchy.reset();
}
