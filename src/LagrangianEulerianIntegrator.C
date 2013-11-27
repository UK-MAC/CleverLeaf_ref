/*
 * Copyright 2013 David Beckingsale.
 * 
 * This file is part of CleverLeaf.
 * 
 * CleverLeaf is free software: you can redistribute it and/or modify it under
 * the terms of the GNU General Public License as published by the Free
 * Software Foundation, either version 3 of the License, or (at your option)
 * any later version.
 * 
 * CleverLeaf is distributed in the hope that it will be useful, but WITHOUT
 * ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
 * FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for
 * more details.
 * 
 * You should have received a copy of the GNU General Public License along with
 * CleverLeaf. If not, see http://www.gnu.org/licenses/.
 */ 

#include "LagrangianEulerianIntegrator.h"

tbox::StartupShutdownManager::Handler
LagrangianEulerianIntegrator::s_initialize_handler(
    LagrangianEulerianIntegrator::initializeCallback,
    0,
    0,
    LagrangianEulerianIntegrator::finalizeCallback,
    tbox::StartupShutdownManager::priorityTimers);

boost::shared_ptr<tbox::Timer>
LagrangianEulerianIntegrator::t_initialize_hierarchy;
boost::shared_ptr<tbox::Timer>
LagrangianEulerianIntegrator::t_advance_hierarchy;
boost::shared_ptr<tbox::Timer>
LagrangianEulerianIntegrator::t_synchronize_levels;
boost::shared_ptr<tbox::Timer>
LagrangianEulerianIntegrator::t_get_min_hierarchy_dt;

LagrangianEulerianIntegrator::LagrangianEulerianIntegrator(
    const boost::shared_ptr<tbox::Database>& input_db,
    const boost::shared_ptr<hier::PatchHierarchy>& hierarchy,
    const boost::shared_ptr<LagrangianEulerianLevelIntegrator>&
      level_integrator,
    const boost::shared_ptr<mesh::GriddingAlgorithmStrategy>&
      gridding_algorithm):
  d_patch_hierarchy(hierarchy),
  d_level_integrator(level_integrator),
  d_gridding_algorithm(gridding_algorithm),
  d_integrator_step(0)
{
  d_regrid_interval = input_db->getIntegerWithDefault("regrid_interval", 1);

  d_start_time = input_db->getDoubleWithDefault("start_time", 0.0);
  d_end_time = input_db->getDoubleWithDefault("end_time", 1.0);
  d_end_step = input_db->getIntegerWithDefault(
      "max_integrator_steps", 100000000);

  d_dt = input_db->getDoubleWithDefault("initial_dt", 0.04);
  d_grow_dt = input_db->getDoubleWithDefault("grow_dt", 1.5);
  d_max_dt = input_db->getDoubleWithDefault(
      "max_dt", tbox::MathUtilities<double>::getMax());
  d_fix_dt = input_db->getBoolWithDefault("fix_dt", false);

  d_integrator_time = d_start_time;
}

double LagrangianEulerianIntegrator::initializeHierarchy()
{
  t_initialize_hierarchy->start();

  d_level_integrator->initializeLevelIntegrator(d_gridding_algorithm);

  d_gridding_algorithm->makeCoarsestLevel(d_start_time);

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

    double regrid_start_time = d_integrator_time - d_dt;

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

  int level_number;
  int finest_level_number = d_patch_hierarchy->getFinestLevelNumber();
  double dt_new = tbox::MathUtilities<double>::getMax();

  for (level_number = 0; level_number <= finest_level_number; level_number++) {
    boost::shared_ptr<hier::PatchLevel> patch_level(
        d_patch_hierarchy->getPatchLevel(level_number));

    d_level_integrator->stampDataTime(patch_level, d_integrator_time);
  }

  for (level_number = 0; level_number <= finest_level_number; level_number++) {
    boost::shared_ptr<hier::PatchLevel> patch_level(
        d_patch_hierarchy->getPatchLevel(level_number));

    d_level_integrator->lagrangianPredictor(patch_level, dt);
  }

  for (level_number = 0; level_number <= finest_level_number; level_number++) {
    boost::shared_ptr<hier::PatchLevel> patch_level(
        d_patch_hierarchy->getPatchLevel(level_number));

    d_level_integrator->halfStepHaloExchange(
        patch_level,
        d_patch_hierarchy,
        d_integrator_time);
  }

  for (level_number = 0; level_number <= finest_level_number; level_number++) {
    boost::shared_ptr<hier::PatchLevel> patch_level(
        d_patch_hierarchy->getPatchLevel(level_number));

    d_level_integrator->lagrangianCorrector(patch_level, dt);
  }

  for (level_number = 0; level_number <= finest_level_number; level_number++) {
    boost::shared_ptr<hier::PatchLevel> patch_level(
        d_patch_hierarchy->getPatchLevel(level_number));

    d_level_integrator->preCellHaloExchange(
        patch_level,
        d_patch_hierarchy,
        d_integrator_time);
  }

  for (level_number = 0; level_number <= finest_level_number; level_number++) {
    boost::shared_ptr<hier::PatchLevel> patch_level(
        d_patch_hierarchy->getPatchLevel(level_number));

    d_level_integrator->advecCellSweep1(patch_level);
  }

  for (level_number = 0; level_number <= finest_level_number; level_number++) {
    boost::shared_ptr<hier::PatchLevel> patch_level(
        d_patch_hierarchy->getPatchLevel(level_number));

    d_level_integrator->preMomSweep1HaloExchange(
        patch_level,
        d_patch_hierarchy,
        d_integrator_time);
  }

  for (level_number = 0; level_number <= finest_level_number; level_number++) {
    boost::shared_ptr<hier::PatchLevel> patch_level(
        d_patch_hierarchy->getPatchLevel(level_number));

    d_level_integrator->advecMomSweep1(patch_level);
  }

  for (level_number = 0; level_number <= finest_level_number; level_number++) {
    boost::shared_ptr<hier::PatchLevel> patch_level(
        d_patch_hierarchy->getPatchLevel(level_number));

    d_level_integrator->advecCellSweep2(patch_level);
  }

  for (level_number = 0; level_number <= finest_level_number; level_number++) {

    boost::shared_ptr<hier::PatchLevel> patch_level(
        d_patch_hierarchy->getPatchLevel(level_number));

    d_level_integrator->preMomSweep2HaloExchange(
        patch_level,
        d_patch_hierarchy,
        d_integrator_time);
  }

  for (level_number = 0; level_number <= finest_level_number; level_number++) {
    boost::shared_ptr<hier::PatchLevel> patch_level(
        d_patch_hierarchy->getPatchLevel(level_number));

    d_level_integrator->advecMomSweep2(patch_level);
  }

  d_level_integrator->swapAdvecDir();

  for (level_number = 0; level_number <= finest_level_number; level_number++) {
    boost::shared_ptr<hier::PatchLevel> patch_level(
        d_patch_hierarchy->getPatchLevel(level_number));

    d_level_integrator->resetField(patch_level);
  }

  getMinHeirarchyDt(false);

  d_integrator_time += dt;
  d_integrator_step++;

  if ((d_integrator_time + d_dt) > d_end_time) {
    d_dt = d_end_time - d_integrator_time;
  }

  t_advance_hierarchy->stop();

  int coarse_level_number = 0;

  if (finest_level_number > 0) {
    t_synchronize_levels->start();

    d_level_integrator->standardLevelSynchronization(
        d_patch_hierarchy,
        0,
        finest_level_number,
        d_integrator_time,
        d_integrator_time - dt);

    t_synchronize_levels->stop();
  }

  /*
   * Are we ready to re-grid??
   */
  bool regrid_now = (d_integrator_step % d_regrid_interval == 0);

  if (regrid_now) {

    int max_levels = d_patch_hierarchy->getMaxNumberOfLevels();

    d_tag_buffer.resize(max_levels);

    for (int i = 0; i < d_tag_buffer.size(); i++) {
      d_tag_buffer[i] = d_regrid_interval;
    }

    d_gridding_algorithm->regridAllFinerLevels(
          coarse_level_number,
          d_tag_buffer,
          d_integrator_step,
          d_integrator_time);

    /*
     * Synchronize data on new levels.
     */
    if (d_patch_hierarchy->getFinestLevelNumber() > 0) {
      const bool initial_time = false;

      t_synchronize_levels->start();

      d_level_integrator->synchronizeNewLevels(
          d_patch_hierarchy,
          coarse_level_number,
          d_patch_hierarchy->getFinestLevelNumber(),
          d_integrator_time,
          initial_time);

      t_synchronize_levels->stop();
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
  t_get_min_hierarchy_dt->start();

  int finest_level_number = d_patch_hierarchy->getFinestLevelNumber();
  int level_number = 0;

  const tbox::SAMRAI_MPI& mpi(tbox::SAMRAI_MPI::getSAMRAIWorld());

  for (level_number = 0; level_number <= finest_level_number; level_number++) {
    boost::shared_ptr<hier::PatchLevel> patch_level(
        d_patch_hierarchy->getPatchLevel(level_number));

    d_level_integrator->timestepEoS(patch_level);
  }

  for (level_number = 0; level_number <= finest_level_number; level_number++) {
    boost::shared_ptr<hier::PatchLevel> patch_level(
        d_patch_hierarchy->getPatchLevel(level_number));

    d_level_integrator->preLagrangeHaloExchange(
        patch_level, d_patch_hierarchy, d_integrator_time);
  }

  for (level_number = 0; level_number <= finest_level_number; level_number++) {
    boost::shared_ptr<hier::PatchLevel> patch_level(
        d_patch_hierarchy->getPatchLevel(level_number));

    d_level_integrator->viscosity(patch_level);
  }

  for (level_number = 0; level_number <= finest_level_number; level_number++) {
    boost::shared_ptr<hier::PatchLevel> patch_level(
        d_patch_hierarchy->getPatchLevel(level_number));

    d_level_integrator->postViscosityHaloExchange(
        patch_level, d_patch_hierarchy, d_integrator_time);
  }

  double dt = tbox::MathUtilities<double>::getMax();

  for (level_number = 0; level_number <= finest_level_number; level_number++) {
    boost::shared_ptr<hier::PatchLevel> patch_level(
        d_patch_hierarchy->getPatchLevel(level_number));

    double level_dt = d_level_integrator->calcDt(patch_level);

    dt = tbox::MathUtilities<double>::Min(dt, level_dt);
  }

  if (mpi.getSize() > 1) {
    mpi.AllReduce(&dt, 1, MPI_MIN);
  }

  double dt_old = d_dt;

  if(!d_fix_dt) {
    d_dt = tbox::MathUtilities<double>::Min(dt, (dt_old*d_grow_dt));
    d_dt = tbox::MathUtilities<double>::Min(d_dt, d_max_dt);
  } else {
    d_dt = d_max_dt;
  }

  t_get_min_hierarchy_dt->stop();
}

double LagrangianEulerianIntegrator::printFieldSummary()
{
  double volume = 0.0;
  double mass = 0.0;
  double pressure = 0.0;
  double internal_energy = 0.0;
  double kinetic_energy = 0.0;

  double global_volume = 0.0;
  double global_mass = 0.0;
  double global_pressure = 0.0;
  double global_internal_energy = 0.0;
  double global_kinetic_energy = 0.0;
  double global_total_energy = 0.0;

  int finest_level_number = d_patch_hierarchy->getFinestLevelNumber();

  for(int level_number = 0; level_number <= finest_level_number; level_number++) {
    boost::shared_ptr<hier::PatchLevel> patch_level(
        d_patch_hierarchy->getPatchLevel(level_number));

    d_level_integrator->getFieldSummary(
        patch_level,
        &volume,
        &mass,
        &pressure,
        &internal_energy,
        &kinetic_energy);

    global_volume += volume;
    global_mass += mass;
    global_pressure += pressure;
    global_internal_energy += internal_energy;
    global_kinetic_energy += kinetic_energy;
  }

  global_total_energy = global_internal_energy + global_kinetic_energy;

  tbox::plog << std::setw(17) << " Volume"
    << std::setw(17) << " Mass"
    << std::setw(17) << " Density"
    << std::setw(17) << " Pressure"
    << std::setw(17) << " Internal Energy"
    << std::setw(17) << " Kinetic Energy"
    << std::setw(17) << " Total Energy" << std::endl;

  tbox::plog << std::scientific << std::setprecision(4)
    <<  std::setw(17) << global_volume
    <<  std::setw(17) << global_mass
    <<  std::setw(17) << global_mass/global_volume
    <<  std::setw(17) << global_pressure/global_volume
    <<  std::setw(17) << global_internal_energy
    <<  std::setw(17) << global_kinetic_energy
    <<  std::setw(17) << global_total_energy << std::endl;

  return global_kinetic_energy;
}

void LagrangianEulerianIntegrator::initializeCallback()
{
  t_initialize_hierarchy = tbox::TimerManager::getManager()->getTimer(
      "LagrangianEulerianIntegrator::initializeHierarchy()");
  t_advance_hierarchy = tbox::TimerManager::getManager()->getTimer(
      "LagrangianEulerianIntegrator::advanceHierarchy()");
  t_synchronize_levels = tbox::TimerManager::getManager()->getTimer(
      "LagrangianEulerianIntegrator::synchronizeLevels()");
  t_get_min_hierarchy_dt = tbox::TimerManager::getManager()->getTimer(
      "LagrangianEulerianIntegrator::getMinHierarchyDt()");
}

void LagrangianEulerianIntegrator::finalizeCallback()
{
  t_initialize_hierarchy.reset();
  t_advance_hierarchy.reset();
  t_synchronize_levels.reset();
  t_get_min_hierarchy_dt.reset();
}
