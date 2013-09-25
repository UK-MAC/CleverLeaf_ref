#include "LagrangianEulerianLevelIntegrator.h"

#include "SAMRAI/hier/Patch.h"
#include "SAMRAI/hier/VariableDatabase.h"
#include "SAMRAI/tbox/MathUtilities.h"
#include "SAMRAI/xfer/CoarsenSchedule.h"
#include "SAMRAI/xfer/RefineSchedule.h"
#include "SAMRAI/xfer/PatchLevelFillPattern.h"
#include "SAMRAI/xfer/PatchLevelBorderFillPattern.h"

#include <string>

tbox::StartupShutdownManager::Handler LagrangianEulerianLevelIntegrator::s_initialize_handler(
        LagrangianEulerianLevelIntegrator::initializeCallback,
        0,
        0,
        LagrangianEulerianLevelIntegrator::finalizeCallback,
        tbox::StartupShutdownManager::priorityTimers);

boost::shared_ptr<tbox::Timer> LagrangianEulerianLevelIntegrator::t_synchronize_levels_create;
boost::shared_ptr<tbox::Timer> LagrangianEulerianLevelIntegrator::t_synchronize_levels_fill;

boost::shared_ptr<tbox::Timer> LagrangianEulerianLevelIntegrator::t_fill_new_levels_create;
boost::shared_ptr<tbox::Timer> LagrangianEulerianLevelIntegrator::t_fill_new_levels_fill;

boost::shared_ptr<tbox::Timer> LagrangianEulerianLevelIntegrator::t_half_step_exchange_create;
boost::shared_ptr<tbox::Timer> LagrangianEulerianLevelIntegrator::t_half_step_exchange_fill;

boost::shared_ptr<tbox::Timer> LagrangianEulerianLevelIntegrator::t_prime_halos_exchange_create;
boost::shared_ptr<tbox::Timer> LagrangianEulerianLevelIntegrator::t_prime_halos_exchange_fill;

boost::shared_ptr<tbox::Timer> LagrangianEulerianLevelIntegrator::t_pre_lagrange_exchange_create;
boost::shared_ptr<tbox::Timer> LagrangianEulerianLevelIntegrator::t_pre_lagrange_exchange_fill;

boost::shared_ptr<tbox::Timer> LagrangianEulerianLevelIntegrator::t_post_viscosity_exchange_create;
boost::shared_ptr<tbox::Timer> LagrangianEulerianLevelIntegrator::t_post_viscosity_exchange_fill;

boost::shared_ptr<tbox::Timer> LagrangianEulerianLevelIntegrator::t_pre_sweep1_cell_exchange_create;
boost::shared_ptr<tbox::Timer> LagrangianEulerianLevelIntegrator::t_pre_sweep1_cell_exchange_fill;

boost::shared_ptr<tbox::Timer> LagrangianEulerianLevelIntegrator::t_pre_sweep1_mom_exchange_create;
boost::shared_ptr<tbox::Timer> LagrangianEulerianLevelIntegrator::t_pre_sweep1_mom_exchange_fill;

boost::shared_ptr<tbox::Timer> LagrangianEulerianLevelIntegrator::t_pre_sweep2_mom_exchange_create;
boost::shared_ptr<tbox::Timer> LagrangianEulerianLevelIntegrator::t_pre_sweep2_mom_exchange_fill;

boost::shared_ptr<tbox::Timer> LagrangianEulerianLevelIntegrator::t_tag_gradient_detector_cells;

LagrangianEulerianLevelIntegrator::LagrangianEulerianLevelIntegrator(
        const std::string& object_name,
        const boost::shared_ptr<tbox::Database>& input_db,
        LagrangianEulerianPatchStrategy* patch_strategy):
    d_dim(patch_strategy->getDim())
{
    d_object_name = object_name;

    d_patch_strategy = patch_strategy;

    /*
     * Default parameter values.
     */

    advect_x = true;

    /*
     * Communication algorithms
     */
    d_bdry_fill_half_step.reset(new xfer::RefineAlgorithm());
    d_bdry_fill_prime_halos.reset(new xfer::RefineAlgorithm());
    d_bdry_fill_pre_lagrange.reset(new xfer::RefineAlgorithm());
    d_bdry_fill_post_viscosity.reset(new xfer::RefineAlgorithm());
    d_bdry_fill_pre_sweep1_cell.reset(new xfer::RefineAlgorithm());
    d_bdry_fill_pre_sweep1_mom.reset(new xfer::RefineAlgorithm());
    d_bdry_fill_pre_sweep2_mom.reset(new xfer::RefineAlgorithm());

    d_fill_new_level.reset(new xfer::RefineAlgorithm());
    d_coarsen_field_data.reset(new xfer::CoarsenAlgorithm(d_dim));

    /*
     * Variable contexts
     *
     * d_current corresponds to the current timelevel
     * d_new corresponds to timelevel 1
     * d_scratch used during coarsen/refine (halo transfer)
     */
    d_current = hier::VariableDatabase::getDatabase()->getContext("CURRENT");
    d_new = hier::VariableDatabase::getDatabase()->getContext("NEW");
    d_scratch = hier::VariableDatabase::getDatabase()->getContext("SCRATCH");
    d_scratch_new = hier::VariableDatabase::getDatabase()->getContext("NEW_SCRATCH");

    /*
     * Pass these contexts up to the patch strategy
     */
    patch_strategy->setCurrentDataContext(d_current);
    patch_strategy->setNewDataContext(d_new);
    patch_strategy->setScratchDataContext(d_scratch);
    patch_strategy->setScratchNewDataContext(d_scratch_new);

    d_plot_context = d_current;
}

LagrangianEulerianLevelIntegrator::~LagrangianEulerianLevelIntegrator()
{
}

void LagrangianEulerianLevelIntegrator::initializeLevelIntegrator(
        const boost::shared_ptr<mesh::GriddingAlgorithmStrategy>& gridding_alg)
{
    /*
     * Register model variables here.
     */
    d_patch_strategy->registerModelVariables(this);
}

double LagrangianEulerianLevelIntegrator::getMaxFinerLevelDt(
        const int finer_level_number,
        const double coarse_dt,
        const hier::IntVector& ratio)
{

    return coarse_dt / double(ratio.max());
}

void LagrangianEulerianLevelIntegrator::standardLevelSynchronization(
        const boost::shared_ptr<hier::PatchHierarchy>& hierarchy,
        const int coarsest_level,
        const int finest_level,
        const double sync_time,
        const double old_time)
{
    boost::shared_ptr<xfer::CoarsenSchedule> coarsen_schedule;

    d_current_hierarchy = hierarchy;

    for (int fine_ln = finest_level; fine_ln > coarsest_level; fine_ln--) {
        const int coarse_ln = fine_ln - 1;

        boost::shared_ptr<hier::PatchLevel> fine_level = hierarchy->getPatchLevel(fine_ln);
        boost::shared_ptr<hier::PatchLevel> coarse_level = hierarchy->getPatchLevel(coarse_ln);

        const boost::shared_ptr<hier::PatchHierarchy> patch_hierarchy = hierarchy;

        d_patch_strategy->setExchangeFlag(FIELD_EXCH);

        t_synchronize_levels_create->start();
        coarsen_schedule = 
            d_coarsen_field_data->createSchedule(
                    coarse_level,
                    fine_level);
        t_synchronize_levels_create->stop();

        t_synchronize_levels_fill->start();
        coarsen_schedule->coarsenData();
        t_synchronize_levels_fill->stop();
    }
}

void LagrangianEulerianLevelIntegrator::synchronizeNewLevels(
        const boost::shared_ptr<hier::PatchHierarchy>& hierarchy,
        const int coarsest_level,
        const int finest_level,
        const double sync_time,
        const bool initial_time)
{

    d_current_hierarchy = hierarchy;
}

void LagrangianEulerianLevelIntegrator::initializeLevelData (
        const boost::shared_ptr<hier::PatchHierarchy>& hierarchy,
        const int level_number,
        const double init_data_time,
        const bool can_be_refined,
        const bool initial_time,
        const boost::shared_ptr<hier::PatchLevel>& old_level,
        const bool allocate_data)
{
    d_current_hierarchy = hierarchy;

    boost::shared_ptr<hier::PatchLevel> level(
            hierarchy->getPatchLevel(level_number));

    boost::shared_ptr<xfer::RefineSchedule> refine_schedule;

    /* 
     * Allocate storage needed to initialize level and fill data
     * from coarser levels in AMR hierarchy, potentially. Since
     * time gets set when we allocate data, re-stamp it to current
     * time if we don't need to allocate.
     */
    //if (allocate_data) {
        level->allocatePatchData(d_var_cur_data, init_data_time); 
        level->allocatePatchData(d_var_new_data, init_data_time); 
    //} else {
        //level->setTime(init_data_time, d_var_cur_data); 
        //level->setTime(init_data_time, d_var_new_data); 
    //}

   level->allocatePatchData(d_var_scratch_data, init_data_time);
   level->allocatePatchData(d_var_scratch_new_data, init_data_time);

   const tbox::SAMRAI_MPI& mpi(level->getBoxLevel()->getMPI());

   if ((level_number > 0) || old_level) {

      d_patch_strategy->setExchangeFlag(FIELD_EXCH);

      const boost::shared_ptr<hier::PatchHierarchy> patch_hierarchy(hierarchy);

      t_fill_new_levels_create->start();
      refine_schedule = 
          d_fill_new_level->createSchedule(level,
                  old_level,
                  level_number-1,
                  hierarchy,
                  d_patch_strategy);
      t_fill_new_levels_create->stop();

      t_fill_new_levels_fill->start();
      refine_schedule->fillData(init_data_time);
      t_fill_new_levels_fill->stop();
   }

    for (hier::PatchLevel::iterator p(level->begin()); p != level->end(); ++p) {
        boost::shared_ptr<hier::Patch> patch(*p);

        d_patch_strategy->initializeDataOnPatch(*patch,
                init_data_time,
                initial_time);
    }


    d_patch_strategy->setExchangeFlag(PRIME_CELLS_EXCH);

    if(initial_time) {
        t_prime_halos_exchange_create->start();
        refine_schedule = 
            d_bdry_fill_prime_halos->createSchedule(level, d_patch_strategy);
        t_prime_halos_exchange_create->stop();

        t_prime_halos_exchange_fill->start();
        refine_schedule->fillData(init_data_time);
        t_prime_halos_exchange_fill->stop();
    } else {
        if ((level_number > 0) || old_level) {
            level->allocatePatchData(d_var_scratch_data, init_data_time);
            level->allocatePatchData(d_var_scratch_new_data, init_data_time);

            const boost::shared_ptr<hier::PatchHierarchy> patch_hierarchy(hierarchy);

            t_prime_halos_exchange_create->start();
            refine_schedule =
                d_bdry_fill_prime_halos->createSchedule(
                        level,
                        old_level,
                        level_number-1,
                        hierarchy,
                        d_patch_strategy);
            t_prime_halos_exchange_create->stop();

            t_prime_halos_exchange_fill->start();
            refine_schedule->fillData(init_data_time);
            t_prime_halos_exchange_fill->stop();
        } else {
            t_prime_halos_exchange_create->start();
            refine_schedule = 
                d_bdry_fill_prime_halos->createSchedule(level, d_patch_strategy);
            t_prime_halos_exchange_create->stop();

            t_prime_halos_exchange_fill->start();
            refine_schedule->fillData(init_data_time);
            t_prime_halos_exchange_fill->stop();
        }
    }
}

void LagrangianEulerianLevelIntegrator::resetHierarchyConfiguration (
        const boost::shared_ptr<hier::PatchHierarchy>& hierarchy,
        const int coarsest_level,
        const int finest_level)
{
    int finest_level_number = hierarchy->getFinestLevelNumber();

    d_half_step_schedules.resize(finest_level_number+1);
    d_prime_halos_schedules.resize(finest_level_number+1);
    d_pre_lagrange_schedules.resize(finest_level_number+1);
    d_post_viscosity_schedules.resize(finest_level_number+1);
    d_pre_sweep1_cell_schedules.resize(finest_level_number+1);
    d_pre_sweep1_mom_schedules.resize(finest_level_number+1);
    d_pre_sweep2_mom_schedules.resize(finest_level_number+1);

    for(int level_number = coarsest_level; level_number <= finest_level_number; level_number++) {
        boost::shared_ptr<hier::PatchLevel> level(hierarchy->getPatchLevel(level_number));

        t_half_step_exchange_create->start();
        d_half_step_schedules[level_number] = 
            d_bdry_fill_half_step->createSchedule(
                    level,
                    level_number-1,
                    hierarchy,
                    d_patch_strategy);
        t_half_step_exchange_create->stop();

        t_prime_halos_exchange_create->start();
        d_prime_halos_schedules[level_number] =
            d_bdry_fill_prime_halos->createSchedule(
                    level,
                    level_number-1,
                    hierarchy,
                    d_patch_strategy);
        t_prime_halos_exchange_create->stop();

        t_pre_lagrange_exchange_create->start();
        d_pre_lagrange_schedules[level_number] = 
            d_bdry_fill_pre_lagrange->createSchedule(
                    level,
                    level_number-1,
                    hierarchy,
                    d_patch_strategy);
        t_pre_lagrange_exchange_create->stop();

        t_post_viscosity_exchange_create->start();
        d_post_viscosity_schedules[level_number] = 
            d_bdry_fill_post_viscosity->createSchedule(
                    level,
                    level_number-1,
                    hierarchy,
                    d_patch_strategy);
        t_post_viscosity_exchange_create->stop();

        t_pre_sweep1_cell_exchange_create->start();
        d_pre_sweep1_cell_schedules[level_number] = 
            d_bdry_fill_pre_sweep1_cell->createSchedule(
                    level,
                    level_number-1,
                    hierarchy,
                    d_patch_strategy);
        t_pre_sweep1_cell_exchange_create->stop();

        t_pre_sweep1_mom_exchange_create->start();
        d_pre_sweep1_mom_schedules[level_number] =
            d_bdry_fill_pre_sweep1_mom->createSchedule(
                    level,
                    level_number-1,
                    hierarchy,
                    d_patch_strategy);
        t_pre_sweep1_mom_exchange_create->stop();

        t_pre_sweep2_mom_exchange_create->start();
        d_pre_sweep2_mom_schedules[level_number] =
            d_bdry_fill_pre_sweep2_mom->createSchedule(
                    level,
                    level_number-1,
                    hierarchy,
                    d_patch_strategy);
        t_pre_sweep2_mom_exchange_create->stop();
    }

}

void LagrangianEulerianLevelIntegrator::applyGradientDetector (
        const boost::shared_ptr<hier::PatchHierarchy>& hierarchy,
        const int level_number,
        const double error_data_time,
        const int tag_index,
        const bool initial_time,
        const bool uses_richardson_extrapolation_too)
{
   boost::shared_ptr<hier::PatchLevel> level(hierarchy->getPatchLevel(level_number));

   const tbox::SAMRAI_MPI& mpi(level->getBoxLevel()->getMPI());

   t_tag_gradient_detector_cells->start();

    for (hier::PatchLevel::iterator ip(level->begin()); ip != level->end(); ++ip) {
      boost::shared_ptr<hier::Patch> patch(*ip);
      d_patch_strategy->
      tagGradientDetectorCells(*patch,
         error_data_time,
         initial_time,
         tag_index);
   }

   t_tag_gradient_detector_cells->stop();
}

void LagrangianEulerianLevelIntegrator::registerVariable(
        const boost::shared_ptr<hier::Variable>& var,
        const int var_type,
        const int var_exchanges,
        const hier::IntVector ghosts,
        const boost::shared_ptr<hier::BaseGridGeometry>& transfer_geom)
{
    const tbox::Dimension dim(ghosts.getDim());

    hier::VariableDatabase* variable_db = hier::VariableDatabase::getDatabase();

    const hier::IntVector& zero_ghosts(hier::IntVector::getZero(dim));

    if((var_type & FIELD) == FIELD) {
        d_field_vars.push_back(var);
    }

    if((var_type & REVERT) == REVERT) {
        //std::cout << "Found a revert var: " << var->getName() << std::endl;
        d_revert_vars.push_back(var);
    }

    int cur_id = variable_db->registerVariableAndContext(var,
            d_current,
            ghosts);

    int new_id = variable_db->registerVariableAndContext(var,
            d_new,
            ghosts);

    int scr_id = variable_db->registerVariableAndContext(var,
            d_scratch,
            ghosts);

    int scr_new_id = variable_db->registerVariableAndContext(var,
            d_scratch_new,
            ghosts);

    boost::shared_ptr<hier::RefineOperator> refine_op;
    boost::shared_ptr<hier::CoarsenOperator> coarsen_op;



//    if (var->getName() == "massflux") {
//        //tbox::pout << "Using CONSERVATIVE_LINEAR_REFINE..." << std::endl;
//        refine_op = transfer_geom->lookupRefineOperator(var, "CONSERVATIVE_LINEAR_REFINE");
//    } else if (var->getName() == "volflux") {
//        //tbox::pout << "Using CONSERVATIVE_LINEAR_REFINE..." << std::endl;
//        refine_op = transfer_geom->lookupRefineOperator(var, "CONSERVATIVE_LINEAR_REFINE");
//    } else if (var->getName() == "velocity") {
//        //tbox::pout << "Using LINEAR_REFINE..." << std::endl;
//        refine_op = transfer_geom->lookupRefineOperator(var, "LINEAR_REFINE");
//    }

    if(var->getName() == "velocity" || var->getName()=="vertexdeltas" || var->getName()=="vertexcoords") {
        refine_op = transfer_geom->lookupRefineOperator(var, "LINEAR_REFINE");
        coarsen_op = transfer_geom->lookupCoarsenOperator(var, "CONSTANT_COARSEN");

    } else {
        refine_op = transfer_geom->lookupRefineOperator(var, "CONSERVATIVE_LINEAR_REFINE");
    }

    if((var_type & FIELD) == FIELD) {

        d_fill_new_level->registerRefine(
                cur_id,
                cur_id,
                scr_id,
                refine_op);

        if (var->getName() == "velocity") {
            d_coarsen_field_data->registerCoarsen(
                    cur_id,
                    cur_id,
                    coarsen_op);
        } else if (var->getName() == "density") {
                d_coarsen_field_data->registerCoarsen(
                        cur_id,
                        cur_id,
                        transfer_geom->lookupCoarsenOperator(var, "VOLUME_WEIGHTED_COARSEN"));
        } else if (var->getName() == "energy") {
            d_coarsen_field_data->registerCoarsen(
                    cur_id,
                    cur_id,
                    transfer_geom->lookupCoarsenOperator(var, "MASS_WEIGHTED_COARSEN"));
        }
    }

    if((var_exchanges & PRIME_CELLS_EXCH) == PRIME_CELLS_EXCH) {
#ifdef DEBUG
        tbox::pout << "Registering " << var->getName() << " for initial exchange..." << std::endl;
#endif

        d_bdry_fill_prime_halos->registerRefine(
                cur_id,
                cur_id,
                scr_id,
                refine_op);

        if(var->getName() == "density" ||
                var->getName() == "energy" ||
                var->getName() == "velocity") {
            d_bdry_fill_prime_halos->registerRefine(
                    new_id,
                    new_id,
                    scr_new_id,
                    refine_op);
        }
    }

    if((var_exchanges & PRE_LAGRANGE_EXCH) == PRE_LAGRANGE_EXCH) {
#ifdef DEBUG
        tbox::pout << "Registering " << var->getName() << " for pre-lagrange exchange..." << std::endl;
#endif
        d_bdry_fill_pre_lagrange->registerRefine(
                cur_id,
                cur_id,
                scr_id,
                refine_op);
    }
               
    if((var_exchanges & POST_VISCOSITY_EXCH) == POST_VISCOSITY_EXCH) {
#ifdef DEBUG
        tbox::pout << "Registering " << var->getName() << " for post-viscosity exchange..." << std::endl;
#endif
        d_bdry_fill_post_viscosity->registerRefine(
                cur_id,
                cur_id,
                scr_id,
                refine_op);
    }

    if((var_exchanges & HALF_STEP_EXCH) == HALF_STEP_EXCH) {
#ifdef DEBUG
        tbox::pout << "Registering " << var->getName() << " for half-step exchange..." << std::endl;
#endif

        d_bdry_fill_half_step->registerRefine(
                cur_id,
                cur_id,
                scr_id,
                refine_op);
    }

    if((var_exchanges & PRE_SWEEP_1_CELL_EXCH) == PRE_SWEEP_1_CELL_EXCH) {

#ifdef DEBUG
        tbox::pout << "Registering " << var->getName() << " for pre-sweep1 cell exchange" << std::endl;
#endif
        if(var->getName() == "volflux") {
            d_bdry_fill_pre_sweep1_cell->registerRefine(
                    cur_id,
                    cur_id,
                    scr_id,
                    refine_op);
        } else {
            d_bdry_fill_pre_sweep1_cell->registerRefine(
                    new_id,
                    new_id,
                    scr_new_id,
                    refine_op);
        }
    }

    if((var_exchanges & PRE_SWEEP_1_MOM_EXCH) == PRE_SWEEP_1_MOM_EXCH) {
#ifdef DEBUG
        tbox::pout << "Registering " << var->getName() << " for pre-sweep1 mom exchange" << std::endl;
#endif
        if(var->getName() == "massflux") {
        d_bdry_fill_pre_sweep1_mom->registerRefine(
                cur_id,
                cur_id,
                scr_id,
                refine_op);
        } else {
            d_bdry_fill_pre_sweep1_mom->registerRefine(
                    new_id,
                    new_id,
                    scr_new_id,
                    refine_op);
        }
    }

    if((var_exchanges & PRE_SWEEP_2_MOM_EXCH) == PRE_SWEEP_2_MOM_EXCH) {
#ifdef DEBUG
        tbox::pout << "Registering " << var->getName() << " for pre-sweep1 mom exchange" << std::endl;
#endif
        if(var->getName() == "massflux") {
            d_bdry_fill_pre_sweep2_mom->registerRefine(
                    cur_id,
                    cur_id,
                    scr_id,
                    refine_op);
        } else {
            d_bdry_fill_pre_sweep2_mom->registerRefine(
                    new_id,
                    new_id,
                    scr_new_id,
                    refine_op);
        }
    }

    d_var_cur_data.setFlag(cur_id);
    d_var_new_data.setFlag(new_id);
    d_var_scratch_data.setFlag(scr_id);
    d_var_scratch_new_data.setFlag(scr_new_id);
}

boost::shared_ptr<hier::VariableContext> LagrangianEulerianLevelIntegrator::getPlotContext()
{
    return d_plot_context;
}


void LagrangianEulerianLevelIntegrator::resetField(
   const boost::shared_ptr<hier::PatchLevel>& level)
{
   for (hier::PatchLevel::iterator ip(level->begin()); ip != level->end(); ++ip) {
      boost::shared_ptr<hier::Patch> patch = *ip;

      std::list<boost::shared_ptr<hier::Variable> >::iterator
         field_var = d_field_vars.begin();

      while (field_var != d_field_vars.end()) {

          //tbox::pout << "Copying " << field_var()->getName() << " back to tl0" << std::endl;

         boost::shared_ptr<hier::PatchData> src_data =
            patch->getPatchData(*field_var, d_new);
         boost::shared_ptr<hier::PatchData> dst_data =
            patch->getPatchData(*field_var, d_current);

         dst_data->copy(*src_data);

         field_var++;
      }
   }
}

void LagrangianEulerianLevelIntegrator::revert(
   const boost::shared_ptr<hier::PatchLevel>& level)
{
   for (hier::PatchLevel::iterator ip(level->begin()); ip != level->end(); ++ip) {
      boost::shared_ptr<hier::Patch> patch = *ip;

      std::list<boost::shared_ptr<hier::Variable> >::iterator
         revert_var = d_revert_vars.begin();

      while (revert_var != d_revert_vars.end()) {
#if DEBUG
          tbox::pout << "Copying " << revert_var()->getName() << " back to tl0" << std::endl;
#endif
         boost::shared_ptr<hier::PatchData> dst_data =
            patch->getPatchData(*revert_var, d_new);
         boost::shared_ptr<hier::PatchData> src_data =
            patch->getPatchData(*revert_var, d_current);

        dst_data->copy(*src_data);
         revert_var++;
      }
   }
}

void LagrangianEulerianLevelIntegrator::printFieldSummary(
        double time,
        int step)
{
    boost::shared_ptr<hier::PatchLevel> level = d_current_hierarchy->getPatchLevel(0);
    const tbox::SAMRAI_MPI& mpi(level->getBoxLevel()->getMPI());

    double vol;
    double mass;
    double press;
    double ie;
    double ke;

    double global_vol;
    double global_mass;
    double global_press;
    double global_ie;
    double global_ke;

    for(hier::PatchLevel::iterator p(level->begin()); p != level->end(); ++p){

        boost::shared_ptr<hier::Patch>patch=*p;

        d_patch_strategy->ideal_gas_knl(*patch, false);

        d_patch_strategy->field_summary(*patch, &vol, &mass, &press, &ie, &ke);
    }

    global_vol = vol;
    global_mass = mass;
    global_press = press;
    global_ie = ie;
    global_ke = ke;

    if (mpi.getSize() > 1) {
        mpi.AllReduce(&global_vol, 1, MPI_SUM);
        mpi.AllReduce(&global_mass, 1, MPI_SUM);
        mpi.AllReduce(&global_press, 1, MPI_SUM);
        mpi.AllReduce(&global_ie, 1, MPI_SUM);
        mpi.AllReduce(&global_ke, 1, MPI_SUM);
    }

    if (mpi.getRank() == 0) {
        printf("%13s%16s %16s %16s %16s %16s %16s %16s\n", " ", "Volume", "Mass", "Density", "Pressure", "Internal Energy", "Kinetic Energy", "Total Energy");
        printf("%6s %7d %16.4E %16.4E %16.4E %16.4E %16.4E %16.4E %16.4E\n", "step:", step, global_vol, global_mass, global_mass/global_vol, global_press/global_vol, global_ie, global_ke, global_ie+global_ke);

    }
}

void LagrangianEulerianLevelIntegrator::lagrangianPredictor(
        const boost::shared_ptr<hier::PatchLevel>& level,
        const double dt)
{
    for (hier::PatchLevel::iterator p(level->begin()); p != level->end(); ++p) {

        boost::shared_ptr<hier::Patch>patch=*p;

        d_patch_strategy->pdv_knl(*patch,dt, true);
    }

    for (hier::PatchLevel::iterator p(level->begin()); p != level->end(); ++p) {

        boost::shared_ptr<hier::Patch>patch=*p;

        d_patch_strategy->ideal_gas_knl(*patch,true);
    }
}

void LagrangianEulerianLevelIntegrator::halfStepHaloExchange(
        const boost::shared_ptr<hier::PatchLevel>& level,
        const boost::shared_ptr<hier::PatchHierarchy>& hierarchy,
        const double current_time)
{
    d_patch_strategy->setExchangeFlag(HALF_STEP_EXCH);

    level->allocatePatchData(d_var_scratch_data, current_time);
    level->allocatePatchData(d_var_scratch_new_data, current_time);

    t_half_step_exchange_fill->start();
    d_half_step_schedules[level->getLevelNumber()]->fillData(current_time);
    t_half_step_exchange_fill->stop();

    level->deallocatePatchData(d_var_scratch_data);
    level->deallocatePatchData(d_var_scratch_new_data);
}

void LagrangianEulerianLevelIntegrator::lagrangianCorrector(
        const boost::shared_ptr<hier::PatchLevel>& level,
        const double dt)
{   
    revert(level);

    for (hier::PatchLevel::iterator p(level->begin()); p != level->end(); ++p) {

        boost::shared_ptr<hier::Patch>patch=*p;

        d_patch_strategy->accelerate(*patch,dt);
    }

    for (hier::PatchLevel::iterator p(level->begin()); p != level->end(); ++p) {

        boost::shared_ptr<hier::Patch>patch=*p;

        d_patch_strategy->pdv_knl(*patch,dt, false);
    }

    for (hier::PatchLevel::iterator p(level->begin()); p != level->end(); ++p) {

        boost::shared_ptr<hier::Patch>patch=*p;

        d_patch_strategy->flux_calc_knl(*patch,dt);
    }
}

void LagrangianEulerianLevelIntegrator::preCellHaloExchange(
        const boost::shared_ptr<hier::PatchLevel>& level,
        const boost::shared_ptr<hier::PatchHierarchy>& hierarchy,
        const double current_time)
{
    d_patch_strategy->setExchangeFlag(PRE_SWEEP_1_CELL_EXCH);

    level->allocatePatchData(d_var_scratch_data, current_time);
    level->allocatePatchData(d_var_scratch_new_data, current_time);

    t_pre_sweep1_cell_exchange_fill->start();
    d_pre_sweep1_cell_schedules[level->getLevelNumber()]->fillData(current_time);
    t_pre_sweep1_cell_exchange_fill->stop();

    level->deallocatePatchData(d_var_scratch_data);
    level->deallocatePatchData(d_var_scratch_new_data);
}

void LagrangianEulerianLevelIntegrator::advecCellSweep1(
        const boost::shared_ptr<hier::PatchLevel>& level)
{
    int sweep_number=1;
    LagrangianEulerianPatchStrategy::ADVEC_DIR direction;

    if(advect_x)  direction = LagrangianEulerianPatchStrategy::X;
    if(!advect_x) direction = LagrangianEulerianPatchStrategy::Y;

    for(hier::PatchLevel::iterator p(level->begin()); p != level->end(); ++p){

        boost::shared_ptr<hier::Patch>patch=*p;

        d_patch_strategy->advec_cell(*patch,sweep_number,direction);
    }
}

void LagrangianEulerianLevelIntegrator::preMomSweep1HaloExchange(
        const boost::shared_ptr<hier::PatchLevel>& level,
        const boost::shared_ptr<hier::PatchHierarchy>& hierarchy,
        const double current_time)
{
    d_patch_strategy->setExchangeFlag(PRE_SWEEP_1_MOM_EXCH);

    level->allocatePatchData(d_var_scratch_data, current_time);
    level->allocatePatchData(d_var_scratch_new_data, current_time);

    t_pre_sweep1_mom_exchange_fill->start();
    d_pre_sweep1_mom_schedules[level->getLevelNumber()]->fillData(current_time);
    t_pre_sweep1_mom_exchange_fill->stop();

    level->deallocatePatchData(d_var_scratch_data);
    level->deallocatePatchData(d_var_scratch_new_data);
}


void LagrangianEulerianLevelIntegrator::advecMomSweep1(
        const boost::shared_ptr<hier::PatchLevel>& level)
{
    int sweep_number=1;
    LagrangianEulerianPatchStrategy::ADVEC_DIR direction;

    if(advect_x)  direction = LagrangianEulerianPatchStrategy::X;
    if(!advect_x) direction = LagrangianEulerianPatchStrategy::Y;

    for(hier::PatchLevel::iterator p(level->begin()); p != level->end(); ++p){

        boost::shared_ptr<hier::Patch>patch=*p;

        /*
         * advection for x and y momentum.
         */
        d_patch_strategy->advec_mom(*patch,sweep_number,direction, LagrangianEulerianPatchStrategy::X);
    }

    for(hier::PatchLevel::iterator p(level->begin()); p != level->end(); ++p){

        boost::shared_ptr<hier::Patch>patch=*p;

        /*
         * advection for x and y momentum.
         */
        d_patch_strategy->advec_mom(*patch,sweep_number,direction, LagrangianEulerianPatchStrategy::Y);
    }
}

void LagrangianEulerianLevelIntegrator::preMomSweep2HaloExchange(
        const boost::shared_ptr<hier::PatchLevel>& level,
        const boost::shared_ptr<hier::PatchHierarchy>& hierarchy,
        const double current_time)
{
    d_patch_strategy->setExchangeFlag(PRE_SWEEP_2_MOM_EXCH);

    level->allocatePatchData(d_var_scratch_data, current_time);
    level->allocatePatchData(d_var_scratch_new_data, current_time);

    t_pre_sweep2_mom_exchange_fill->start();
    d_pre_sweep2_mom_schedules[level->getLevelNumber()]->fillData(current_time);
    t_pre_sweep2_mom_exchange_fill->stop();

    level->deallocatePatchData(d_var_scratch_data);
    level->deallocatePatchData(d_var_scratch_new_data);
}

void LagrangianEulerianLevelIntegrator::advecCellSweep2(
        const boost::shared_ptr<hier::PatchLevel>& level)
{
    int sweep_number=2;
    LagrangianEulerianPatchStrategy::ADVEC_DIR direction;

    if(advect_x)  direction = LagrangianEulerianPatchStrategy::Y;
    if(!advect_x) direction = LagrangianEulerianPatchStrategy::X;

    for(hier::PatchLevel::iterator p(level->begin()); p != level->end(); ++p){

        boost::shared_ptr<hier::Patch>patch=*p;

        d_patch_strategy->advec_cell(*patch,sweep_number,direction);
    }
}

void LagrangianEulerianLevelIntegrator::advecMomSweep2(
        const boost::shared_ptr<hier::PatchLevel>& level)
{
  int sweep_number=2;
  LagrangianEulerianPatchStrategy::ADVEC_DIR direction;

  if(advect_x)  direction = LagrangianEulerianPatchStrategy::Y;
  if(!advect_x) direction = LagrangianEulerianPatchStrategy::X;

    for(hier::PatchLevel::iterator p(level->begin()); p != level->end(); ++p){

        boost::shared_ptr<hier::Patch>patch=*p;

        /*
         * advection for x and y momentum.
         */
        d_patch_strategy->advec_mom(*patch,sweep_number, direction, LagrangianEulerianPatchStrategy::X);
    }

    for(hier::PatchLevel::iterator p(level->begin()); p != level->end(); ++p){

        boost::shared_ptr<hier::Patch>patch=*p;

        /*
         * advection for x and y momentum.
         */
        d_patch_strategy->advec_mom(*patch,sweep_number, direction, LagrangianEulerianPatchStrategy::Y);
    }
}

void LagrangianEulerianLevelIntegrator::swapAdvecDir()
{
   advect_x = !advect_x;
}

void LagrangianEulerianLevelIntegrator::timestepEoS(
        const boost::shared_ptr<hier::PatchLevel>& level)
{
    for (hier::PatchLevel::iterator ip(level->begin()); ip != level->end(); ip++) {
        boost::shared_ptr<hier::Patch> patch = *ip;

        d_patch_strategy->ideal_gas_knl(*patch, false);
    }
}

void LagrangianEulerianLevelIntegrator::preLagrangeHaloExchange(
        const boost::shared_ptr<hier::PatchLevel>& level,
        const boost::shared_ptr<hier::PatchHierarchy>& hierarchy,
        const double current_time)
{
    d_patch_strategy->setExchangeFlag(PRE_LAGRANGE_EXCH);

    level->allocatePatchData(d_var_scratch_data, current_time);
    level->allocatePatchData(d_var_scratch_new_data, current_time);

    t_pre_lagrange_exchange_fill->start();
    d_pre_lagrange_schedules[level->getLevelNumber()]->fillData(current_time);
    t_pre_lagrange_exchange_fill->stop();

    level->deallocatePatchData(d_var_scratch_data);
    level->deallocatePatchData(d_var_scratch_new_data);
}

void LagrangianEulerianLevelIntegrator::primeBoundaryHaloExchange(
        const boost::shared_ptr<hier::PatchLevel>& level,
        const boost::shared_ptr<hier::PatchHierarchy>& hierarchy,
        const double current_time)
{
    d_patch_strategy->setExchangeFlag(PRIME_CELLS_EXCH);

    level->allocatePatchData(d_var_scratch_data, current_time);
    level->allocatePatchData(d_var_scratch_new_data, current_time);

    t_prime_halos_exchange_fill->start();
    d_prime_halos_schedules[level->getLevelNumber()]->fillData(current_time);
    t_prime_halos_exchange_fill->stop();

    level->deallocatePatchData(d_var_scratch_data);
    level->deallocatePatchData(d_var_scratch_new_data);
}

void LagrangianEulerianLevelIntegrator::viscosity(
        const boost::shared_ptr<hier::PatchLevel>& level)
{
    for (hier::PatchLevel::iterator ip(level->begin()); ip != level->end(); ++ip) {
        boost::shared_ptr<hier::Patch> patch = *ip;

        d_patch_strategy->viscosity_knl(*patch);
    }
}

void LagrangianEulerianLevelIntegrator::postViscosityHaloExchange(
        const boost::shared_ptr<hier::PatchLevel>& level,
        const boost::shared_ptr<hier::PatchHierarchy>& hierarchy,
        const double current_time)
{
    d_patch_strategy->setExchangeFlag(POST_VISCOSITY_EXCH);

    level->allocatePatchData(d_var_scratch_data, current_time);
    level->allocatePatchData(d_var_scratch_new_data, current_time);

    t_post_viscosity_exchange_fill->start();
    d_post_viscosity_schedules[level->getLevelNumber()]->fillData(current_time);
    t_post_viscosity_exchange_fill->stop();

    level->deallocatePatchData(d_var_scratch_data);
    level->deallocatePatchData(d_var_scratch_new_data);
}

double LagrangianEulerianLevelIntegrator::calcDt(
        const boost::shared_ptr<hier::PatchLevel>& level)
{
    double dt = tbox::MathUtilities<double>::getMax();
    double patch_dt;

#if LOOPPRINT
    tbox::pout << "LagrangianEulerianLevelIntegrator: getLevelDt: calc_dt {{{" << std::endl;
#endif
    for (hier::PatchLevel::iterator ip(level->begin()); ip != level->end(); ++ip) {
        boost::shared_ptr<hier::Patch> patch = *ip;

        patch_dt = d_patch_strategy->
            calc_dt_knl(*patch);


        dt = tbox::MathUtilities<double>::Min(dt, patch_dt);
    }

    return dt;
}

void LagrangianEulerianLevelIntegrator::stampDataTime(
        const boost::shared_ptr<hier::PatchLevel>& level,
        const double current_time)
{
    level->allocatePatchData(d_var_cur_data, current_time);
    level->allocatePatchData(d_var_new_data, current_time);
}

void LagrangianEulerianLevelIntegrator::debugLevel(
                const boost::shared_ptr<hier::PatchLevel>& level)
{
    for (hier::PatchLevel::iterator ip(level->begin()); ip != level->end(); ++ip) {
        boost::shared_ptr<hier::Patch> patch = *ip;

        d_patch_strategy->debug_knl(*patch);
    }
}

void LagrangianEulerianLevelIntegrator::initializeCallback()
{
    t_synchronize_levels_create = tbox::TimerManager::getManager()->
        getTimer("LagrangianEulerianLevelIntegrator::t_synchronize_levels_create");
    t_synchronize_levels_fill = tbox::TimerManager::getManager()->
        getTimer("LagrangianEulerianLevelIntegrator::t_synchronize_levels_fill");

    t_fill_new_levels_create = tbox::TimerManager::getManager()->
        getTimer("LagrangianEulerianLevelIntegrator::t_fill_new_levels_create");
    t_fill_new_levels_fill = tbox::TimerManager::getManager()->
        getTimer("LagrangianEulerianLevelIntegrator::t_fill_new_levels_fill");

    t_half_step_exchange_create = tbox::TimerManager::getManager()->
        getTimer("LagrangianEulerianLevelIntegrator::t_half_step_exchange_create");
    t_half_step_exchange_fill = tbox::TimerManager::getManager()->
        getTimer("LagrangianEulerianLevelIntegrator::t_half_step_exchange_fill");

    t_prime_halos_exchange_create = tbox::TimerManager::getManager()->
        getTimer("LagrangianEulerianLevelIntegrator::t_prime_halos_exchange_create");
    t_prime_halos_exchange_fill = tbox::TimerManager::getManager()->
        getTimer("LagrangianEulerianLevelIntegrator::t_prime_halos_exchange_fill");

    t_pre_lagrange_exchange_create = tbox::TimerManager::getManager()->
        getTimer("LagrangianEulerianLevelIntegrator::t_pre_lagrange_exchange_create");
    t_pre_lagrange_exchange_fill = tbox::TimerManager::getManager()->
        getTimer("LagrangianEulerianLevelIntegrator::t_pre_lagrange_exchange_fill");

    t_post_viscosity_exchange_create = tbox::TimerManager::getManager()->
        getTimer("LagrangianEulerianLevelIntegrator::t_post_viscosity_exchange_create");
    t_post_viscosity_exchange_fill = tbox::TimerManager::getManager()->
        getTimer("LagrangianEulerianLevelIntegrator::t_post_viscosity_exchange_fill");

    t_pre_sweep1_cell_exchange_create = tbox::TimerManager::getManager()->
        getTimer("LagrangianEulerianLevelIntegrator::t_pre_sweep1_cell_exchange_create");
    t_pre_sweep1_cell_exchange_fill = tbox::TimerManager::getManager()->
        getTimer("LagrangianEulerianLevelIntegrator::t_pre_sweep1_cell_exchange_fill");

    t_pre_sweep1_mom_exchange_create = tbox::TimerManager::getManager()->
        getTimer("LagrangianEulerianLevelIntegrator::t_pre_sweep1_mom_exchange_create");
    t_pre_sweep1_mom_exchange_fill = tbox::TimerManager::getManager()->
        getTimer("LagrangianEulerianLevelIntegrator::t_pre_sweep1_mom_exchange_fill");

    t_pre_sweep2_mom_exchange_create = tbox::TimerManager::getManager()->
        getTimer("LagrangianEulerianLevelIntegrator::t_pre_sweep2_mom_exchange_create");
    t_pre_sweep2_mom_exchange_fill = tbox::TimerManager::getManager()->
        getTimer("LagrangianEulerianLevelIntegrator::t_pre_sweep2_mom_exchange_fill");

    t_tag_gradient_detector_cells = tbox::TimerManager::getManager()->
        getTimer("LagrangianEulerianLevelIntegrator::t_tag_gradient_detector_cells");
}

void LagrangianEulerianLevelIntegrator::finalizeCallback()
{
    t_synchronize_levels_create.reset();
    t_synchronize_levels_fill.reset();

    t_fill_new_levels_create.reset();
    t_fill_new_levels_fill.reset();

    t_half_step_exchange_create.reset();
    t_half_step_exchange_fill.reset();

    t_prime_halos_exchange_create.reset();
    t_prime_halos_exchange_fill.reset();

    t_pre_lagrange_exchange_create.reset();
    t_pre_lagrange_exchange_fill.reset();

    t_post_viscosity_exchange_create.reset();
    t_post_viscosity_exchange_fill.reset();

    t_pre_sweep1_cell_exchange_create.reset();
    t_pre_sweep1_cell_exchange_fill.reset();

    t_pre_sweep1_mom_exchange_create.reset();
    t_pre_sweep1_mom_exchange_fill.reset();

    t_pre_sweep2_mom_exchange_create.reset();
    t_pre_sweep2_mom_exchange_fill.reset();

    t_tag_gradient_detector_cells.reset();
}
