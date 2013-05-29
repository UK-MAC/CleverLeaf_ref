#include "LagrangianEulerianLevelIntegrator.h"

#include "SAMRAI/tbox/MathUtilities.h"
#include "SAMRAI/hier/VariableDatabase.h"
#include "SAMRAI/hier/Patch.h"

#include <string>

#define LOOPPRINT 0

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
    d_scratch_new = hier::VariableDatabase::getDatabase()->getContext("NSCRATCH");

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
double LagrangianEulerianLevelIntegrator::getLevelDt(
        const boost::shared_ptr<hier::PatchLevel>& level,
        const double dt_time,
        const bool initial_time)
{

    const tbox::SAMRAI_MPI& mpi(level->getBoxLevel()->getMPI());

    double dt = tbox::MathUtilities<double>::getMax();
    double patch_dt;

    level->allocatePatchData(d_var_cur_data, dt_time);
    level->allocatePatchData(d_var_new_data, dt_time);

#if LOOPPRINT
    tbox::pout << "LagrangianEulerianLevelIntegrator: getLevelDt: ideal_gas corrector {{{" << std::endl;
#endif
    for (hier::PatchLevel::iterator ip(level->begin()); ip != level->end(); ip++) {
        boost::shared_ptr<hier::Patch> patch = *ip;

        d_patch_strategy->ideal_gas_knl(*patch, false);
    }
#if LOOPPRINT
    tbox::pout << "}}}" << std::endl;
#endif

    /* 
     * TODO: update_halos pressure, energy, density, velocity0
     */

    d_patch_strategy->setCurrentDataContext(d_scratch);
    level->allocatePatchData(d_var_scratch_data, dt_time);
    level->allocatePatchData(d_var_scratch_new_data, dt_time);
    

    d_bdry_fill_prime_halos->createSchedule(
                level,
                level->getLevelNumber()-1,
                d_current_hierarchy,
                d_patch_strategy)->fillData(dt_time);

    d_bdry_fill_pre_lagrange->createSchedule(level, d_patch_strategy)->fillData(dt_time);

    level->deallocatePatchData(d_var_scratch_data);
    level->deallocatePatchData(d_var_scratch_new_data);
    d_patch_strategy->setCurrentDataContext(d_current);



#if LOOPPRINT
    tbox::pout << "LagrangianEulerianLevelIntegrator: getLevelDt: viscosity {{{" << std::endl;
#endif
    for (hier::PatchLevel::iterator ip(level->begin()); ip != level->end(); ++ip) {
        boost::shared_ptr<hier::Patch> patch = *ip;

        d_patch_strategy->viscosity_knl(*patch);
    }
#if LOOPPRINT
    tbox::pout << "}}}" << std::endl;
#endif

    /*
     * TODO: update_halos viscosity
     */
    d_patch_strategy->setCurrentDataContext(d_scratch);
    level->allocatePatchData(d_var_scratch_data, dt_time);
    level->allocatePatchData(d_var_scratch_new_data, dt_time);

    //d_bdry_fill_post_viscosity->createSchedule(level, d_patch_strategy)->fillData(dt_time);

    d_bdry_fill_post_viscosity->createSchedule(
                level,
                level->getLevelNumber()-1,
                d_current_hierarchy,
                d_patch_strategy)->fillData(dt_time);

    level->deallocatePatchData(d_var_scratch_data);
    level->deallocatePatchData(d_var_scratch_new_data);
    d_patch_strategy->setCurrentDataContext(d_current);

#if LOOPPRINT
    tbox::pout << "LagrangianEulerianLevelIntegrator: getLevelDt: calc_dt {{{" << std::endl;
#endif
    for (hier::PatchLevel::iterator ip(level->begin()); ip != level->end(); ++ip) {
        boost::shared_ptr<hier::Patch> patch = *ip;

        patch_dt = d_patch_strategy->
            calc_dt_knl(*patch);


        dt = tbox::MathUtilities<double>::Min(dt, patch_dt);
    }
#if LOOPPRINT
    tbox::pout << "}}}" << std::endl;
#endif

    double global_dt = dt;

    if (mpi.getSize() > 1) {
        mpi.AllReduce(&global_dt, 1, MPI_MIN);
    }

    /*
     * TODO: Hard code the max_timestep here for now...
     */
//    if (global_dt > 0.04) {
//        return 0.04;
//    } else {
        return global_dt;
//    }
}


double LagrangianEulerianLevelIntegrator::getMaxFinerLevelDt(
        const int finer_level_number,
        const double coarse_dt,
        const hier::IntVector& ratio)
{
    return coarse_dt / double(ratio.max());
}

double LagrangianEulerianLevelIntegrator::advanceLevel(
        const boost::shared_ptr<hier::PatchLevel>& level,
        const boost::shared_ptr<hier::PatchHierarchy>& hierarchy,
        const double current_time,
        const double new_time,
        const bool first_step,
        const bool last_step,
        const bool regrid_advance)
{

    double dt = new_time - current_time;

    d_current_hierarchy = hierarchy;

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
    level->allocatePatchData(d_var_cur_data, current_time);
    level->allocatePatchData(d_var_new_data, current_time);

    /*
     * PdV kernel, predictor.
     */
#if LOOPPRINT
    tbox::pout << "LagrangianEulerianLevelIntegrator: advanceLevel: PdV predictor {{{" << std::endl;
#endif
    for (hier::PatchLevel::iterator p(level->begin()); p != level->end(); ++p) {

        boost::shared_ptr<hier::Patch>patch=*p;

        d_patch_strategy->pdv_knl(*patch,dt, true);
    }
#if LOOPPRINT
    tbox::pout << "}}}" << std::endl;
#endif

    /*
     * PdV kernel, predictor needs ideal gas call.
     */
#if LOOPPRINT
    tbox::pout << "LagrangianEulerianLevelIntegrator: advanceLevel: ideal_gas predictor {{{" << std::endl;
#endif
    for (hier::PatchLevel::iterator p(level->begin()); p != level->end(); ++p) {

        boost::shared_ptr<hier::Patch>patch=*p;

        d_patch_strategy->ideal_gas_knl(*patch,true);
    }
#if LOOPPRINT
    tbox::pout << "}}}" << std::endl;
#endif

    /*
     * TODO: Update pressure halos!
     */
    d_patch_strategy->setCurrentDataContext(d_scratch);
    level->allocatePatchData(d_var_scratch_data, current_time);
    level->allocatePatchData(d_var_scratch_new_data, current_time);
    //d_bdry_fill_half_step->createSchedule(level, d_patch_strategy)->fillData(current_time);

    d_bdry_fill_half_step->createSchedule(
                level,
                level->getLevelNumber()-1,
                hierarchy,
                d_patch_strategy)->fillData(current_time);

    //tbox::pout << "Pressure halo updated!" << std::endl;
    d_patch_strategy->setCurrentDataContext(d_current);
    

    /*
     * Call revert to reset density and energy
     */
    revert(level);

    /*
     * Acceleration due to pressure/velocity
     */ 
#if LOOPPRINT
    tbox::pout << "LagrangianEulerianLevelIntegrator: advanceLevel: acceleration {{{" << std::endl;
#endif
    for (hier::PatchLevel::iterator p(level->begin()); p != level->end(); ++p) {

        boost::shared_ptr<hier::Patch>patch=*p;

        d_patch_strategy->accelerate(*patch,dt);
    }
#if LOOPPRINT
    tbox::pout << "}}}" << std::endl;
#endif

   /*
    * PdV kernel, corrector.
    */
#if LOOPPRINT
    tbox::pout << "LagrangianEulerianLevelIntegrator: advanceLevel: PdV corrector {{{" << std::endl;
#endif
    for (hier::PatchLevel::iterator p(level->begin()); p != level->end(); ++p) {

        boost::shared_ptr<hier::Patch>patch=*p;

        d_patch_strategy->pdv_knl(*patch,dt, false);
    }
#if LOOPPRINT
    tbox::pout << "}}}" << std::endl;
#endif

   /*
    * flux calculations...
    */
    for (hier::PatchLevel::iterator p(level->begin()); p != level->end(); ++p) {

        boost::shared_ptr<hier::Patch>patch=*p;

        d_patch_strategy->flux_calc_knl(*patch,dt);
    }

   /*
    * advection here...
    */
   advection(level, hierarchy, current_time);

   advect_x = !advect_x;

    /*
     * reset_field is used to copy density, energy and velocity
     * timelevel 1 values back to timelevel 0.
     */

    //level->setTime(new_time, d_var_cur_data);

    resetField(level);

    /*
     * Compute our next dt.
     */
    return getLevelDt(level,dt,false);
}

void LagrangianEulerianLevelIntegrator::standardLevelSynchronization(
        const boost::shared_ptr<hier::PatchHierarchy>& hierarchy,
        const int coarsest_level,
        const int finest_level,
        const double sync_time,
        const tbox::Array<double>& old_times)
{
    d_current_hierarchy = hierarchy;

    double old_time = old_times[0];

    for (int fine_ln = finest_level; fine_ln > coarsest_level; fine_ln--) {
       const int coarse_ln = fine_ln - 1;

       boost::shared_ptr<hier::PatchLevel> fine_level = hierarchy->getPatchLevel(fine_ln);
       boost::shared_ptr<hier::PatchLevel> coarse_level = hierarchy->getPatchLevel(coarse_ln);

       const boost::shared_ptr<hier::PatchHierarchy> patch_hierarchy = hierarchy;

 /*
  * Sync
  */
       d_coarsen_field_data->createSchedule(
               coarse_level,
               fine_level)->coarsenData();


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

void LagrangianEulerianLevelIntegrator::resetTimeDependentData(
        const boost::shared_ptr<hier::PatchLevel>& level,
        const double new_time,
        const bool can_be_refined){}

void LagrangianEulerianLevelIntegrator::resetDataToPreadvanceState(
        const boost::shared_ptr<hier::PatchLevel>& level){}

bool LagrangianEulerianLevelIntegrator::usingRefinedTimestepping() const
{
    return false;
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


    /* 
     * Allocate storage needed to initialize level and fill data
     * from coarser levels in AMR hierarchy, potentially. Since
     * time gets set when we allocate data, re-stamp it to current
     * time if we don't need to allocate.
     */
    if (allocate_data) {
        level->allocatePatchData(d_var_cur_data, init_data_time); 
        level->allocatePatchData(d_var_new_data, init_data_time); 
    } else {
        level->setTime(init_data_time, d_var_cur_data); 
        level->setTime(init_data_time, d_var_new_data); 
    }

   level->allocatePatchData(d_var_scratch_data, init_data_time);
    level->allocatePatchData(d_var_scratch_new_data, init_data_time);

   const tbox::SAMRAI_MPI& mpi(level->getBoxLevel()->getMPI());

   if ((level_number > 0) || old_level) {

      const boost::shared_ptr<hier::PatchHierarchy> patch_hierarchy(hierarchy);

      d_patch_strategy->setCurrentDataContext(d_scratch);

      boost::shared_ptr<xfer::RefineSchedule> sched(
         d_fill_new_level->createSchedule(level,
            old_level,
            level_number-1,
            hierarchy,
            d_patch_strategy));


      sched->fillData(init_data_time);

      d_patch_strategy->setCurrentDataContext(d_current);
   }

    for (hier::PatchLevel::iterator p(level->begin()); p != level->end(); ++p) {
        boost::shared_ptr<hier::Patch> patch(*p);

        d_patch_strategy->initializeDataOnPatch(*patch,
                init_data_time,
                initial_time);
    }

//    } else {
//
    boost::shared_ptr<xfer::RefineSchedule> refine_sched;

    d_patch_strategy->setCurrentDataContext(d_scratch);

    refine_sched = d_bdry_fill_prime_halos->createSchedule(level, d_patch_strategy);
    //refine_sched->printClassData(tbox::pout);
    //tbox::pout << "Exchanged initial data" << std::endl;
    refine_sched->fillData(init_data_time);
    //level->deallocatePatchData(d_var_scratch_data);
    d_patch_strategy->setCurrentDataContext(d_current);
//   }

    if ((level_number > 0) || old_level) {

        d_patch_strategy->setCurrentDataContext(d_scratch);
        level->allocatePatchData(d_var_scratch_data, init_data_time);
        level->allocatePatchData(d_var_scratch_new_data, init_data_time);

        const boost::shared_ptr<hier::PatchHierarchy> patch_hierarchy(hierarchy);

        d_patch_strategy->setCurrentDataContext(d_scratch);

        boost::shared_ptr<xfer::RefineSchedule> sched(
                d_bdry_fill_prime_halos->createSchedule(level,
                    old_level,
                    level_number-1,
                    hierarchy,
                    d_patch_strategy));

        sched->fillData(init_data_time);

        d_patch_strategy->setCurrentDataContext(d_current);
    }

    //printFieldSummary(init_data_time, 0.04);


}

void LagrangianEulerianLevelIntegrator::resetHierarchyConfiguration (
        const boost::shared_ptr<hier::PatchHierarchy>& hierarchy,
        const int coarsest_level,
        const int finest_level){}

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

   //d_bdry_sched_advance[level_number]->fillData(error_data_time);

    for (hier::PatchLevel::iterator ip(level->begin()); ip != level->end(); ++ip) {
      boost::shared_ptr<hier::Patch> patch(*ip);
      d_patch_strategy->
      tagGradientDetectorCells(*patch,
         error_data_time,
         initial_time,
         tag_index);
   }
}

void LagrangianEulerianLevelIntegrator::registerVariable(
        const boost::shared_ptr<hier::Variable>& var,
        const int var_type,
        const int var_exchanges,
        const hier::IntVector ghosts,
        const boost::shared_ptr<hier::BaseGridGeometry>& transfer_geom)
{
    //tbox::pout << "Registering variable: " << var->getName() << std::endl;

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

        if(var->getName() == "denisty" ||
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
            d_bdry_fill_pre_sweep1_cell->registerRefine(
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

void LagrangianEulerianLevelIntegrator::advection(
        const boost::shared_ptr<hier::PatchLevel>& level,
        const boost::shared_ptr<hier::PatchHierarchy>& hierarchy,
        double current_time)
{

    int sweep_number=1;
    LagrangianEulerianPatchStrategy::ADVEC_DIR direction;

    if(advect_x)  direction = LagrangianEulerianPatchStrategy::X;
    if(!advect_x) direction = LagrangianEulerianPatchStrategy::Y;

  /*
   * TODO: update energy, density and volflux halos
   */
    //d_bdry_fill_pre_sweep1_cell->createSchedule(level, d_patch_strategy)->fillData(current_time);

    level->allocatePatchData(d_var_scratch_data, current_time);
    level->allocatePatchData(d_var_scratch_new_data, current_time);

    d_bdry_fill_pre_sweep1_cell->createSchedule(
                level,
                level->getLevelNumber()-1,
                hierarchy,
                d_patch_strategy)->fillData(current_time);

    level->deallocatePatchData(d_var_scratch_data);
    level->deallocatePatchData(d_var_scratch_new_data);

#if LOOPPRINT
    tbox::pout << "LagrangianEulerianLevelIntegrator: advection: advec_cell {{{" << std::endl;
#endif
    for(hier::PatchLevel::iterator p(level->begin()); p != level->end(); ++p){

        boost::shared_ptr<hier::Patch>patch=*p;

        d_patch_strategy->advec_cell(*patch,sweep_number,direction);
    }

#if LOOPPRINT
    tbox::pout << "}}}" << std::endl;
#endif

  /*
   * TODO: update density1, energy1, velocity1 and mass_flux halos
   */
    //d_bdry_fill_pre_sweep1_mom->createSchedule(level, d_patch_strategy)->fillData(current_time);
    level->allocatePatchData(d_var_scratch_data, current_time);
    level->allocatePatchData(d_var_scratch_new_data, current_time);

    d_bdry_fill_pre_sweep1_mom->createSchedule(
                level,
                level->getLevelNumber()-1,
                hierarchy,
                d_patch_strategy)->fillData(current_time);

    level->deallocatePatchData(d_var_scratch_data);
    level->deallocatePatchData(d_var_scratch_new_data);
#if LOOPPRINT
    tbox::pout << "LagrangianEulerianLevelIntegrator: advection: advec_mom x {{{" << std::endl;
#endif
    for(hier::PatchLevel::iterator p(level->begin()); p != level->end(); ++p){

        boost::shared_ptr<hier::Patch>patch=*p;

        /*
         * advection for x and y momentum.
         */
        d_patch_strategy->advec_mom(*patch,sweep_number,direction, LagrangianEulerianPatchStrategy::X);
    }
#if LOOPPRINT
    tbox::pout << "}}}" << std::endl;
#endif

#if LOOPPRINT
    tbox::pout << "LagrangianEulerianLevelIntegrator: advection: advec_mom y {{{" << std::endl;
#endif
    for(hier::PatchLevel::iterator p(level->begin()); p != level->end(); ++p){

        boost::shared_ptr<hier::Patch>patch=*p;

        /*
         * advection for x and y momentum.
         */
        d_patch_strategy->advec_mom(*patch,sweep_number,direction, LagrangianEulerianPatchStrategy::Y);
    }
#if LOOPPRINT
    tbox::pout << "}}}" << std::endl;
#endif

  sweep_number=2;

  if(advect_x)  direction = LagrangianEulerianPatchStrategy::Y;
  if(!advect_x) direction = LagrangianEulerianPatchStrategy::X;

#if LOOPPRINT
    tbox::pout << "LagrangianEulerianLevelIntegrator: advection sweep 2: advec_cell {{{" << std::endl;
#endif
    for(hier::PatchLevel::iterator p(level->begin()); p != level->end(); ++p){

        boost::shared_ptr<hier::Patch>patch=*p;

        d_patch_strategy->advec_cell(*patch,sweep_number,direction);
    }
#if LOOPPRINT
    tbox::pout << "}}}" << std::endl;
#endif

  /*
   * TODO: Update density1, energy1, vel1 and mass flux halos
   */
    //d_bdry_fill_pre_sweep2_mom->createSchedule(level, d_patch_strategy)->fillData(current_time);
    level->allocatePatchData(d_var_scratch_data, current_time);
    level->allocatePatchData(d_var_scratch_new_data, current_time);

    d_bdry_fill_pre_sweep2_mom->createSchedule(
                level,
                level->getLevelNumber()-1,
                hierarchy,
                d_patch_strategy)->fillData(current_time);

    level->deallocatePatchData(d_var_scratch_data);
    level->deallocatePatchData(d_var_scratch_new_data);

#if LOOPPRINT
    tbox::pout << "LagrangianEulerianLevelIntegrator: advection: advec_mom x {{{" << std::endl;
#endif
    for(hier::PatchLevel::iterator p(level->begin()); p != level->end(); ++p){

        boost::shared_ptr<hier::Patch>patch=*p;

        /*
         * advection for x and y momentum.
         */
        d_patch_strategy->advec_mom(*patch,sweep_number, direction, LagrangianEulerianPatchStrategy::X);
    }
#if LOOPPRINT
    tbox::pout << "}}}" << std::endl;
#endif

#if LOOPPRINT
    tbox::pout << "LagrangianEulerianLevelIntegrator: advection: advec_mom y {{{" << std::endl;
#endif
    for(hier::PatchLevel::iterator p(level->begin()); p != level->end(); ++p){

        boost::shared_ptr<hier::Patch>patch=*p;

        /*
         * advection for x and y momentum.
         */
        d_patch_strategy->advec_mom(*patch,sweep_number, direction, LagrangianEulerianPatchStrategy::Y);
    }
#if LOOPPRINT
    tbox::pout << "}}}" << std::endl;
#endif
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


