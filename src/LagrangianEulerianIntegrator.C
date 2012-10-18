#include "LagrangianEulerianIntegrator.h"

#include "SAMRAI/tbox/MathUtilities.h"
#include "SAMRAI/hier/VariableDatabase.h"
#include "SAMRAI/hier/Patch.h"

#define LOOPPRINT 0

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

    advect_x = true;

    /*
     * Communication algorithms
     */
    d_bdry_fill_half_step = new xfer::RefineAlgorithm(d_dim);
    d_bdry_fill_prime_halos = new xfer::RefineAlgorithm(d_dim);
    d_bdry_fill_pre_lagrange = new xfer::RefineAlgorithm(d_dim);
    d_bdry_fill_post_viscosity = new xfer::RefineAlgorithm(d_dim);
    d_bdry_fill_pre_sweep1_cell = new xfer::RefineAlgorithm(d_dim);
    d_bdry_fill_pre_sweep1_mom = new xfer::RefineAlgorithm(d_dim);
    d_bdry_fill_pre_sweep2_mom = new xfer::RefineAlgorithm(d_dim);

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

    /*
     * Pass these contexts up to the patch strategy
     */
    patch_strategy->setCurrentDataContext(d_current);
    patch_strategy->setNewDataContext(d_new);
    patch_strategy->setScratchDataContext(d_scratch);

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
    const tbox::SAMRAI_MPI& mpi(level->getBoxLevel()->getMPI());

    double dt = tbox::MathUtilities<double>::getMax();
    double patch_dt;

    level->allocatePatchData(d_var_cur_data, dt_time);
    //level->allocatePatchData(d_var_new_data, dt_time);

#if LOOPPRINT
    tbox::pout << "LagrangianEulerianIntegrator: getLevelDt: ideal_gas corrector {{{" << std::endl;
#endif
    for (hier::PatchLevel::Iterator ip(level); ip; ip++) {
        tbox::Pointer<hier::Patch> patch = *ip;

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
    d_bdry_fill_pre_lagrange->createSchedule(level, d_patch_strategy)->fillData(dt_time);
    level->deallocatePatchData(d_var_scratch_data);
    d_patch_strategy->setCurrentDataContext(d_current);


#if LOOPPRINT
    tbox::pout << "LagrangianEulerianIntegrator: getLevelDt: viscosity {{{" << std::endl;
#endif
    for (hier::PatchLevel::Iterator ip(level); ip; ip++) {
        tbox::Pointer<hier::Patch> patch = *ip;

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
    d_bdry_fill_post_viscosity->createSchedule(level, d_patch_strategy)->fillData(dt_time);
    level->deallocatePatchData(d_var_scratch_data);
    d_patch_strategy->setCurrentDataContext(d_current);

#if LOOPPRINT
    tbox::pout << "LagrangianEulerianIntegrator: getLevelDt: calc_dt {{{" << std::endl;
#endif
    for (hier::PatchLevel::Iterator ip(level); ip; ip++) {
        tbox::Pointer<hier::Patch> patch = *ip;

        patch_dt = d_patch_strategy->
            calc_dt_knl(*patch);

        std::cout << " patch_dt = " << patch_dt << std::endl;

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
    if (global_dt > 0.04) {
        return 0.04;
    } else {
        return global_dt;
    }
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
    level->allocatePatchData(d_var_new_data, new_time);
    level->allocatePatchData(d_var_cur_data, current_time);

    /*
     * PdV kernel, predictor.
     */
#if LOOPPRINT
    tbox::pout << "LagrangianEulerianIntegrator: advanceLevel: PdV predictor {{{" << std::endl;
#endif
    for(hier::PatchLevel::Iterator p(level);p;p++){

        tbox::Pointer<hier::Patch>patch=*p;

        d_patch_strategy->pdv_knl(*patch,dt, true);
    }
#if LOOPPRINT
    tbox::pout << "}}}" << std::endl;
#endif

    /*
     * PdV kernel, predictor needs ideal gas call.
     */
#if LOOPPRINT
    tbox::pout << "LagrangianEulerianIntegrator: advanceLevel: ideal_gas predictor {{{" << std::endl;
#endif
    for(hier::PatchLevel::Iterator p(level);p;p++){

        tbox::Pointer<hier::Patch>patch=*p;

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
    d_bdry_fill_half_step->createSchedule(level, d_patch_strategy)->fillData(current_time);
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
    tbox::pout << "LagrangianEulerianIntegrator: advanceLevel: acceleration {{{" << std::endl;
#endif
    for(hier::PatchLevel::Iterator p(level);p;p++){

        tbox::Pointer<hier::Patch>patch=*p;

        d_patch_strategy->accelerate(*patch,dt);
    }
#if LOOPPRINT
    tbox::pout << "}}}" << std::endl;
#endif

   /*
    * PdV kernel, corrector.
    */
#if LOOPPRINT
    tbox::pout << "LagrangianEulerianIntegrator: advanceLevel: PdV corrector {{{" << std::endl;
#endif
   for(hier::PatchLevel::Iterator p(level);p;p++){

        tbox::Pointer<hier::Patch>patch=*p;

        d_patch_strategy->pdv_knl(*patch,dt, false);
    }
#if LOOPPRINT
    tbox::pout << "}}}" << std::endl;
#endif

   /*
    * flux calculations...
    */
   for(hier::PatchLevel::Iterator p(level);p;p++){

        tbox::Pointer<hier::Patch>patch=*p;

        d_patch_strategy->flux_calc_knl(*patch,dt);
    }

   /*
    * advection here...
    */
   advection(level, hierarchy, new_time);

   advect_x = !advect_x;

    /*
     * reset_field is used to copy density, energy and velocity
     * timelevel 1 values back to timelevel 0.
     */

    level->setTime(new_time, d_var_cur_data);

    resetField(level);

    /*
     * Compute our next dt.
     */
    return getLevelDt(level,dt,false);
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

    d_current_hierarchy = hierarchy;

    tbox::Pointer<hier::PatchLevel> level(
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
    }

    for (hier::PatchLevel::Iterator p(level); p; p++) {
        tbox::Pointer<hier::Patch> patch(*p);

        d_patch_strategy->initializeDataOnPatch(*patch,
                init_data_time,
                initial_time);
    }

   // TODO: prime_halos_exch here. 
    d_patch_strategy->setCurrentDataContext(d_scratch);
    level->allocatePatchData(d_var_scratch_data, init_data_time);
    tbox::Pointer<xfer::RefineSchedule> refine_sched = d_bdry_fill_prime_halos->createSchedule(level, d_patch_strategy);
    //refine_sched->printClassData(tbox::pout);
    refine_sched->fillData(init_data_time);
    //tbox::pout << "Exchanged initial data" << std::endl;
    level->deallocatePatchData(d_var_scratch_data);
    d_patch_strategy->setCurrentDataContext(d_current);

    printFieldSummary(init_data_time, 0.04);


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
        const int var_type,
        const int var_exchanges,
        const hier::IntVector ghosts,
        const tbox::Pointer<hier::GridGeometry> transfer_geom)
{
    //tbox::pout << "Registering variable: " << var->getName() << std::endl;

    const tbox::Dimension dim(ghosts.getDim());

    hier::VariableDatabase* variable_db = hier::VariableDatabase::getDatabase();

    const hier::IntVector& zero_ghosts(hier::IntVector::getZero(dim));

    if((var_type & FIELD) == FIELD) {
        d_field_vars.appendItem(var);
    }

    if((var_type & REVERT) == REVERT) {
        //std::cout << "Found a revert var: " << var->getName() << std::endl;
        d_revert_vars.appendItem(var);
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

    tbox::Pointer<hier::RefineOperator> refine_op;

    if (var->getName() == "massflux") {
        //tbox::pout << "Using CONSERVATIVE_LINEAR_REFINE..." << std::endl;
        refine_op = transfer_geom->lookupRefineOperator(var, "CONSERVATIVE_LINEAR_REFINE");
    } else if (var->getName() == "volflux") {
        //tbox::pout << "Using CONSERVATIVE_LINEAR_REFINE..." << std::endl;
        refine_op = transfer_geom->lookupRefineOperator(var, "CONSERVATIVE_LINEAR_REFINE");
    } else {
        //tbox::pout << "Using LINEAR_REFINE..." << std::endl;
        refine_op = transfer_geom->lookupRefineOperator(var, "LINEAR_REFINE");
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
                    scr_id,
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
                    scr_id,
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
                    scr_id,
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
                    scr_id,
                    refine_op);
        }
    }

    d_var_cur_data.setFlag(cur_id);
    d_var_new_data.setFlag(new_id);
    d_var_scratch_data.setFlag(scr_id);
}

tbox::Pointer<hier::VariableContext> LagrangianEulerianIntegrator::getPlotContext()
{
    return d_plot_context;
}


void LagrangianEulerianIntegrator::resetField(
   const tbox::Pointer<hier::PatchLevel> level)
{
   for (hier::PatchLevel::Iterator ip(level); ip; ip++) {
      tbox::Pointer<hier::Patch> patch = *ip;

      tbox::List<tbox::Pointer<hier::Variable> >::Iterator
         field_var = d_field_vars.listStart();
      while (field_var) {

          //tbox::pout << "Copying " << field_var()->getName() << " back to tl0" << std::endl;

         tbox::Pointer<hier::PatchData> src_data =
            patch->getPatchData(field_var(), d_new);
         tbox::Pointer<hier::PatchData> dst_data =
            patch->getPatchData(field_var(), d_current);

         dst_data->copy(*src_data);
         field_var++;
      }
   }
}

void LagrangianEulerianIntegrator::revert(
   const tbox::Pointer<hier::PatchLevel> level)
{
   for (hier::PatchLevel::Iterator ip(level); ip; ip++) {
      tbox::Pointer<hier::Patch> patch = *ip;

      tbox::List<tbox::Pointer<hier::Variable> >::Iterator
         revert_var = d_revert_vars.listStart();
      while (revert_var) {
#if DEBUG
          tbox::pout << "Copying " << revert_var()->getName() << " back to tl0" << std::endl;
#endif
         tbox::Pointer<hier::PatchData> dst_data =
            patch->getPatchData(revert_var(), d_new);
         tbox::Pointer<hier::PatchData> src_data =
            patch->getPatchData(revert_var(), d_current);

        dst_data->copy(*src_data);
         revert_var++;
      }
   }
}

void LagrangianEulerianIntegrator::advection(
        const tbox::Pointer<hier::PatchLevel> level,
        const tbox::Pointer<hier::PatchHierarchy> hierarchy,
        double current_time)
{

    int sweep_number=1;
    LagrangianEulerianPatchStrategy::ADVEC_DIR direction;

    if(advect_x)  direction = LagrangianEulerianPatchStrategy::X;
    if(!advect_x) direction = LagrangianEulerianPatchStrategy::Y;

  /*
   * TODO: update energy, density and volflux halos
   */
    d_bdry_fill_pre_sweep1_cell->createSchedule(level, d_patch_strategy)->fillData(current_time);

#if LOOPPRINT
    tbox::pout << "LagrangianEulerianIntegrator: advection: advec_cell {{{" << std::endl;
#endif
    for(hier::PatchLevel::Iterator p(level);p;p++){

        tbox::Pointer<hier::Patch>patch=*p;

        d_patch_strategy->advec_cell(*patch,sweep_number,direction);
    }

#if LOOPPRINT
    tbox::pout << "}}}" << std::endl;
#endif

  /*
   * TODO: update density1, energy1, velocity1 and mass_flux halos
   */
    d_bdry_fill_pre_sweep1_mom->createSchedule(level, d_patch_strategy)->fillData(current_time);

#if LOOPPRINT
    tbox::pout << "LagrangianEulerianIntegrator: advection: advec_mom x {{{" << std::endl;
#endif
    for(hier::PatchLevel::Iterator p(level);p;p++){

        tbox::Pointer<hier::Patch>patch=*p;

        /*
         * advection for x and y momentum.
         */
        d_patch_strategy->advec_mom(*patch,sweep_number,direction, LagrangianEulerianPatchStrategy::X);
    }
#if LOOPPRINT
    tbox::pout << "}}}" << std::endl;
#endif

#if LOOPPRINT
    tbox::pout << "LagrangianEulerianIntegrator: advection: advec_mom y {{{" << std::endl;
#endif
    for(hier::PatchLevel::Iterator p(level);p;p++){

        tbox::Pointer<hier::Patch>patch=*p;

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
    tbox::pout << "LagrangianEulerianIntegrator: advection sweep 2: advec_cell {{{" << std::endl;
#endif
    for(hier::PatchLevel::Iterator p(level);p;p++){

        tbox::Pointer<hier::Patch>patch=*p;

        d_patch_strategy->advec_cell(*patch,sweep_number,direction);
    }
#if LOOPPRINT
    tbox::pout << "}}}" << std::endl;
#endif

  /*
   * TODO: Update density1, energy1, vel1 and mass flux halos
   */
    d_bdry_fill_pre_sweep2_mom->createSchedule(level, d_patch_strategy)->fillData(current_time);

#if LOOPPRINT
    tbox::pout << "LagrangianEulerianIntegrator: advection: advec_mom x {{{" << std::endl;
#endif
    for(hier::PatchLevel::Iterator p(level);p;p++){

        tbox::Pointer<hier::Patch>patch=*p;

        /*
         * advection for x and y momentum.
         */
        d_patch_strategy->advec_mom(*patch,sweep_number, direction, LagrangianEulerianPatchStrategy::X);
    }
#if LOOPPRINT
    tbox::pout << "}}}" << std::endl;
#endif

#if LOOPPRINT
    tbox::pout << "LagrangianEulerianIntegrator: advection: advec_mom y {{{" << std::endl;
#endif
    for(hier::PatchLevel::Iterator p(level);p;p++){

        tbox::Pointer<hier::Patch>patch=*p;

        /*
         * advection for x and y momentum.
         */
        d_patch_strategy->advec_mom(*patch,sweep_number, direction, LagrangianEulerianPatchStrategy::Y);
    }
#if LOOPPRINT
    tbox::pout << "}}}" << std::endl;
#endif
}

void LagrangianEulerianIntegrator::printFieldSummary(
        double time,
        int step)
{
    tbox::Pointer<hier::PatchLevel> level = d_current_hierarchy->getPatchLevel(0);
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

    for(hier::PatchLevel::Iterator p(level);p;p++){

        tbox::Pointer<hier::Patch>patch=*p;

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


