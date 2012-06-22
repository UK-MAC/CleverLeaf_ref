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

    /*
     * Pass these contexts up to the patch strategy
     */
    patch_strategy->setCurrentDataContext(d_current);
    patch_strategy->setNewDataContext(d_new);

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
    if (initial_time) {
        return 0.04;
    }

    const tbox::SAMRAI_MPI& mpi(level->getBoxLevel()->getMPI());

    double dt = tbox::MathUtilities<double>::getMax();
    double patch_dt;

    level->allocatePatchData(d_temp_var_cur_data, dt_time);
    level->allocatePatchData(d_temp_var_new_data, dt_time);

    for (hier::PatchLevel::Iterator ip(level); ip; ip++) {
        tbox::Pointer<hier::Patch> patch = *ip;

        d_patch_strategy->ideal_gas_knl(*patch, false);
    }
    
    /* 
     * TODO: update_halos pressure, energy, density, velocity0
     */

    for (hier::PatchLevel::Iterator ip(level); ip; ip++) {
        tbox::Pointer<hier::Patch> patch = *ip;

        d_patch_strategy->viscosity_knl(*patch);
    }

    /*
     * TODO: update_halos viscosity
     */

    for (hier::PatchLevel::Iterator ip(level); ip; ip++) {
        tbox::Pointer<hier::Patch> patch = *ip;

        patch_dt = d_patch_strategy->
            calc_dt_knl(*patch);

        dt = tbox::MathUtilities<double>::Min(dt, patch_dt);
    }




    double global_dt = dt;

    if (mpi.getSize() > 1) {
        mpi.AllReduce(&global_dt, 1, MPI_MIN);
    }

    return global_dt;
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
    level->allocatePatchData(d_temp_var_new_data, new_time);
    level->allocatePatchData(d_temp_var_cur_data, current_time);

    /*
     * TODO: PdV kernel.
     */
    for(hier::PatchLevel::Iterator p(level);p;p++){

        tbox::Pointer<hier::Patch>patch=*p;

        d_patch_strategy->pdv_knl(*patch,dt, true);
    }


    /*
     * reset_field is used to copy density, energy and velocity
     * timelevel 1 values back to timelevel 0.
     */

    level->setTime(new_time, d_temp_var_cur_data);

    resetField(level);


    /*
     * Compute our next dt.
     */
//    const tbox::SAMRAI_MPI& mpi(level->getBoxLevel()->getMPI());
//    double dt_next = tbox::MathUtilities<double>::getMax();
//    double patch_dt;
//
//    for (hier::PatchLevel::Iterator ip(level); ip; ip++) {
//        tbox::Pointer<hier::Patch> patch = *ip;
//
//        patch_dt = d_patch_strategy->
//            computeStableDtOnPatch(*patch,
//                    false,
//                    new_time);
//
//        dt_next = tbox::MathUtilities<double>::Min(dt_next, patch_dt);
//    }
//
//    double next_dt = dt_next;
//
//    if (mpi.getSize() > 1) {
//        mpi.AllReduce(&next_dt, 1, MPI_MIN);
//    }

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


    for (hier::PatchLevel::Iterator p(level); p; p++) {
        tbox::Pointer<hier::Patch> patch(*p);

        d_patch_strategy->initializeDataOnPatch(*patch,
                init_data_time,
                initial_time);
    }
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
        const VAR_TYPE var_type,
        const hier::IntVector ghosts,
        const tbox::Pointer<hier::GridGeometry> transfer_geom)
{
#ifdef DEBUG
    tbox::pout << "Registering variable: " << var->getName() << std::endl;
#endif

    const tbox::Dimension dim(ghosts.getDim());

    hier::VariableDatabase* variable_db = hier::VariableDatabase::getDatabase();

    const hier::IntVector& zero_ghosts(hier::IntVector::getZero(dim));

    if(var_type == FIELD)
        d_field_vars.appendItem(var);

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


void LagrangianEulerianIntegrator::resetField(
   const tbox::Pointer<hier::PatchLevel> level)
{
   for (hier::PatchLevel::Iterator ip(level); ip; ip++) {
      tbox::Pointer<hier::Patch> patch = *ip;

      tbox::List<tbox::Pointer<hier::Variable> >::Iterator
         field_var = d_field_vars.listStart();
      while (field_var) {
#ifdef DEBUG
          tbox::pout << "Copying " << field_var()->getName() << " back to tl0" << std::endl;
#endif
         tbox::Pointer<hier::PatchData> src_data =
            patch->getPatchData(field_var(), d_new);
         tbox::Pointer<hier::PatchData> dst_data =
            patch->getPatchData(field_var(), d_current);

         dst_data->copy(*src_data);
         field_var++;
      }
   }
}
