#ifndef CLEVERLEAF_H 
#define CLEVERLEAF_H 

#include "LagrangianEulerianPatchStrategy.h"

#include "SAMRAI/tbox/Pointer.h"
#include "SAMRAI/appu/VisItDataWriter.h"
#include "SAMRAI/mesh/StandardTagAndInitStrategy.h"
#include "SAMRAI/hier/PatchHierarchy.h"
#include "SAMRAI/hier/PatchLevel.h"
#include "SAMRAI/pdat/CellVariable.h"

using namespace SAMRAI;

class Cleverleaf:
    public LagrangianEulerianPatchStrategy
{
    public:
        Cleverleaf(tbox::Pointer<hier::PatchHierarchy>);

        void initializeLevelData(
                const tbox::Pointer<hier::PatchHierarchy> hierarchy,
                const int level_number,
                const double init_data_time,
                const bool can_be_refined,
                const bool initial_time,
                const tbox::Pointer<hier::PatchLevel> old_level,
                const bool allocate_data);

        void resetHierarchyConfiguration(
                tbox::Pointer<hier::PatchHierarchy> new_hierarchy,
                int coarsest_level,
                int finest_level);

        void registerVisItDataWriter(tbox::Pointer<appu::VisItDataWriter>);

        void registerModelVariables();
    private:
        tbox::Pointer<hier::PatchHierarchy> d_hierarchy;
        tbox::Pointer<appu::VisItDataWriter> d_visit_writer;

        /*
         * Variables
         */
        pdat::CellVariable<double> d_velocity0;
        pdat::CellVariable<double> d_massflux;
        pdat::CellVariable<double> d_volflux;
        pdat::CellVariable<double> d_pressure;
        pdat::CellVariable<double> d_viscosity;
        pdat::CellVariable<double> d_soundspeed;
        pdat::CellVariable<double> d_density;
        pdat::CellVariable<double> d_energy;

        /*
         * Variable contexts
         */
        
};
#endif
