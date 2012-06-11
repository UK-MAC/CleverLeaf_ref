#ifndef CLAMOUR_H_
#define CLAMOUR_H_

#include "SAMRAI/tbox/Pointer.h"
#include "SAMRAI/appu/VisItDataWriter.h"
#include "SAMRAI/mesh/StandardTagAndInitStrategy.h"
#include "SAMRAI/hier/PatchHierarchy.h"
#include "SAMRAI/hier/PatchLevel.h"

using namespace SAMRAI;

class Clamour:
    public mesh::StandardTagAndInitStrategy
{
    public:
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
    private:
        tbox::Pointer<appu::VisItDataWriter> d_visit_writer;

};
#endif
