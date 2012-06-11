#include "Clamour.h"

void Clamour::initializeLevelData(
        const tbox::Pointer<hier::PatchHierarchy> hierarchy,
        const int level_number,
        const double init_data_time,
        const bool can_be_refined,
        const bool initial_time,
        const tbox::Pointer<hier::PatchLevel> old_level,
        const bool allocate_data) {}

void Clamour::resetHierarchyConfiguration(
        tbox::Pointer<hier::PatchHierarchy> new_hierarchy,
        int coarsest_level,
        int finest_level){}

void Clamour::registerVisItDataWriter(tbox::Pointer<appu::VisItDataWriter> writer) {
    d_visit_writer = writer;
}
