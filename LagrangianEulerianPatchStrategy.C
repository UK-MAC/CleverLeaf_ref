#include "LagrangianEulerianPatchStrategy.h"

void LagrangianEulerianPatchStrategy::tagGradientDetectorCells(
                hier::Patch& patch,
                const double regrid_time,
                const bool initial_error,
                const int tag_index)
{
    /*
     * Do nothing!
     */
}

tbox::Pointer<hier::VariableContext> LagrangianEulerianPatchStrategy::getDataContext()
{
    return d_data_context;
}

void LagrangianEulerianPatchStrategy::setDataContext(
        tbox::Pointer<hier::VariableContext> context)
{
    d_data_context = context;
}

const tbox::Dimension& LagrangianEulerianPatchStrategy::getDim()
{
    return d_dim;
}
