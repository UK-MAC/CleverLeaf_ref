#include "LagrangianEulerianPatchStrategy.h"

LagrangianEulerianPatchStrategy::LagrangianEulerianPatchStrategy(
        const tbox::Dimension& dim):
    d_dim(dim)
{
}

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

void LagrangianEulerianPatchStrategy::clearDataContext()
{
    d_data_context = NULL;
}

const tbox::Dimension& LagrangianEulerianPatchStrategy::getDim() const
{
    return d_dim;
}
