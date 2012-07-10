#include "LagrangianEulerianPatchStrategy.h"

LagrangianEulerianPatchStrategy::LagrangianEulerianPatchStrategy(
        const tbox::Dimension& dim):
    RefinePatchStrategy(dim),
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

void LagrangianEulerianPatchStrategy::setCurrentDataContext(
        tbox::Pointer<hier::VariableContext> context)
{
    d_current_data_context = context;
}

void LagrangianEulerianPatchStrategy::setNewDataContext(
        tbox::Pointer<hier::VariableContext> context)
{
    d_new_data_context = context;
}

tbox::Pointer<hier::VariableContext> LagrangianEulerianPatchStrategy::getCurrentDataContext()
{
    return d_current_data_context;
}

tbox::Pointer<hier::VariableContext> LagrangianEulerianPatchStrategy::getNewDataContext()
{
    return d_new_data_context;;
}

const tbox::Dimension& LagrangianEulerianPatchStrategy::getDim() const
{
    return d_dim;
}
