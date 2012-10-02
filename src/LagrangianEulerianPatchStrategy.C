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

void LagrangianEulerianPatchStrategy::setScratchDataContext(
        tbox::Pointer<hier::VariableContext> context)
{
    d_scratch_data_context = context;
}

tbox::Pointer<hier::VariableContext> LagrangianEulerianPatchStrategy::getCurrentDataContext()
{
    return d_current_data_context;
}

tbox::Pointer<hier::VariableContext> LagrangianEulerianPatchStrategy::getNewDataContext()
{
    return d_new_data_context;
}

tbox::Pointer<hier::VariableContext> LagrangianEulerianPatchStrategy::getScratchDataContext()
{
    return d_scratch_data_context;
}

const tbox::Dimension& LagrangianEulerianPatchStrategy::getDim() const
{
    return d_dim;
}

hier::IntVector LagrangianEulerianPatchStrategy::getRefineOpStencilWidth() const 
{
    return hier::IntVector(d_dim,0);
}

void LagrangianEulerianPatchStrategy::preprocessRefine(
        hier::Patch& fine,
        const hier::Patch& coarse,
        const hier::Box& fine_box,
        const hier::IntVector& ratio)
{
    /*
     * Dummy implementation!
     */
}

void LagrangianEulerianPatchStrategy::postprocessRefine(
        hier::Patch& fine,
        const hier::Patch& coarse,
        const hier::Box& fine_box,
        const hier::IntVector& ratio)
{
    /*
     * Dummy implementation!
     */
}
