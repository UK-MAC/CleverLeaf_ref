#include "LagrangianEulerianPatchStrategy.h"

LagrangianEulerianPatchStrategy::LagrangianEulerianPatchStrategy(
        const tbox::Dimension& dim):
    RefinePatchStrategy(),
    d_dim(dim)
{
}

void LagrangianEulerianPatchStrategy::setCurrentDataContext(
        boost::shared_ptr<hier::VariableContext> context)
{
    d_current_data_context = context;
}

void LagrangianEulerianPatchStrategy::setNewDataContext(
        boost::shared_ptr<hier::VariableContext> context)
{
    d_new_data_context = context;
}

void LagrangianEulerianPatchStrategy::setScratchDataContext(
        boost::shared_ptr<hier::VariableContext> context)
{
    d_scratch_data_context = context;
}

boost::shared_ptr<hier::VariableContext> LagrangianEulerianPatchStrategy::getCurrentDataContext()
{
    return d_current_data_context;
}

boost::shared_ptr<hier::VariableContext> LagrangianEulerianPatchStrategy::getNewDataContext()
{
    return d_new_data_context;
}

boost::shared_ptr<hier::VariableContext> LagrangianEulerianPatchStrategy::getScratchDataContext()
{
    return d_scratch_data_context;
}

const tbox::Dimension& LagrangianEulerianPatchStrategy::getDim() const
{
    return d_dim;
}

hier::IntVector LagrangianEulerianPatchStrategy::getRefineOpStencilWidth( const tbox::Dimension &dim ) const 
{
    return hier::IntVector(dim,0);
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
