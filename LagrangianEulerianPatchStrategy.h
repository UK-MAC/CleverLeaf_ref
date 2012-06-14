#ifndef included_LagrangianEulerianPatchStrategy
#define included_LagrangianEulerianPatchStrategy

#include "SAMRAI/hier/Patch.h"
#include "SAMRAI/hier/VariableContext.h"
#include "SAMRAI/tbox/Dimension.h"
#include "SAMRAI/tbox/Pointer.h"

using namespace SAMRAI;

class LagrangianEulerianIntegrator;

class LagrangianEulerianPatchStrategy
{
    public:
        LagrangianEulerianPatchStrategy(const tbox::Dimension& dim);

        virtual void registerModelVariables(
                LagrangianEulerianIntegrator* integrator) = 0;

        virtual void initializeDataOnPatch(
                hier::Patch&,
                double init_data_time,
                bool initial_time) = 0;

        virtual double computeStableDtOnPatch(
                hier::Patch& patch,
                const bool initial_time,
                const double dt_time) = 0;

        virtual void tagGradientDetectorCells(
                hier::Patch& patch,
                const double regrid_time,
                const bool initial_error,
                const int tag_index);

        tbox::Pointer<hier::VariableContext> getDataContext();

        void setDataContext(
                tbox::Pointer<hier::VariableContext> context);

        void clearDataContext();

        const tbox::Dimension& getDim() const;

    private:
        const tbox::Dimension d_dim;

        tbox::Pointer<hier::VariableContext> d_data_context;
};
#endif
