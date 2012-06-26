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
        enum ADVEC_DIR {
            X = 1,
            Y = 2 };

        LagrangianEulerianPatchStrategy(const tbox::Dimension& dim);

        virtual void registerModelVariables(
                LagrangianEulerianIntegrator* integrator) = 0;

        virtual void initializeDataOnPatch(
                hier::Patch&,
                double init_data_time,
                bool initial_time) = 0;

        virtual void accelerate(
                hier::Patch& patch,
                double dt) = 0;

        virtual double calc_dt_knl(
                hier::Patch& patch) = 0;

        virtual void ideal_gas_knl(
                hier::Patch& patch,
                bool predict) = 0;

        virtual void viscosity_knl(
                hier::Patch& patch) = 0;

        virtual void pdv_knl(
                hier::Patch& patch,
                double dt,
                bool predict) = 0;

        virtual void flux_calc_knl(
                hier::Patch& patch,
                double dt) = 0;

        virtual void advec_cell(hier::Patch& patch,
                int sweep_number,
                ADVEC_DIR direction) = 0;

        virtual void advec_mom(hier::Patch& patch,
                int sweep_number,
                ADVEC_DIR direction) = 0;

        virtual void tagGradientDetectorCells(
                hier::Patch& patch,
                const double regrid_time,
                const bool initial_error,
                const int tag_index);

        tbox::Pointer<hier::VariableContext> getCurrentDataContext();
        tbox::Pointer<hier::VariableContext> getNewDataContext();

        void setCurrentDataContext(
                tbox::Pointer<hier::VariableContext> context);

        void setNewDataContext(
                tbox::Pointer<hier::VariableContext> context);

        const tbox::Dimension& getDim() const;

    private:
        const tbox::Dimension d_dim;

        tbox::Pointer<hier::VariableContext> d_new_data_context;

        tbox::Pointer<hier::VariableContext> d_current_data_context;
};
#endif
