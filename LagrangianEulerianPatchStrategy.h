#ifndef included_LagrangianEulerianPatchStrategy
#define included_LagrangianEulerianPatchStrategy

class LagrangianEulerianPatchStrategy
{
    public:
        LagrangianEulerianPatchStrategy();

        virtual void registerModelVariables(
                LagrangianEulerianIntegrator* integrator) = 0;

        virtual double computeStableDtOnPatch(
                hier::Patch& patch,
                const bool initial_time,
                const double dt_time);

    private:
};
#endif
