#include "Cleverleaf.h"

#include "SAMRAI/hier/VariableDatabase.h"

#include <iostream>

Cleverleaf::Cleverleaf(
        tbox::Pointer<hier::PatchHierarchy> hierarchy,
        const tbox::Dimension& dim,
        tbox::Pointer<geom::CartesianGridGeometry> grid_geometry):
        LagrangianEulerianPatchStrategy(dim),
        d_dim(dim),
        d_nghosts(dim)
{

    d_hierarchy = hierarchy;
    d_grid_geometry = grid_geometry;

    d_nghosts = hier::IntVector(d_dim, 2);

    /*
     * Register variables
     */
    d_velocity = new pdat::CellVariable<double>(d_dim, "velocity", d_dim.getValue());
    d_massflux  = new pdat::CellVariable<double>(d_dim, "massflux", d_dim.getValue());
    d_volflux   = new pdat::CellVariable<double>(d_dim, "volflux", d_dim.getValue());
    d_pressure  = new pdat::CellVariable<double>(d_dim, "pressure", 1);
    d_viscosity  = new pdat::CellVariable<double>(d_dim, "viscosity", 1);
    d_soundspeed  = new pdat::CellVariable<double>(d_dim, "soundspeed", 1);
    d_density  = new pdat::CellVariable<double>(d_dim, "density", 1);
    d_energy  = new pdat::CellVariable<double>(d_dim, "energy", 1);
}

void Cleverleaf::registerModelVariables(LagrangianEulerianIntegrator* integrator) {

   integrator->registerVariable(d_density, d_nghosts, d_grid_geometry);
   integrator->registerVariable(d_velocity, d_nghosts, d_grid_geometry);

   hier::VariableDatabase* vardb = hier::VariableDatabase::getDatabase();

   d_plot_context = integrator->getPlotContext();

   if (!(d_visit_writer.isNull())) {
      d_visit_writer->registerPlotQuantity("Density",
         "SCALAR",
         vardb->mapVariableAndContextToIndex(
            d_density, d_plot_context));

      d_visit_writer->registerPlotQuantity("Velocity",
         "VECTOR",
         vardb->mapVariableAndContextToIndex(
            d_velocity, d_plot_context));
    }
}

void Cleverleaf::registerVisItDataWriter(tbox::Pointer<appu::VisItDataWriter> writer) {
    d_visit_writer = writer;
}

void Cleverleaf::initializeDataOnPatch(
        hier::Patch& patch,
        double init_data_time,
        bool initial_time)
{
        tbox::Pointer<pdat::CellData<double> > velocity = patch.getPatchData(d_velocity, getDataContext());
        tbox::Pointer<pdat::CellData<double> > density = patch.getPatchData(d_density, getDataContext());

        cout << "FILLING DATA" << endl;

        velocity->fillAll(4.4);
        density->fillAll(4.4);
}

double Cleverleaf::computeStableDtOnPatch(
        hier::Patch& patch,
        const bool initial_time,
        const double dt_time)
{
    return 0.04;
}

