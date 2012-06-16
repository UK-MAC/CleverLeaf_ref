#include "Cleverleaf.h"

#include "SAMRAI/hier/VariableDatabase.h"

#include <iostream>

#define POLY2(i, j, imin, jmin, nx) ((i - imin) + (j-jmin) * nx)

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

void Cleverleaf::registerModelVariables(LagrangianEulerianIntegrator* integrator) 
{

    integrator->registerVariable(d_velocity, d_nghosts, d_grid_geometry);
    integrator->registerVariable(d_massflux, d_nghosts, d_grid_geometry);
    integrator->registerVariable(d_volflux, d_nghosts, d_grid_geometry);
    integrator->registerVariable(d_pressure, d_nghosts, d_grid_geometry);
    integrator->registerVariable(d_viscosity, d_nghosts, d_grid_geometry);
    integrator->registerVariable(d_soundspeed, d_nghosts, d_grid_geometry);
    integrator->registerVariable(d_density, d_nghosts, d_grid_geometry);
    integrator->registerVariable(d_energy, d_nghosts, d_grid_geometry);

    hier::VariableDatabase* vardb = hier::VariableDatabase::getDatabase();

    d_plot_context = integrator->getPlotContext();

    if (!(d_visit_writer.isNull())) {
        /*
         * Register scalar variables with the VisIt writer.
         */
        d_visit_writer->registerPlotQuantity("Pressure",
                "SCALAR",
                vardb->mapVariableAndContextToIndex(
                    d_pressure, d_plot_context));

        d_visit_writer->registerPlotQuantity("Viscosity",
                "SCALAR",
                vardb->mapVariableAndContextToIndex(
                    d_viscosity, d_plot_context));

        d_visit_writer->registerPlotQuantity("Soundspeed",
                "SCALAR",
                vardb->mapVariableAndContextToIndex(
                    d_soundspeed, d_plot_context));

        d_visit_writer->registerPlotQuantity("Density",
                "SCALAR",
                vardb->mapVariableAndContextToIndex(
                    d_density, d_plot_context));

        d_visit_writer->registerPlotQuantity("Energy",
                "SCALAR",
                vardb->mapVariableAndContextToIndex(
                    d_energy, d_plot_context));

        /*
         * Register vectors with the VisIt writer.
         */
        d_visit_writer->registerPlotQuantity("Velocity",
                "VECTOR",
                vardb->mapVariableAndContextToIndex(
                    d_velocity, d_plot_context));

        d_visit_writer->registerPlotQuantity("Massflux",
                "VECTOR",
                vardb->mapVariableAndContextToIndex(
                    d_massflux, d_plot_context));

        d_visit_writer->registerPlotQuantity("Volflux",
                "VECTOR",
                vardb->mapVariableAndContextToIndex(
                    d_volflux, d_plot_context));
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
    tbox::Pointer<pdat::CellData<double> > massflux = patch.getPatchData(d_massflux, getDataContext());
    tbox::Pointer<pdat::CellData<double> > volflux = patch.getPatchData(d_volflux, getDataContext());
    tbox::Pointer<pdat::CellData<double> > pressure = patch.getPatchData(d_pressure, getDataContext());
    tbox::Pointer<pdat::CellData<double> > viscosity = patch.getPatchData(d_viscosity, getDataContext());
    tbox::Pointer<pdat::CellData<double> > soundspeed = patch.getPatchData(d_soundspeed, getDataContext());
    tbox::Pointer<pdat::CellData<double> > density = patch.getPatchData(d_density, getDataContext());
    tbox::Pointer<pdat::CellData<double> > energy = patch.getPatchData(d_energy, getDataContext());

    double* pressure_data = pressure->getPointer();
    double*  u = velocity->getPointer(0);
    double*  v = velocity->getPointer(1);


    const hier::Index ifirst = patch.getBox().lower();
    const hier::Index ilast = patch.getBox().upper();

    hier::IntVector pressure_ghosts = pressure->getGhostCellWidth();

    int imin = ifirst(0) - pressure_ghosts(0);
    int imax = ilast(0) + pressure_ghosts(0);
    int jmin = ifirst(1) - pressure_ghosts(1);
    int jmax = ilast(1) + pressure_ghosts(1);

    int nx = imax - imin + 1;
    int ny = jmax - jmin + 1;

    for(int j = jmin; j <= jmax; j++) {
        for(int i = imin; i <= imax; i++) {
            int n1 = POLY2(i,j,imin,jmin, nx);

            if (((i >= imin + 2) && (i <= imax - 2)) &&
                    ((j >= jmin + 2) && ( j <= jmax - 2))) {
                pressure_data[n1] = 0.1;
            } else {
                pressure_data[n1] = 1.0;
            }

            u[n1] = 1.0;
            v[n1] = 1.0;
        }
    }

    cout << "FILLING DATA" << endl;
}

double Cleverleaf::computeStableDtOnPatch(
        hier::Patch& patch,
        const bool initial_time,
        const double dt_time)
{
    return 0.04;
}

