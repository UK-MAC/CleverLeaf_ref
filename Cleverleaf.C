#include "Cleverleaf.h"

#include "SAMRAI/hier/VariableDatabase.h"
#include "SAMRAI/geom/CartesianPatchGeometry.h"

#include <iostream>
#include <cmath>

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
    d_volume  = new pdat::CellVariable<double>(d_dim, "volume", 1);
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
    integrator->registerVariable(d_volume, d_nghosts, d_grid_geometry);

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

        d_visit_writer->registerPlotQuantity("Volume",
                "SCALAR",
                vardb->mapVariableAndContextToIndex(
                    d_volume, d_plot_context));

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
    if (initial_time) {
        tbox::Pointer<pdat::CellData<double> > velocity = patch.getPatchData(d_velocity, getCurrentDataContext());
        tbox::Pointer<pdat::CellData<double> > massflux = patch.getPatchData(d_massflux, getCurrentDataContext());
        tbox::Pointer<pdat::CellData<double> > volflux = patch.getPatchData(d_volflux, getCurrentDataContext());
        tbox::Pointer<pdat::CellData<double> > pressure = patch.getPatchData(d_pressure, getCurrentDataContext());
        tbox::Pointer<pdat::CellData<double> > viscosity = patch.getPatchData(d_viscosity, getCurrentDataContext());
        tbox::Pointer<pdat::CellData<double> > soundspeed = patch.getPatchData(d_soundspeed, getCurrentDataContext());
        tbox::Pointer<pdat::CellData<double> > density = patch.getPatchData(d_density, getCurrentDataContext());
        tbox::Pointer<pdat::CellData<double> > energy = patch.getPatchData(d_energy, getCurrentDataContext());
        tbox::Pointer<pdat::CellData<double> > volume = patch.getPatchData(d_volume, getCurrentDataContext());

        /*
         * Fill density and energy with some data, these are our initial conditions.
         */
        double* density_data = density->getPointer();
        double* energy_data = energy->getPointer();

        const hier::Index ifirst = patch.getBox().lower();
        const hier::Index ilast = patch.getBox().upper();

        hier::IntVector density_ghosts = density->getGhostCellWidth();

        int imin = ifirst(0) - density_ghosts(0);
        int imax = ilast(0) + density_ghosts(0);
        int jmin = ifirst(1) - density_ghosts(1);
        int jmax = ilast(1) + density_ghosts(1);

        int nx = imax - imin + 1;
        int ny = jmax - jmin + 1;

        for(int j = jmin; j <= jmax; j++) {
            for(int i = imin; i <= imax; i++) {
                int n1 = POLY2(i,j,imin,jmin,nx);

                if (((i >= imin + 2) && (i <= imax - 2)) &&
                        ((j >= jmin + 2) && ( j <= jmax - 2))) {
                    density_data[n1] = 0.1;
                    energy_data[n1] = 1.0;
                } else {
                    density_data[n1] = 1.0;
                    energy_data[n1] = 2.5;
                }
            }
        }

        /*
         * Fill in the volume array...
         */
        const tbox::Pointer<geom::CartesianPatchGeometry> pgeom = patch.getPatchGeometry();
        const double* dxs = pgeom->getDx();
        double dx = dxs[0];
        double dy = dxs[1];
        double vol = dx*dy;

        /*
         * Use the fillAll() methods to initialise other variables for now...
         */
        velocity->fillAll(0.0);
        massflux->fillAll(0.0);
        volflux->fillAll(0.0);
        viscosity->fillAll(0.0);
        soundspeed->fillAll(0.0);
        volume->fillAll(vol);
    }
}

double Cleverleaf::computeStableDtOnPatch(
        hier::Patch& patch,
        const bool initial_time,
        const double dt_time)
{

    tbox::Pointer<pdat::CellData<double> > density0 = patch.getPatchData(d_density, getCurrentDataContext());
    tbox::Pointer<pdat::CellData<double> > density1 = patch.getPatchData(d_density, getNewDataContext());


    tbox::Pointer<pdat::CellData<double> > energy0 = patch.getPatchData(d_energy, getCurrentDataContext());
    tbox::Pointer<pdat::CellData<double> > energy1 = patch.getPatchData(d_energy, getNewDataContext());

    tbox::Pointer<pdat::CellData<double> > pressure = patch.getPatchData(d_pressure, getCurrentDataContext());

    tbox::Pointer<pdat::CellData<double> > soundspeed = patch.getPatchData(d_soundspeed, getCurrentDataContext());

    const hier::Index ifirst = patch.getBox().lower();
    const hier::Index ilast = patch.getBox().upper();

    int xmin = ifirst(0); 
    int xmax = ilast(0); 
    int ymin = ifirst(1); 
    int ymax = ilast(1); 

    /*
     * ideal_gas() predictor
     */
    ideal_gas(xmin,
            xmax,
            ymin,
            ymax,
            density0->getPointer(),
            energy0->getPointer(),
            pressure->getPointer(),
            soundspeed->getPointer());

    /* 
     * update_halos pressure, energy, density, velocity0
     *
     * viscosity
     *
     * update_halos viscosity
     *
     * calc_dt()
     *
     *
     */
    return 0.04;
}

void Cleverleaf::accelerate(
        hier::Patch& patch,
        double dt)
{
}

void Cleverleaf::ideal_gas(
        int xmin,
        int xmax,
        int ymin,
        int ymax,
        double* density,
        double* energy,
        double* pressure,
        double* soundspeed)
{
    double pressurebyenergy = 0;
    double pressurebyvolume = 0;
    double sound_speed_squared = 0;
    double v = 0;

    for(int j = ymin; j <= ymax; j++) {
        for(int i = xmin; i <= xmax; i++){
            int ind = POLY2(i,j,xmin,ymin,(xmax-xmin +1));

            v=1.0/density[ind];

            pressure[ind]=(1.4-1.0)*density[ind]*energy[ind];
            pressurebyenergy=(1.4-1.0)*density[ind];
            pressurebyvolume=-density[ind]*pressure[ind];

            sound_speed_squared=v*v*(pressure[ind]*pressurebyenergy-pressurebyvolume);

            soundspeed[ind]=sqrt(sound_speed_squared);
#ifdef DEBUG
            tbox::pout << "Updating pressure[" << ind << "]=" << pressure[ind] << ", soundspeed[" << ind << "]=" << soundspeed[ind] << std::endl;
#endif
        }
    }
}

void Cleverleaf::viscosity(
        int xmin,
        int xmax,
        int ymin,
        int ymax,
        double* density,
        double* pressure,
        double* viscosity,
        double* xvel0,
        double* yvel0)
{

  double ugrad,vgrad,grad2,pgradx,pgrady,
         pgradx2,pgrady2,grad,ygrad,pgrad,
         xgrad,div,strain2,limiter;

  for(int j = ymin; j <= ymax; j++) {
      for(int i = xmin; i <= xmax; i++){

          int n1 = POLY2(i,j,imin,jmin, (xmax-xmin+1));
          int n2 = POLY2(i+1,j,imin,jmin, (xmax-xmin+1));
          int n3 = POLY2(i,j-1,imin,jmin, (xmax-xmin+1));
          int n4 = POLY2(i-1,j,imin,jmin, (xmax-xmin+1));
          int n5 = POLY2(i,j+1,imin,jmin, (xmax-xmin+1));
          int n6 = POLY2(i+1,j+1,imin,jmin, (xmax-xmin+1));
          int n7 = POLY2(i+1,j-1,imin,jmin, (xmax-xmin+1));

          ugrad=(xvel0(j+1,k)+xvel0(j+1,k+1))-(xvel0(j,k)+xvel0(j,k+1));

          vgrad=(yvel0(j,k+1)+yvel0(j+1,k+1))-(yvel0(j,k)+yvel0(j+1,k));

          div=(celldx(j)*(ugrad)+celldy(k)*(vgrad));
          strain2=0.5*(xvel0(j,k+1)+xvel0(j+1,k+1)-xvel0(j,k)-xvel0(j+1,k))/celldy(k)&
                  +0.5*(yvel0(j+1,k)+yvel0(j+1,k+1)-yvel0(j,k)-yvel0(j,k+1))/celldx(j)

              pgradx=(pressure(j+1,k)-pressure(j-1,k))/(celldx(j)+celldx(j+1))
              pgrady=(pressure(j,k+1)-pressure(j,k-1))/(celldy(k)+celldy(k+1))

              pgradx2 = pgradx**2
              pgrady2 = pgrady**2

              limiter = ((0.5*(ugrad)/celldx(j))*pgradx2+(0.5*(vgrad)/celldy(k))*pgrady2+strain2*pgradx*pgrady)  &
              /MAX(pgradx2+pgrady2,1.0e-16_8)

              pgradx = SIGN(MAX(1.0e-16_8,ABS(pgradx)),pgradx)
              pgrady = SIGN(MAX(1.0e-16_8,ABS(pgrady)),pgrady)
              pgrad = SQRT(pgradx**2+pgrady**2)
              xgrad = ABS(celldx(j)*pgrad/pgradx)
              ygrad = ABS(celldy(k)*pgrad/pgrady)
              grad  = MIN(xgrad,ygrad)
              grad2 = grad*grad

              IF (.NOT.((limiter.GT.0.0).OR.(div.GE.0.0)))THEN
              viscosity(j,k)=2.0_8*density0(j,k)*grad2*limiter**2
              ELSE
              viscosity(j,k) = 0.0
              ENDIF

      }
  }
}
