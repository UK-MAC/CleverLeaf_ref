#include "Cleverleaf.h"

#include "SAMRAI/hier/VariableDatabase.h"
#include "SAMRAI/geom/CartesianPatchGeometry.h"

#include "SAMRAI/pdat/CellData.h"
#include "SAMRAI/pdat/NodeData.h"
#include "SAMRAI/pdat/EdgeData.h"

#include <iostream>
#include <cmath>

#define POLY2(i, j, imin, jmin, nx) ((i - imin) + (j-jmin) * nx)

// Arrays are defined up here so we can access them with fortran notation
#define density(i,j) density[((i-xmin)) + (j-ymin)*nx]

#define energy(i,j) energy[((i-xmin)) + (j-ymin)*nx]

#define xarea(i,j) xarea[((i-xmin)) + (j-ymin)*nx]
#define yarea(i,j) yarea[((i-xmin)) + (j-ymin)*nx]
#define volume(i,j) volume[((i-xmin)) + (j-ymin)*nx]
#define density0(i,j) density0[((i-xmin)) + (j-ymin)*nx]
#define density1(i,j) density1[((i-xmin)) + (j-ymin)*nx]
#define energy0(i,j) energy0[((i-xmin)) + (j-ymin)*nx]
#define energy1(i,j) energy1[((i-xmin)) + (j-ymin)*nx]
#define pressure(i,j) pressure[((i-xmin)) + (j-ymin)*nx]
#define viscosity(i,j) viscosity[((i-xmin)) + (j-ymin)*nx]
#define celldx(i,j) celldx[((i-xmin)) + (j-ymin)*nx]
#define celldy(i,j) celldy[((i-xmin)) + (j-ymin)*nx]
#define soundspeed(i,j) soundspeed[((i-xmin)) + (j-ymin)*nx]

#define xvel0(j,k) xvel0[((j-xmin)) + (k-ymin)*(nx+1)]
#define yvel0(j,k) yvel0[((j-xmin)) + (k-ymin)*(nx+1)]
#define xvel1(j,k) xvel1[((j-xmin)) + (k-ymin)*(nx+1)]
#define yvel1(j,k) yvel1[((j-xmin)) + (k-ymin)*(nx+1)]
#define stepbymass(j,k) stepbymass[((j)) + (k)*(nx+1)]

#define vol_flux_x(j,k) vol_flux_x[((j-xmin)) + (k-ymin)*(nx+1)]
#define vol_flux_y(j,k) vol_flux_y[((j-xmin)) + (k-ymin)*(nx+1)]

#define volume_change(i,j) volume_change[((i)) + (j)*nx]

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
    d_velocity = new pdat::NodeVariable<double>(d_dim, "velocity", d_dim.getValue());

    d_massflux  = new pdat::EdgeVariable<double>(d_dim, "massflux", 1);
    d_volflux   = new pdat::EdgeVariable<double>(d_dim, "volflux", 1);

    d_pressure  = new pdat::CellVariable<double>(d_dim, "pressure", 1);
    d_viscosity  = new pdat::CellVariable<double>(d_dim, "viscosity", 1);
    d_soundspeed  = new pdat::CellVariable<double>(d_dim, "soundspeed", 1);
    d_density  = new pdat::CellVariable<double>(d_dim, "density", 1);
    d_energy  = new pdat::CellVariable<double>(d_dim, "energy", 1);
    d_volume  = new pdat::CellVariable<double>(d_dim, "volume", 1);

    d_celldeltas = new pdat::CellVariable<double>(d_dim, "celldelta", d_dim.getValue());
    d_cellcoords = new pdat::CellVariable<double>(d_dim, "cellcoords", d_dim.getValue());

    d_vertexdeltas = new pdat::NodeVariable<double>(d_dim, "vertexdeltas", d_dim.getValue());
    d_vertexcoords = new pdat::NodeVariable<double>(d_dim, "vertexcoords", d_dim.getValue());
}

void Cleverleaf::registerModelVariables(
        LagrangianEulerianIntegrator* integrator) 
{
    integrator->registerVariable(d_velocity, LagrangianEulerianIntegrator::FIELD,  d_nghosts, d_grid_geometry);
    integrator->registerVariable(d_massflux, LagrangianEulerianIntegrator::NORMAL,  d_nghosts, d_grid_geometry);
    integrator->registerVariable(d_volflux, LagrangianEulerianIntegrator::NORMAL,  d_nghosts, d_grid_geometry);

    integrator->registerVariable(d_pressure, LagrangianEulerianIntegrator::NORMAL,  d_nghosts, d_grid_geometry);
    integrator->registerVariable(d_viscosity, LagrangianEulerianIntegrator::NORMAL,  d_nghosts, d_grid_geometry);
    integrator->registerVariable(d_soundspeed, LagrangianEulerianIntegrator::NORMAL,  d_nghosts, d_grid_geometry);
    integrator->registerVariable(d_density, LagrangianEulerianIntegrator::FIELD,  d_nghosts, d_grid_geometry);
    integrator->registerVariable(d_energy, LagrangianEulerianIntegrator::FIELD,  d_nghosts, d_grid_geometry);
    integrator->registerVariable(d_volume, LagrangianEulerianIntegrator::NORMAL,  d_nghosts, d_grid_geometry);

    integrator->registerVariable(d_celldeltas, LagrangianEulerianIntegrator::NORMAL,  d_nghosts, d_grid_geometry);
    integrator->registerVariable(d_cellcoords, LagrangianEulerianIntegrator::NORMAL,  d_nghosts, d_grid_geometry);

    integrator->registerVariable(d_vertexdeltas, LagrangianEulerianIntegrator::NORMAL,  d_nghosts, d_grid_geometry);
    integrator->registerVariable(d_vertexcoords, LagrangianEulerianIntegrator::NORMAL,  d_nghosts, d_grid_geometry);

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

//        d_visit_writer->registerPlotQuantity("Massflux",
//                "SCALAR",
//                vardb->mapVariableAndContextToIndex(
//                    d_massflux, d_plot_context));

//        d_visit_writer->registerPlotQuantity("Volflux",
//                "SCALAR",
//                vardb->mapVariableAndContextToIndex(
//                    d_volflux, d_plot_context));

        d_visit_writer->registerPlotQuantity("Vertexcoords",
                "VECTOR",
                vardb->mapVariableAndContextToIndex(
                    d_vertexcoords, d_plot_context));

        d_visit_writer->registerPlotQuantity("Cellcoords",
                "VECTOR",
                vardb->mapVariableAndContextToIndex(
                    d_cellcoords, d_plot_context));
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
        tbox::Pointer<pdat::NodeData<double> > velocity = patch.getPatchData(d_velocity, getCurrentDataContext());
        tbox::Pointer<pdat::EdgeData<double> > massflux = patch.getPatchData(d_massflux, getCurrentDataContext());
        tbox::Pointer<pdat::EdgeData<double> > volflux = patch.getPatchData(d_volflux, getCurrentDataContext());
        tbox::Pointer<pdat::CellData<double> > pressure = patch.getPatchData(d_pressure, getCurrentDataContext());
        tbox::Pointer<pdat::CellData<double> > viscosity = patch.getPatchData(d_viscosity, getCurrentDataContext());
        tbox::Pointer<pdat::CellData<double> > soundspeed = patch.getPatchData(d_soundspeed, getCurrentDataContext());
        tbox::Pointer<pdat::CellData<double> > v_density = patch.getPatchData(d_density, getCurrentDataContext());
        tbox::Pointer<pdat::CellData<double> > v_energy = patch.getPatchData(d_energy, getCurrentDataContext());
        tbox::Pointer<pdat::CellData<double> > volume = patch.getPatchData(d_volume, getCurrentDataContext());

        tbox::Pointer<pdat::CellData<double> > celldeltas = patch.getPatchData(
                d_celldeltas,
                getCurrentDataContext());
        tbox::Pointer<pdat::CellData<double> > cellcoords = patch.getPatchData(
                d_cellcoords,
                getCurrentDataContext());

        tbox::Pointer<pdat::NodeData<double> > vertexdeltas = patch.getPatchData(
                d_vertexdeltas,
                getCurrentDataContext());
        tbox::Pointer<pdat::NodeData<double> > vertexcoords = patch.getPatchData(
                d_vertexcoords,
                getCurrentDataContext());

        /*
         * Fill density and energy with some data, these are our initial conditions.
         */
        double* density = v_density->getPointer();
        double* energy = v_energy->getPointer();

        const hier::Index ifirst = patch.getBox().lower();
        const hier::Index ilast = patch.getBox().upper();

        hier::IntVector density_ghosts = v_density->getGhostCellWidth();

        int xmin = ifirst(0) - density_ghosts(0);
        int xmax = ilast(0) + density_ghosts(0);
        int ymin = ifirst(1) - density_ghosts(1);
        int ymax = ilast(1) + density_ghosts(1);

        int nx = xmax - xmin + 1;
        int ny = ymax - ymin + 1;

        for(int j = ymin; j <= ymax; j++) {
            for(int i = xmin; i <= xmax; i++) {
                //int n1 = POLY2(i,j,imin,jmin,nx);

                if (((i >= xmin + 4) && (i <= xmax - 4)) &&
                        ((j >= ymin + 4) && ( j <= ymax - 4))) {
                    density(i,j) = 1.0;
                    energy(i,j) = 2.5;
                } else {
                    density(i,j) = 0.1;
                    energy(i,j) = 1.0;
                }
            }
        }

        /*
         * Fill in the volume array...
         */
        const tbox::Pointer<geom::CartesianPatchGeometry> pgeom = patch.getPatchGeometry();
        const double* dxs = pgeom->getDx();
        const double* coords = pgeom->getXLower();
        double rxmin = coords[0];
        double rymin = coords[1];
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
        pressure->fillAll(0.0);

        /*
         * Fill in arrays of dx/dy
         */
        celldeltas->fill(dx, 0);
        celldeltas->fill(dy, 1);

        vertexdeltas->fill(dx, 0);
        vertexdeltas->fill(dy, 1);
        
        double* vertexx = vertexcoords->getPointer(0);
        double* vertexy = vertexcoords->getPointer(1);

        int xcount = 0;
        int ycount = 0;

        hier::IntVector vertex_ghosts = vertexcoords->getGhostCellWidth();

        int vimin = ifirst(0) - vertex_ghosts(0);
        int vimax = ilast(0) + vertex_ghosts(0);
        int vjmin = ifirst(1) - vertex_ghosts(1);
        int vjmax = ilast(1) + vertex_ghosts(1);

        vimax+=1;
        vjmax+=1;

        int vnx = vimax - vimin + 1;

        for(int j = vjmin; j <= vjmax; j++) {

            xcount = 0;

            for(int i = vimin; i <= vimax; i++) {
                int ind = POLY2(i,j,vimin,vjmin,vnx);

                vertexx[ind] = rxmin + dx*(xcount-rxmin);
                vertexy[ind] = rymin + dy*(ycount-rymin);

                xcount++;
            }

            ycount++;
        }

        double* cellx = cellcoords->getPointer(0);
        double* celly = cellcoords->getPointer(1);

        for(int j = ymin; j <= ymax; j++) {
            for(int i = xmin; i <= xmax; i++) {
                int vind = POLY2(i,j,vimin,vjmin,vnx);
                int vind2 = POLY2(i+1,j,vimin,vjmin,vnx);
                int vind3 = POLY2(i,j+1,vimin,vjmin,vnx);

                int ind = POLY2(i,j,xmin,ymin,nx);

                cellx[ind] = 0.5*(vertexx[vind]+vertexx[vind2]);
                celly[ind] = 0.5*(vertexy[vind]+vertexy[vind3]);
            }
        }
    }
}

void Cleverleaf::accelerate(
        hier::Patch& patch,
        double dt)
{
    double nodal_mass;

    const hier::Index ifirst = patch.getBox().lower();
    const hier::Index ilast = patch.getBox().upper();

    int xmin = ifirst(0); 
    int xmax = ilast(0); 
    int ymin = ifirst(1); 
    int ymax = ilast(1); 

    int nx = xmax - xmin + 1;

    tbox::Pointer<pdat::CellData<double> > v_density = patch.getPatchData(d_density, getCurrentDataContext());
    tbox::Pointer<pdat::CellData<double> > v_volume = patch.getPatchData(d_volume, getCurrentDataContext());
    tbox::Pointer<pdat::CellData<double> > v_celldeltas = patch.getPatchData(d_celldeltas, getCurrentDataContext());
    tbox::Pointer<pdat::CellData<double> > v_pressure = patch.getPatchData(d_pressure, getCurrentDataContext());
    tbox::Pointer<pdat::CellData<double> > v_viscosity = patch.getPatchData(d_viscosity, getCurrentDataContext());
    tbox::Pointer<pdat::NodeData<double> > v_vel0 = patch.getPatchData(d_velocity, getCurrentDataContext());
    tbox::Pointer<pdat::NodeData<double> > v_vel1 = patch.getPatchData(d_velocity, getNewDataContext());

    double* xarea = v_celldeltas->getPointer(1);
    double* yarea = v_celldeltas->getPointer(0);
    double* volume = v_volume->getPointer();
    double* density0 = v_density->getPointer();
    double* pressure = v_pressure->getPointer();
    double* viscosity = v_viscosity->getPointer();
    double* xvel0 = v_vel0->getPointer(0);
    double* yvel0 = v_vel0->getPointer(1);
    double* xvel1 = v_vel1->getPointer(0);
    double* yvel1 = v_vel1->getPointer(1);

    double* stepbymass;
    stepbymass = (double*) malloc(((xmax-xmin+2)*(ymax-ymin+2))*sizeof(double));

    for(int k = ymin; k <= ymax+1; k++) {
        for(int j = xmin; j <= xmax+1; j++ ) {

            nodal_mass=(density0(j-1,k-1)*volume(j-1,k-1)
                    +density0(j,k-1)*volume(j,k-1)
                    +density0(j,k)*volume(j,k)
                    +density0(j-1,k)*volume(j-1,k))
                    *0.25;

                stepbymass(j,k)=0.5*dt/nodal_mass;
        }
    }

    for(int k = ymin; k <= ymax+1; k++) {
        for(int j = xmin; j <= xmax+1; j++ ) {

            xvel1(j,k)=xvel0(j,k)-stepbymass(j,k)*(xarea(j,k)*(pressure(j,k)-pressure(j-1,k))
                    +xarea(j,k-1)*(pressure(j,k-1)-pressure(j-1,k-1)));
        }
    }

    for(int k = ymin; k <= ymax+1; k++) {
        for(int j = xmin; j <= xmax+1; j++ ) {

            yvel1(j,k)=yvel0(j,k)-stepbymass(j,k)*(yarea(j,k)*(pressure(j,k)-pressure(j,k-1))
                    +yarea(j-1,k)*(pressure(j-1,k)-pressure(j-1,k-1)));
        }
    }

    for(int k = ymin; k <= ymax+1; k++) {
        for(int j = xmin; j <= xmax+1; j++ ) {

            xvel1(j,k)=xvel1(j,k)-stepbymass(j,k)*(xarea(j,k)*(viscosity(j,k)-viscosity(j-1,k))
                    +xarea(j,k-1)*(viscosity(j,k-1)-viscosity(j-1,k-1)));

        }
    }

    for(int k = ymin; k <= ymax+1; k++) {
        for(int j = xmin; j <= xmax+1; j++ ) {

            yvel1(j,k)=yvel1(j,k)-stepbymass(j,k)*(yarea(j,k)*(viscosity(j,k)-viscosity(j,k-1))
                    +yarea(j-1,k)*(viscosity(j-1,k)-viscosity(j-1,k-1)));
        }
    }

    delete stepbymass;
}

void Cleverleaf::ideal_gas_knl(
        hier::Patch& patch,
        bool predict)
{

    const hier::Index ifirst = patch.getBox().lower();
    const hier::Index ilast = patch.getBox().upper();

    int xmin = ifirst(0); 
    int xmax = ilast(0); 
    int ymin = ifirst(1); 
    int ymax = ilast(1); 

    int nx = xmax - xmin + 1;

    tbox::Pointer<pdat::CellData<double> > v_density;
    tbox::Pointer<pdat::CellData<double> > v_energy;

    if (predict) {
        v_density = patch.getPatchData(d_density, getNewDataContext());
    } else {
        v_density = patch.getPatchData(d_density, getCurrentDataContext());
    }

    if (predict) {
        v_energy = patch.getPatchData(d_energy, getNewDataContext());
    } else {
        v_energy = patch.getPatchData(d_energy, getCurrentDataContext());
    }

    tbox::Pointer<pdat::CellData<double> > v_pressure = patch.getPatchData(d_pressure, getCurrentDataContext());


    tbox::Pointer<pdat::CellData<double> > v_soundspeed = patch.getPatchData(d_soundspeed, getCurrentDataContext());

    double* density = v_density->getPointer();
    double* energy = v_energy->getPointer();
    double* pressure = v_pressure->getPointer();
    double* soundspeed = v_soundspeed->getPointer();

    double pressurebyenergy = 0;
    double pressurebyvolume = 0;
    double sound_speed_squared = 0;
    double v = 0;

    for(int k = ymin; k <= ymax; k++) {
        for(int j = xmin; j <= xmax; j++){

            v=1.0/density(j,k);

            tbox::pout << "density: " << density(j,k) << std::endl;
            tbox::pout << "energy: " << energy(j,k) << std::endl;

            pressure(j,k)=(1.4-1.0)*density(j,k)*energy(j,k);

            pressurebyenergy=(1.4-1.0)*density(j,k);

            pressurebyvolume=-density(j,k)*pressure(j,k);

            sound_speed_squared=v*v*(pressure(j,k)*pressurebyenergy-pressurebyvolume);

            soundspeed(j,k)=sqrt(sound_speed_squared);

            tbox::pout << "Updating pressure[" << j << "][" << k << "]=" << pressure(j,k) << ", soundspeed[" << j << "][" << k << "]=" << soundspeed(j,k) << std::endl;
        }
    }
}

void Cleverleaf::viscosity_knl(
        hier::Patch& patch)
{

  double ugrad,vgrad,grad2,pgradx,pgrady,
         pgradx2,pgrady2,grad,ygrad,pgrad,
         xgrad,div,strain2,limiter;

    const hier::Index ifirst = patch.getBox().lower();
    const hier::Index ilast = patch.getBox().upper();

    int xmin = ifirst(0); 
    int xmax = ilast(0); 
    int ymin = ifirst(1); 
    int ymax = ilast(1); 

    int nx = xmax - xmin + 1;

    tbox::Pointer<pdat::CellData<double> > v_density0 = patch.getPatchData(d_density, getCurrentDataContext());

    tbox::Pointer<pdat::CellData<double> > v_pressure = patch.getPatchData(d_pressure, getCurrentDataContext());

    tbox::Pointer<pdat::CellData<double> > v_viscosity = patch.getPatchData(d_viscosity, getCurrentDataContext());
    tbox::Pointer<pdat::NodeData<double> > v_vel0 = patch.getPatchData(d_velocity, getCurrentDataContext());
    tbox::Pointer<pdat::CellData<double> > v_celldeltas = patch.getPatchData(d_celldeltas, getCurrentDataContext());

    double* density = v_density0->getPointer();
    double* pressure = v_pressure->getPointer();
    double* viscosity = v_viscosity->getPointer();
    double* xvel0 = v_vel0->getPointer(0);
    double* yvel0 = v_vel0->getPointer(1);
    double* celldx = v_celldeltas->getPointer(0);
    double* celldy = v_celldeltas->getPointer(1);

    for(int k = ymin; k <= ymax; k++) {
        for(int j = xmin; j <= xmax; j++){

            ugrad=(xvel0(j+1,k)+xvel0(j+1,k+1))-(xvel0(j,k)+xvel0(j,k+1));

#ifdef DEBUG
            tbox::pout << "xvel(j+1,k) = " << xvel0(j+1,k) << std::endl;
#endif
            vgrad=(yvel0(j,k+1)+yvel0(j+1,k+1))-(yvel0(j,k)+yvel0(j+1,k));

            div=(celldx(j,k)*(ugrad)+celldy(j,k)*(vgrad));

#ifdef DEBUG
            tbox::pout << "ugrad = " << ugrad << ",vgrad = " << vgrad << std::endl;
#endif

            strain2=0.5*(xvel0(j,k+1)+xvel0(j+1,k+1)-xvel0(j,k)-xvel0(j+1,k))/celldy(j,k)
                +0.5*(yvel0(j+1,k)+yvel0(j+1,k+1)-yvel0(j,k)-yvel0(j,k+1))/celldx(j,k);

            pgradx=(pressure(j+1,k)-pressure(j-1,k))/(celldx(j,k)+celldx(j+1,k));
            pgrady=(pressure(j,k+1)-pressure(j,k-1))/(celldy(j,k)+celldy(j,k+1));

            pgradx2=pgradx*pgradx;
            pgrady2=pgrady*pgrady;

            limiter=((0.5*(ugrad)/celldx(j,k))*pgradx2+(0.5*(vgrad)/celldy(j,k))*pgrady2+strain2*pgradx*pgrady)/max(pgradx2+pgrady2,1.0e-16);

            pgradx = copysign(max(1.0e-16,abs(pgradx)),pgradx);
            pgrady = copysign(max(1.0e-16,abs(pgrady)),pgrady);
            pgrad = sqrt((pgradx*pgradx)+(pgrady*pgrady));
            xgrad = abs(celldx(j,k)*pgrad/pgradx);
            ygrad = abs(celldy(j,k)*pgrad/pgrady);
            grad  = min(xgrad,ygrad);
            grad2 = grad*grad;

#ifdef DEBUG
            tbox::pout << "limiter = " << limiter << ",div = " << div << std::endl;
#endif
            if(!((limiter > 0.0) || (div >= 0.0))) {
#ifdef DEBUG
                tbox::pout << "Updating viscosity[" << j-xmin << "][" << k-ymin << "]=" << 2.0*density(j,k)*grad2*(limiter*limiter) << std::endl;
#endif
                viscosity(j,k)=2.0*density(j,k)*grad2*(limiter*limiter);
            } else {
#ifdef DEBUG
                tbox::pout << "[ZERO] Updating viscosity[" << j-xmin << "][" << k-ymin << "]=" << 0.0 << std::endl;
#endif
                viscosity(j,k) = 0.0;
            }
        }
  }
}

double Cleverleaf::calc_dt_knl(
        hier::Patch& patch)
{
    double div,dsx,dsy,dtut,dtvt,dtct,dtdivt,cc,dv1,dv2,kldt,jldt;
    double xl_pos,
           yl_pos;

    const hier::Index ifirst = patch.getBox().lower();
    const hier::Index ilast = patch.getBox().upper();

    int xmin = ifirst(0); 
    int xmax = ilast(0); 
    int ymin = ifirst(1); 
    int ymax = ilast(1); 

    tbox::Pointer<pdat::CellData<double> > v_density0 = patch.getPatchData(d_density, getCurrentDataContext());
    tbox::Pointer<pdat::CellData<double> > v_density1 = patch.getPatchData(d_density, getNewDataContext());


    tbox::Pointer<pdat::CellData<double> > v_energy0 = patch.getPatchData(d_energy, getCurrentDataContext());
    tbox::Pointer<pdat::CellData<double> > v_energy1 = patch.getPatchData(d_energy, getNewDataContext());

    tbox::Pointer<pdat::CellData<double> > v_pressure = patch.getPatchData(d_pressure, getCurrentDataContext());

    tbox::Pointer<pdat::CellData<double> > v_soundspeed = patch.getPatchData(d_soundspeed, getCurrentDataContext());

    tbox::Pointer<pdat::CellData<double> > v_viscosity = patch.getPatchData(d_viscosity, getCurrentDataContext());
    tbox::Pointer<pdat::NodeData<double> > v_vel0 = patch.getPatchData(d_velocity, getCurrentDataContext());
    tbox::Pointer<pdat::CellData<double> > v_celldeltas = patch.getPatchData(d_celldeltas, getCurrentDataContext());

    tbox::Pointer<pdat::CellData<double> > v_volume = patch.getPatchData(d_volume, getCurrentDataContext());
    tbox::Pointer<pdat::CellData<double> > v_cellcoords = patch.getPatchData(d_cellcoords, getCurrentDataContext());
    tbox::Pointer<pdat::CellData<double> > v_energy = patch.getPatchData(d_energy, getCurrentDataContext());

    double* celldx = v_celldeltas->getPointer(0);
    double* celldy = v_celldeltas->getPointer(1);
    double* soundspeed = v_soundspeed->getPointer();
    double* viscosity = v_viscosity->getPointer();
    double* pressure = v_pressure->getPointer();
    double* xvel0 = v_vel0->getPointer(0);
    double* yvel0 = v_vel0->getPointer(1);
    double* density0 = v_density0->getPointer();
    double* energy = v_energy->getPointer();
    double* xarea = v_celldeltas->getPointer(1);
    double* yarea = v_celldeltas->getPointer(0);
    double* volume = v_volume->getPointer();
    double* cellx = v_cellcoords->getPointer(0);
    double* celly = v_cellcoords->getPointer(1);

    double dt_min_val = 1.0e+21;
    double small=0;
    double dtc_safe = 0.9;
    int dtl_control;
    double dtu_safe = 0.5;
    double dtv_safe = 0.5;
    double dtdiv_safe = 0.7;

    for(int k = ymin; k <= ymax; k++) {
        for(int j = xmin; j <= xmax; j++){

            int n1 = POLY2(j,k,xmin,ymin, (xmax-xmin+1));
            int n2 = POLY2(j+1,k,xmin,ymin, (xmax-xmin+1));
            int n3 = POLY2(j,k-1,xmin,ymin, (xmax-xmin+1));
            int n4 = POLY2(j-1,k,xmin,ymin, (xmax-xmin+1));
            int n5 = POLY2(j,k+1,xmin,ymin, (xmax-xmin+1));
            int n6 = POLY2(j+1,k+1,xmin,ymin, (xmax-xmin+1));
            int n7 = POLY2(j+1,k-1,xmin,ymin, (xmax-xmin+1));

            dsx=celldx[n1];
            dsy=celldy[n1];

            cc=soundspeed[n1]*soundspeed[n1];
            cc=cc+2.0*viscosity[n1]/density0[n1];
            cc=max(sqrt(cc),1.0e-16);

            dtct=min(dsx,dsy)/cc;

            div=0.0;

            dv1=(xvel0[n1]+xvel0[n5])*xarea[n1];
            dv2=(xvel0[n2]+xvel0[n6])*xarea[n2];

            div=div+dv2-dv1;

            dtut=2.0*volume[n1]/max(abs(dv1), max(abs(dv2),1.0e-16*volume[n1]));

            dv1=(yvel0[n1]+yvel0[n2])*yarea[n1];
            dv2=(yvel0[n5]+yvel0[n6])*yarea[n5];

            div=div+dv2-dv1;

            dtvt=2.0*volume[n1]/max(abs(dv1),max(abs(dv2),1.0e-16*volume[n1]));

            div=div/(2.0*volume[n1]);

            if (div < -1.0e-16) {
                dtdivt=-1.0/div;
            } else {
                dtdivt=1.0e+21;
            }

            dtct=dtct*dtc_safe;
            if (dtct < dt_min_val) {
                jldt=j;
                kldt=k;
                dt_min_val=dtct;
                dtl_control=1;
                xl_pos=cellx[n1];
                yl_pos=celly[n1];
            }

            dtut=dtut*dtu_safe;
            if (dtut < dt_min_val) {
                jldt=j;
                kldt=k;
                dt_min_val=dtut;
                dtl_control=2;
                xl_pos=cellx[n1];
                yl_pos=celly[n1];
            }

            dtvt=dtvt*dtv_safe;
            if (dtvt < dt_min_val) {
                jldt=j;
                kldt=k;
                dt_min_val=dtvt;
                dtl_control=3;
                xl_pos=cellx[n1];
                yl_pos=celly[n1];
            }

            dtdivt=dtdivt*dtdiv_safe;
            if (dtdivt < dt_min_val) {
                jldt=j;
                kldt=k;
                dt_min_val=dtdivt;
                dtl_control=4;
                xl_pos=cellx[n1];
                yl_pos=celly[n1];
            }
        }

    }

    cout << "RETURNING " << dt_min_val << " FROM calc_dt_knl()" << endl;
    return dt_min_val;
}

void Cleverleaf::pdv_knl(
        hier::Patch& patch,
        double dt,
        bool predict)
{
    double left_flux, right_flux, bottom_flux, top_flux, total_flux;
    double recip_volume, energy_change, min_cell_volume;

    /*
     * Sort out x and y bounds of patch
     */
    const hier::Index ifirst = patch.getBox().lower();
    const hier::Index ilast = patch.getBox().upper();

    int xmin = ifirst(0); 
    int xmax = ilast(0); 

    int ymin = ifirst(1); 
    int ymax = ilast(1); 

    int nx = xmax - xmin + 1;

    /*
     * Get necessary variables
     */
    tbox::Pointer<pdat::CellData<double> > area = patch.getPatchData(d_celldeltas, getCurrentDataContext());
    tbox::Pointer<pdat::CellData<double> > vol = patch.getPatchData(d_volume, getCurrentDataContext());
    tbox::Pointer<pdat::CellData<double> > dens0 = patch.getPatchData(d_density, getCurrentDataContext());
    tbox::Pointer<pdat::CellData<double> > dens1 = patch.getPatchData(d_density, getNewDataContext());
    tbox::Pointer<pdat::CellData<double> > en0 = patch.getPatchData(d_energy, getCurrentDataContext());
    tbox::Pointer<pdat::CellData<double> > en1 = patch.getPatchData(d_energy, getNewDataContext());
    tbox::Pointer<pdat::CellData<double> > pres = patch.getPatchData(d_pressure, getCurrentDataContext());
    tbox::Pointer<pdat::CellData<double> > visc = patch.getPatchData(d_viscosity, getCurrentDataContext());
    tbox::Pointer<pdat::NodeData<double> > v0 = patch.getPatchData(d_velocity, getCurrentDataContext());
    tbox::Pointer<pdat::NodeData<double> > v1 = patch.getPatchData(d_velocity, getNewDataContext());

    double* xarea = area->getPointer(1);   
    double* yarea = area->getPointer(0);   
    double* volume = vol->getPointer(); 
    double* density0 = dens0->getPointer();
    double* density1  = dens1->getPointer();
    double* energy0 = en0->getPointer();
    double* energy1 = en1->getPointer();
    double* pressure = pres->getPointer();
    double* viscosity= visc->getPointer();
    double* xvel0 = v0->getPointer(0);
    double* xvel1 = v1->getPointer(0);
    double* yvel0 = v0->getPointer(1);
    double* yvel1 = v1->getPointer(1);

    double* volume_change = (double*) malloc(((xmax-xmin+1)*(ymax-ymin+1))*sizeof(double));

    if (predict) {

        for (int  k = ymin; k <= ymax; k++) {
            for (int j = xmin; j <= xmax; j++) {

                left_flux=(xarea(j,k)*(xvel0(j,k)+xvel0(j,k+1)
                            +xvel0(j,k)+xvel0(j,k+1)))*0.25*dt*0.5;

                right_flux=(xarea(j+1,k)*(xvel0(j+1,k)+xvel0(j+1,k+1)
                            +xvel0(j+1,k)+xvel0(j+1,k+1)))*0.25*dt*0.5;

                bottom_flux=(yarea(j,k)*(yvel0(j,k)+yvel0(j+1,k)
                            +yvel0(j,k)+yvel0(j+1,k)))*0.25*dt*0.5;

                top_flux=(yarea(j,k+1)*(yvel0(j,k+1)+yvel0(j+1,k+1)
                            +yvel0(j,k+1)+yvel0(j+1,k+1)))*0.25*dt*0.5;

                total_flux=right_flux-left_flux+top_flux-bottom_flux;

                volume_change(j,k)=volume(j,k)/(volume(j,k)+total_flux);

                min_cell_volume=min(volume(j,k)+right_flux-left_flux+top_flux-bottom_flux,
                        min(volume(j,k)+right_flux-left_flux,volume(j,k)+top_flux-bottom_flux));

                //        IF(volume_change(j,k).LE.0.0) THEN ! Perhaps take these tests out so it will vectorise
                //          error_condition=1 ! Do I need atomic for OpenMP?
                //        ENDIF
                //        IF(min_cell_volume.LE.0.0) THEN
                //          error_condition=2
                //        ENDIF

                recip_volume=1.0/volume(j,k) ;

                energy_change=(pressure(j,k)/density0(j,k)+viscosity(j,k)/density0(j,k))*total_flux*recip_volume;

                tbox::pout << "energy change = " << energy_change << ", volume change = " << volume_change(j,k) << std::endl;
                tbox::pout << "[" << j << "][" << k << "]" << std::endl;

                energy1(j,k)=energy0(j,k)-energy_change;

                density1(j,k)=density0(j,k)*volume_change(j,k);
            }
        }

    } else {

        for (int  k = ymin; k <= ymax; k++) {
            for (int j = xmin; j <= xmax; j++) {

                left_flux=(xarea(j,k)*(xvel0(j,k)+xvel0(j,k+1)
                            +xvel1(j,k)+xvel1(j,k+1)))*0.25*dt;

                right_flux=(xarea(j+1,k)*(xvel0(j+1,k)+xvel0(j+1,k+1)
                            +xvel1(j+1,k)+xvel1(j+1,k+1)))*0.25*dt;

                bottom_flux=(yarea(j,k)*(yvel0(j,k)+yvel0(j+1,k)
                            +yvel1(j,k)+yvel1(j+1,k)))*0.25*dt;

                top_flux=(yarea(j,k+1)*(yvel0(j,k+1)+yvel0(j+1,k+1)
                            +yvel1(j,k+1)+yvel1(j+1,k+1)))*0.25*dt;

                total_flux=right_flux-left_flux+top_flux-bottom_flux;

                volume_change(j,k)=volume(j,k)/(volume(j,k)+total_flux);

                min_cell_volume=min(volume(j,k)+right_flux-left_flux+top_flux-bottom_flux,
                        min(volume(j,k)+right_flux-left_flux,volume(j,k)+top_flux-bottom_flux));

                //        IF(volume_change(j,k).LE.0.0) THEN ! Perhaps take these tests out so it will vectorise
                //          error_condition=1 ! Do I need atomic for OpenMP?
                //        ENDIF
                //        IF(min_cell_volume.LE.0.0) THEN
                //          error_condition=2
                //        ENDIF

                recip_volume=1.0/volume(j,k);

                energy_change=(pressure(j,k)/density0(j,k)+viscosity(j,k)/density0(j,k))*total_flux*recip_volume;

                energy1(j,k)=energy0(j,k)-energy_change;

                density1(j,k)=density0(j,k)*volume_change(j,k);
            }
        }

    }

    delete volume_change;
}

void Cleverleaf::flux_calc_knl(
        hier::Patch& patch,
        double dt)
{

    const hier::Index ifirst = patch.getBox().lower();
    const hier::Index ilast = patch.getBox().upper();

    int xmin = ifirst(0); 
    int xmax = ilast(0); 
    int ymin = ifirst(1); 
    int ymax = ilast(1); 

    int nx = xmax - xmin + 1;

    tbox::Pointer<pdat::NodeData<double> > v_vel0 = patch.getPatchData(d_velocity, getCurrentDataContext());
    tbox::Pointer<pdat::NodeData<double> > v_vel1 = patch.getPatchData(d_velocity, getNewDataContext());

    tbox::Pointer<pdat::EdgeData<double> > v_volflux = patch.getPatchData(d_volflux, getCurrentDataContext());

    tbox::Pointer<pdat::CellData<double> > v_celldeltas = patch.getPatchData(d_celldeltas, getCurrentDataContext());

    double* xvel0 = v_vel0->getPointer(0);
    double* xvel1 = v_vel1->getPointer(0);
    double* xarea = v_celldeltas->getPointer(1);
    double* vol_flux_x = v_volflux->getPointer(1);

    double* yvel0 = v_vel0->getPointer(1);
    double* yvel1 = v_vel1->getPointer(1);
    double* yarea = v_celldeltas->getPointer(0);
    double* vol_flux_y = v_volflux->getPointer(0);

    for (int k = ymin; k <= ymax; k++) {
        for (int j = xmin; j <= xmax+1; j++) {
            vol_flux_x(j,k)=0.25*dt*xarea(j,k)
                *(xvel0(j,k)+xvel0(j,k+1)+xvel1(j,k)+xvel1(j,k+1));

            tbox::pout << "vol_flux_x(" << j << "," << k << ") = " << vol_flux_x(j,k) << std::endl;
        }
    }

    for (int k = ymin; k <= ymax+1; k++) {
        for (int j = xmin; j <= xmax; j++) {
            vol_flux_y(j,k)=0.25*dt*yarea(j,k)
                *(yvel0(j,k)+yvel0(j+1,k)+yvel1(j,k)+yvel1(j+1,k));
        }
    }
}
