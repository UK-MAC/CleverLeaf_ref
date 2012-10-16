#include "Cleverleaf.h"

#include "SAMRAI/hier/VariableDatabase.h"
#include "SAMRAI/geom/CartesianPatchGeometry.h"

#include "SAMRAI/pdat/CellData.h"
#include "SAMRAI/pdat/NodeData.h"
#include "SAMRAI/pdat/EdgeData.h"

#include "SAMRAI/appu/CartesianBoundaryUtilities2.h"
#include "SAMRAI/appu/CartesianBoundaryDefines.h"

#include <iostream>
#include <cmath>


#define POLY2(i, j, imin, jmin, nx) ((i - imin) + (j-jmin) * nx)

// Arrays are defined up here so we can access them with fortran notation
#define density(i,j) density[((i-xmin)) + (j-ymin)*nx]
#define energy(i,j) energy[((i-xmin)) + (j-ymin)*nx]

#define xarea(i,j) xarea[((i-xmin)) + (j-ymin)*nx]
#define yarea(i,j) yarea[((i-xmin)) + (j-ymin)*nx]
#define volume(i,j) volume[((i-xmin)) + (j-ymin)*nx]
#define density0(i,j) density0[((i)-xmin) + ((j)-ymin)*nx]
#define density1(i,j) density1[((i-xmin)) + (j-ymin)*nx]
#define energy0(i,j) energy0[((i-xmin)) + (j-ymin)*nx]
#define energy1(i,j) energy1[((i-xmin)) + (j-ymin)*nx]
#define pressure(i,j) pressure[((i-xmin)) + (j-ymin)*nx]
#define viscosity(i,j) viscosity[((i-xmin)) + (j-ymin)*nx]
#define celldx(i,j) celldx[((i-xmin)) + (j-ymin)*nx]
#define celldy(i,j) celldy[((i-xmin)) + (j-ymin)*nx]
#define soundspeed(i,j) soundspeed[((i-xmin)) + (j-ymin)*nx]
#define cellx(i,j) cellx[((i-xmin)) + (j-ymin)*nx]
#define celly(i,j) celly[((i-xmin)) + (j-ymin)*nx]

#define xvel0(j,k) xvel0[((j-xmin)) + (k-ymin)*(nx+1)]
#define yvel0(j,k) yvel0[((j-xmin)) + (k-ymin)*(nx+1)]
#define xvel1(j,k) xvel1[((j-xmin)) + (k-ymin)*(nx+1)]
#define yvel1(j,k) yvel1[((j-xmin)) + (k-ymin)*(nx+1)]
#define stepbymass(j,k) stepbymass[((j-ifirst(0))) + (k-ifirst(1))*(sbmnx+1)]

#define vertexdx(j) vertexdx[((j-xmin)) + (k-ymin)*(nx+1)]
#define vertexdy(k) vertexdy[((j-xmin)) + (k-ymin)*(nx+1)]

#define vol_flux_x(j,k) vol_flux_x[((j-xmin)) + (k-ymin)*(nx+1)]
#define vol_flux_y(j,k) vol_flux_y[((j-xmin)) + (k-ymin)*(nx)]

#define mass_flux_x(j,k) mass_flux_x[((j-xmin)) + (k-ymin)*(nx+1)]
#define mass_flux_y(j,k) mass_flux_y[((j-xmin)) + (k-ymin)*(nx)]

#define volume_change(i,j) volume_change[((i)-ifirst(0)) + (j-ifirst(1))*vnx]

#define pre_vol(j,k) pre_vol[(j - xmin) + (k - ymin)*(nx)]
#define post_vol(j,k) post_vol[(j - xmin) + (k - ymin)*(nx)]
#define pre_mass(j,k) pre_mass[(j - xmin) + (k - ymin)*(nx)]
#define post_mass(j,k) post_mass[(j - xmin) + (k - ymin)*(nx)]
#define advec_vol(j,k) advec_vol[(j - xmin) + (k - ymin)*(nx)]
#define post_ener(j,k) post_ener[(j - xmin) + (k - ymin)*(nx)]
#define ener_flux(j,k) ener_flux[(j - xmin) + (k - ymin)*(nx)]

#define node_flux(j,k) node_flux[(j-xmin) + (k - ymin)*(nx+1)]
#define node_mass_post(j,k) node_mass_post[(j-xmin) + (k - ymin)*(nx+1)]
#define node_mass_pre(j,k) node_mass_pre[(j-xmin) + (k - ymin)*(nx+1)]
#define advec_vel(j,k) advec_vel[(j-xmin) + (k - ymin)*(nx+1)]
#define mom_flux(j,k) mom_flux[(j-xmin) + (k - ymin)*(nx+1)]

#define vel1(j,k) vel1[((j-xmin)) + (k-ymin)*(nx+1)]

#define data(i,j) data[((i-xmin)) + (j-ymin)*nx]

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
    integrator->registerVariable(
            d_velocity,
            LagrangianEulerianIntegrator::FIELD,
            LagrangianEulerianIntegrator::PRIME_CELLS_EXCH |
                LagrangianEulerianIntegrator::PRE_LAGRANGE_EXCH |
                LagrangianEulerianIntegrator::PRE_SWEEP_1_MOM_EXCH |
                LagrangianEulerianIntegrator::PRE_SWEEP_2_MOM_EXCH,
            d_nghosts,
            d_grid_geometry);

    integrator->registerVariable(
            d_massflux,
            LagrangianEulerianIntegrator::NORMAL,
            LagrangianEulerianIntegrator::PRE_SWEEP_1_MOM_EXCH |
                LagrangianEulerianIntegrator::PRE_SWEEP_2_MOM_EXCH,
            d_nghosts,
            d_grid_geometry);

    integrator->registerVariable(
            d_volflux,
            LagrangianEulerianIntegrator::NORMAL,
            LagrangianEulerianIntegrator::PRE_SWEEP_1_CELL_EXCH,
            d_nghosts,
            d_grid_geometry);

    integrator->registerVariable(
            d_pressure,
            LagrangianEulerianIntegrator::NORMAL,
            LagrangianEulerianIntegrator::PRIME_CELLS_EXCH |
                LagrangianEulerianIntegrator::PRE_LAGRANGE_EXCH |
                LagrangianEulerianIntegrator::HALF_STEP_EXCH,
            d_nghosts,
            d_grid_geometry);

    integrator->registerVariable(
            d_viscosity,
            LagrangianEulerianIntegrator::NORMAL,
            LagrangianEulerianIntegrator::PRIME_CELLS_EXCH |
                LagrangianEulerianIntegrator::POST_VISCOSITY_EXCH,
            d_nghosts,
            d_grid_geometry);

    integrator->registerVariable(
            d_soundspeed,
            LagrangianEulerianIntegrator::NORMAL,
            LagrangianEulerianIntegrator::NO_EXCH,
            d_nghosts,
            d_grid_geometry);

    integrator->registerVariable(
            d_density,
            LagrangianEulerianIntegrator::FIELD | LagrangianEulerianIntegrator::REVERT,
            LagrangianEulerianIntegrator::PRIME_CELLS_EXCH |
                LagrangianEulerianIntegrator::PRE_LAGRANGE_EXCH |
                LagrangianEulerianIntegrator::PRE_SWEEP_1_CELL_EXCH |
                LagrangianEulerianIntegrator::PRE_SWEEP_1_MOM_EXCH |
                LagrangianEulerianIntegrator::PRE_SWEEP_2_MOM_EXCH,
            d_nghosts, 
            d_grid_geometry);

    integrator->registerVariable(
            d_energy,
            LagrangianEulerianIntegrator::FIELD | LagrangianEulerianIntegrator::REVERT,
            LagrangianEulerianIntegrator::PRIME_CELLS_EXCH |
                LagrangianEulerianIntegrator::PRE_LAGRANGE_EXCH |
                LagrangianEulerianIntegrator::PRE_SWEEP_1_CELL_EXCH |
                LagrangianEulerianIntegrator::PRE_SWEEP_1_MOM_EXCH |
                LagrangianEulerianIntegrator::PRE_SWEEP_2_MOM_EXCH,
            d_nghosts,
            d_grid_geometry);

    integrator->registerVariable(
            d_volume,
            LagrangianEulerianIntegrator::NORMAL,
            LagrangianEulerianIntegrator::NO_EXCH,
            d_nghosts,
            d_grid_geometry);

    integrator->registerVariable(
            d_celldeltas,
            LagrangianEulerianIntegrator::NORMAL,
            LagrangianEulerianIntegrator::NO_EXCH,
            d_nghosts,
            d_grid_geometry);
    integrator->registerVariable(
            d_cellcoords,
            LagrangianEulerianIntegrator::NORMAL,
            LagrangianEulerianIntegrator::NO_EXCH,
            d_nghosts,
            d_grid_geometry);

    integrator->registerVariable(
            d_vertexdeltas,
            LagrangianEulerianIntegrator::NORMAL,
            LagrangianEulerianIntegrator::NO_EXCH,
            d_nghosts,
            d_grid_geometry);

    integrator->registerVariable(
            d_vertexcoords,
            LagrangianEulerianIntegrator::NORMAL,
            LagrangianEulerianIntegrator::NO_EXCH,
            d_nghosts,
            d_grid_geometry);

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
         * Fill in the volume array...
         */
        const hier::Index ifirst = patch.getBox().lower();
        const hier::Index ilast = patch.getBox().upper();

        hier::IntVector density_ghosts = v_density->getGhostCellWidth();

        int xmin = ifirst(0) - density_ghosts(0);
        int xmax = ilast(0) + density_ghosts(0);
        int ymin = ifirst(1) - density_ghosts(1);
        int ymax = ilast(1) + density_ghosts(1);

        int xminng = ifirst(0);
        int xmaxng = ilast(0);
        int yminng = ifirst(1);
        int ymaxng = ilast(1);

        int nx = xmax - xmin + 1;
        int ny = ymax - ymin + 1;

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

                /*
                 * Start at rxmin-2dx because we have 2 ghost cells!
                 */
                vertexx[ind] = (rxmin-2*dx) + dx*(xcount);
                vertexy[ind] = (rymin-2*dy) + dy*(ycount);

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

        /*
         * Fill density and energy with some data, these are our initial conditions.
         */
        v_energy->fillAll(50);
        v_density->fillAll(50);

        double* density = v_density->getPointer();
        double* energy = v_energy->getPointer();

        for(int j = ymin; j <= ymax; j++) {
            for(int i = xmin; i <= xmax; i++) {
                int n1 = POLY2(i,j,xmin,ymin,nx);
                int v1 = POLY2(i,j,vimin,vjmin,vnx);

                /*
                 * Produces square of size 60x60 in the centre of the domain
                 */
//                if ((cellx[n1] >= 0.0 && cellx[n1] <= 5.0)
//                        && (celly[n1] >= 0.0 && celly[n1] <= 2.0)) {
//                    density(i,j) = 1.0;
//                    energy(i,j) = 2.5;
//                } else {
//                    density(i,j) = 0.2;
//                    energy(i,j) = 1.0;
//                }


                if ((vertexx[v1] >= 0.0 && vertexx[v1] < 2.0) && 
                    (vertexy[v1] >= 0.0 && vertexy[v1] < 5.0)) {
                        density(i,j) = 1.0;
                        energy(i,j) = 2.5;
                } else {
                    density(i,j) = 0.2;
                    energy(i,j) = 1.0;
                }


                /*
                 * Produces square of size 60x60 in the centre of the domain
                 */
//                if ((cellx[n1] <= -10.0))
//                {
//                    density(i,j) = 1.0;
//                    energy(i,j) = 2.5;
//                } else {
//                    density(i,j) = 0.1;
//                    energy(i,j) = 1.0;
//                }

                /*
                 * Use this loop to set square size per patch explicitly
                 */
//                if ((j >= yminng+1 && j <= ymaxng-1)
//                        && (i >= xminng+1 && i <= xmaxng-1)) {
//                    density(i,j) = 1.0;
//                    energy(i,j) = 2.5;
//                } else {
//                    density(i,j) = 0.1;
//                    energy(i,j) = 1.0;
//                }
                
            }
        }

        ideal_gas_knl(patch, false);
    }
}

void Cleverleaf::accelerate(
        hier::Patch& patch,
        double dt)
{
    double nodal_mass;

    tbox::Pointer<pdat::CellData<double> > v_density = patch.getPatchData(d_density, getCurrentDataContext());
    tbox::Pointer<pdat::CellData<double> > v_volume = patch.getPatchData(d_volume, getCurrentDataContext());
    tbox::Pointer<pdat::CellData<double> > v_celldeltas = patch.getPatchData(d_celldeltas, getCurrentDataContext());
    tbox::Pointer<pdat::CellData<double> > v_pressure = patch.getPatchData(d_pressure, getCurrentDataContext());
    tbox::Pointer<pdat::CellData<double> > v_viscosity = patch.getPatchData(d_viscosity, getCurrentDataContext());
    tbox::Pointer<pdat::NodeData<double> > v_vel0 = patch.getPatchData(d_velocity, getCurrentDataContext());
    tbox::Pointer<pdat::NodeData<double> > v_vel1 = patch.getPatchData(d_velocity, getNewDataContext());
        
    hier::IntVector ghosts = v_pressure->getGhostCellWidth();

    const hier::Index ifirst = patch.getBox().lower();
    const hier::Index ilast = patch.getBox().upper();

    int xmin = ifirst(0) - ghosts(0); 
    int xmax = ilast(0) + ghosts(0); 
    int ymin = ifirst(1) - ghosts(1); 
    int ymax = ilast(1) + ghosts(1); 

    /**
     * nx needs to account for the number of ghosts or things go very wrong!
     */
    int nx = (xmax) - (xmin) + 1;

    /**
     * vnx is used by the volume_change array define at the top of this file
     */
    int sbmnx = ilast(0) - ifirst(0) + 1;

    //tbox::pout << "xmin: " << nx << std::endl;

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

    pdat::NodeData<double> v_stepbymass(patch.getBox(), 1, hier::IntVector(d_dim, 0));
    double* stepbymass = v_stepbymass.getPointer();

    for(int k = ifirst(1); k <= ilast(1)+1; k++) {
        for(int j = ifirst(0); j <= ilast(0)+1; j++ ) {

            nodal_mass=(density0(j-1,k-1)*volume(j-1,k-1)
                    +density0(j,k-1)*volume(j,k-1)
                    +density0(j,k)*volume(j,k)
                    +density0(j-1,k)*volume(j-1,k))
                    *0.25;

                stepbymass(j,k)=0.5*dt/nodal_mass;
        }
    }

    for(int k = ifirst(1); k <= ilast(1)+1; k++) {
        for(int j = ifirst(0); j <= ilast(0)+1; j++ ) {

            xvel1(j,k)=xvel0(j,k)-stepbymass(j,k)*(xarea(j,k)*(pressure(j,k)-pressure(j-1,k))
                    +xarea(j,k-1)*(pressure(j,k-1)-pressure(j-1,k-1)));
#ifdef DEBUG
            cout << xvel0(j,k) << std::endl;
            cout << stepbymass(j,k) << std::endl;
            cout << xarea(j,k) << std::endl;
            cout << pressure(j,k) << std::endl;
            cout << pressure(j-1,k) << std::endl;
            cout << xarea(j,k-1) << std::endl;
            cout << pressure(j,k-1) << std::endl;
            cout << pressure(j-1,k-1) << std::endl;
            cout << "xvel1(" << j << "," << k << ") = " << xvel1(j,k) << std::endl;
#endif
        }
    }

    for(int k = ifirst(1); k <= ilast(1)+1; k++) {
        for(int j = ifirst(0); j <= ilast(0)+1; j++ ) {

            yvel1(j,k)=yvel0(j,k)-stepbymass(j,k)*(yarea(j,k)*(pressure(j,k)-pressure(j,k-1))
                    +yarea(j-1,k)*(pressure(j-1,k)-pressure(j-1,k-1)));
#ifdef DEBUG
            cout << "yvel1(" << j << "," << k << ") = " << yvel1(j,k) << std::endl;
#endif
        }
    }

    for(int k = ifirst(1); k <= ilast(1)+1; k++) {
        for(int j = ifirst(0); j <= ilast(0)+1; j++ ) {

            xvel1(j,k)=xvel1(j,k)-stepbymass(j,k)*(xarea(j,k)*(viscosity(j,k)-viscosity(j-1,k))
                    +xarea(j,k-1)*(viscosity(j,k-1)-viscosity(j-1,k-1)));

#ifdef DEBUG
            cout << "xvel1(" << j << "," << k << ") = " << xvel1(j,k) << std::endl;
#endif
        }
    }

    for(int k = ifirst(1); k <= ilast(1)+1; k++) {
        for(int j = ifirst(0); j <= ilast(0)+1; j++ ) {

            yvel1(j,k)=yvel1(j,k)-stepbymass(j,k)*(yarea(j,k)*(viscosity(j,k)-viscosity(j,k-1))
                    +yarea(j-1,k)*(viscosity(j-1,k)-viscosity(j-1,k-1)));
#ifdef DEBUG
            cout << "yvel1(" << j << "," << k << ") = " << yvel1(j,k) << std::endl;
#endif
        }
    }
}

void Cleverleaf::ideal_gas_knl(
        hier::Patch& patch,
        bool predict)
{

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


    hier::IntVector ghosts = v_pressure->getGhostCellWidth();

    const hier::Index ifirst = patch.getBox().lower();
    const hier::Index ilast = patch.getBox().upper();

    int xmin = ifirst(0) - ghosts(0); 
    int xmax = ilast(0) + ghosts(0); 
    int ymin = ifirst(1) - ghosts(1); 
    int ymax = ilast(1) + ghosts(1); 

    /**
     * nx needs to account for the number of ghosts or things go very wrong!
     */
    int nx = (xmax) - (xmin) + 1;

    //tbox::pout << "xmin: " << nx << std::endl;

    double* density = v_density->getPointer();
    double* energy = v_energy->getPointer();
    double* pressure = v_pressure->getPointer();
    double* soundspeed = v_soundspeed->getPointer();

    double pressurebyenergy = 0;
    double pressurebyvolume = 0;
    double sound_speed_squared = 0;
    double v = 0;

    for(int k = ifirst(1); k <= ilast(1); k++) {
        for(int j = ifirst(0); j <= ilast(0); j++){

            v=1.0/density(j,k);

            //tbox::pout << "density: " << density(j,k) << std::endl;
            //tbox::pout << "energy: " << energy(j,k) << std::endl;

            pressure(j,k)=(1.4-1.0)*density(j,k)*energy(j,k);

            pressurebyenergy=(1.4-1.0)*density(j,k);

            pressurebyvolume=-density(j,k)*pressure(j,k);

            sound_speed_squared=v*v*(pressure(j,k)*pressurebyenergy-pressurebyvolume);

            soundspeed(j,k)=sqrt(sound_speed_squared);

#ifdef DEBUG
            tbox::pout << "Updating pressure[" << j << "][" << k << "]=" << pressure(j,k) << ", soundspeed[" << j << "][" << k << "]=" << soundspeed(j,k) << std::endl;
#endif
        }
    }
}

void Cleverleaf::viscosity_knl(
        hier::Patch& patch)
{

  double ugrad,vgrad,grad2,pgradx,pgrady,
         pgradx2,pgrady2,grad,ygrad,pgrad,
         xgrad,div,strain2,limiter;

    tbox::Pointer<pdat::CellData<double> > v_density0 = patch.getPatchData(d_density, getCurrentDataContext());

    tbox::Pointer<pdat::CellData<double> > v_pressure = patch.getPatchData(d_pressure, getCurrentDataContext());

    tbox::Pointer<pdat::CellData<double> > v_viscosity = patch.getPatchData(d_viscosity, getCurrentDataContext());
    tbox::Pointer<pdat::NodeData<double> > v_vel0 = patch.getPatchData(d_velocity, getCurrentDataContext());
    tbox::Pointer<pdat::CellData<double> > v_celldeltas = patch.getPatchData(d_celldeltas, getCurrentDataContext());


    hier::IntVector ghosts = v_pressure->getGhostCellWidth();

    const hier::Index ifirst = patch.getBox().lower();
    const hier::Index ilast = patch.getBox().upper();

    int xmin = ifirst(0) - ghosts(0); 
    int xmax = ilast(0) + ghosts(0); 
    int ymin = ifirst(1) - ghosts(1); 
    int ymax = ilast(1) + ghosts(1); 

    /**
     * nx needs to account for the number of ghosts or things go very wrong!
     */
    int nx = (xmax) - (xmin) + 1;

    double* density = v_density0->getPointer();
    double* pressure = v_pressure->getPointer();
    double* viscosity = v_viscosity->getPointer();
    double* xvel0 = v_vel0->getPointer(0);
    double* yvel0 = v_vel0->getPointer(1);
    double* celldx = v_celldeltas->getPointer(0);
    double* celldy = v_celldeltas->getPointer(1);

    for(int k = ifirst(1); k <= ilast(1); k++) {
        for(int j = ifirst(0); j <= ilast(0); j++){

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
    double div,dsx,dsy,dtut,dtvt,dtct,dtdivt,cc,dv1,dv2;
    double xl_pos,
           yl_pos;

    int kldt;
    int jldt;

    tbox::Pointer<pdat::CellData<double> > v_density1 = patch.getPatchData(d_density, getNewDataContext());
    tbox::Pointer<pdat::CellData<double> > v_density0 = patch.getPatchData(d_density, getCurrentDataContext());

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

    hier::IntVector ghosts = v_pressure->getGhostCellWidth();

    const hier::Index ifirst = patch.getBox().lower();
    const hier::Index ilast = patch.getBox().upper();

    int xmin = ifirst(0) - ghosts(0);
    int xmax = ilast(0) + ghosts(0);
    int ymin = ifirst(1) - ghosts(1);
    int ymax = ilast(1) + ghosts(1);

    /**
     * nx needs to account for the number of ghosts or things go very wrong!
     */
    int nx = xmax - xmin + 1;

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

    /*
     * Original timestep
     */
//    for(int k = ifirst(1); k <= ilast(1); k++) {
//        for(int j = ifirst(0); j <= ilast(0); j++){
//
//            dsx=celldx(j,k);
//            dsy=celldy(j,k);
//
//            cc=soundspeed(j,k)*soundspeed(j,k);
//            cc=cc+2.0*viscosity(j,k)/density0(j,k);
//            cc=max(sqrt(cc),1.0e-16);
//
//            dtct=min(dsx,dsy)/cc;
//
//            div=0.0;
//
//            dv1=(xvel0(j  ,k)+xvel0(j  ,k+1))*xarea(j  ,k);
//            dv2=(xvel0(j+1,k)+xvel0(j+1,k+1))*xarea(j+1,k);
//
//            div=div+dv2-dv1;
//
//            dtut=2.0*volume(j,k)/max(abs(dv1),max(abs(dv2),1.0e-16*volume(j,k)));
//
//            dv1=(yvel0(j,k  )+yvel0(j+1,k  ))*yarea(j,k  );
//            dv2=(yvel0(j,k+1)+yvel0(j+1,k+1))*yarea(j,k+1);
//
//            div=div+dv2-dv1;
//
//            dtvt=2.0*volume(j,k)/max(abs(dv1),max(abs(dv2),1.0e-16*volume(j,k)));
//
//            div=div/(2.0*volume(j,k));
//
//            if (div < -1.0e-16) {
//                dtdivt=-1.0/div;
//            } else {
//                dtdivt=1.0e+21;
//            } 
//
//            if (dtct*dtc_safe < dt_min_val) {
//                jldt=j;
//                kldt=k;
//                dt_min_val=dtct*dtc_safe;
//                dtl_control=1;
//                xl_pos=cellx(j,k);
//                yl_pos=celly(j,k);
//            } 
//
//            if (dtut*dtu_safe < dt_min_val) {
//                jldt=j;
//                kldt=k;
//                dt_min_val=dtut*dtu_safe;
//                dtl_control=2;
//                xl_pos=cellx(j,k);
//                yl_pos=celly(j,k);
//            } 
//
//            if (dtvt*dtv_safe < dt_min_val) {
//                jldt=j;
//                kldt=k;
//                dt_min_val=dtvt*dtv_safe;
//                dtl_control=3;
//                xl_pos=cellx(j,k);
//                yl_pos=celly(j,k);
//            } 
//
//            if (dtdivt*dtdiv_safe < dt_min_val) {
//                jldt=j;
//                kldt=k;
//                dt_min_val=dtdivt*dtdiv_safe;
//                dtl_control=4;
//                xl_pos=cellx(j,k);
//                yl_pos=celly(j,k);
//            } 
//        }
//    }

    /*
     * new timestep
     */
    for(int k = ifirst(1); k <= ilast(1); k++) {
        for(int j = ifirst(0); j <= ilast(0); j++){

       dsx=celldx(j,k);
       dsy=celldy(j,k);

       cc=soundspeed(j,k)*soundspeed(j,k);
       cc=cc+2.0*viscosity(j,k)/density0(j,k);
       cc=max(sqrt(cc),1.0e-16);

       dtct=dtc_safe*min(dsx,dsy)/cc;

       div=0.0;

       dv1=(xvel0(j  ,k)+xvel0(j  ,k+1))*xarea(j  ,k);
       dv2=(xvel0(j+1,k)+xvel0(j+1,k+1))*xarea(j+1,k);

       div=div+dv2-dv1;

       dtut=dtu_safe*2.0*volume(j,k)/max(abs(dv1),max(abs(dv2),1.0e-16*volume(j,k)));

       dv1=(yvel0(j,k  )+yvel0(j+1,k  ))*yarea(j,k  );
       dv2=(yvel0(j,k+1)+yvel0(j+1,k+1))*yarea(j,k+1);

       div=div+dv2-dv1;

       dtvt=dtv_safe*2.0*volume(j,k)/max(abs(dv1),max(abs(dv2),1.0e-16*volume(j,k)));

       div=div/(2.0*volume(j,k));

       if (div < -1.0e-16) {
         dtdivt=dtdiv_safe*(-1.0/div);
       } else {
         dtdivt=1.0e+21;
       }

       dt_min_val=min(dtct,min(dtut,min(dtvt,min(dtdivt, dt_min_val))));

        }
    }

//    for(int k = ifirst(1); k <= ilast(1); k++) {
//        for(int j = ifirst(0); j <= ilast(0); j++){
//          if(dt_min(j,k) < dt_min_val)
//              dt_min_val=dt_min(j,k);
//        }
//    }

    tbox::pout << "Timestep information:" << std::endl;
//    tbox::pout << "\tj, k  : " << jldt << "," << kldt << std::endl;
//    tbox::pout << "\tx, y  : " << cellx(jldt,kldt) << "," << celly(jldt,kldt) << std::endl;
    tbox::pout << "\ttimestep : " << dt_min_val << std::endl;
//    tbox::pout << "\tCell velocities:" << std::endl;
//    tbox::pout << "\t\t" << xvel0(jldt  ,kldt  ) << "," << yvel0(jldt  ,kldt  ) << std::endl;
//    tbox::pout << "\t\t" << xvel0(jldt+1,kldt  ) << "," << yvel0(jldt+1,kldt  ) << std::endl;
//    tbox::pout << "\t\t" << xvel0(jldt+1,kldt+1) << "," << yvel0(jldt+1,kldt+1) << std::endl;
//    tbox::pout << "\t\t" << xvel0(jldt  ,kldt+1) << "," << yvel0(jldt  ,kldt+1) << std::endl;
//    tbox::pout << "\tdensity, energy, pressure, soundspeed " << std::endl;
//    tbox::pout << "\t" << density0(jldt, kldt) << "," << energy(jldt,kldt) << "," << pressure(jldt,kldt) << "," << soundspeed(jldt,kldt) << std::endl;

    //cout << "RETURNING " << dt_min_val << " FROM calc_dt_knl()" << endl;
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

    hier::IntVector ghosts = pres->getGhostCellWidth();

    const hier::Index ifirst = patch.getBox().lower();
    const hier::Index ilast = patch.getBox().upper();

    int xmin = ifirst(0) - ghosts(0);
    int xmax = ilast(0) + ghosts(0);
    int ymin = ifirst(1) - ghosts(1);
    int ymax = ilast(1) + ghosts(1);

    /**
     * nx needs to account for the number of ghosts or things go very wrong!
     */
    int nx = xmax- xmin + 1;

    /**
     * vnx is used by the volume_change array define at the top of this file
     */
    int vnx = ilast(0) - ifirst(0) + 1;


    double* xarea = area->getPointer(1);   
    double* yarea = area->getPointer(0);   
    double* volume = vol->getPointer(); 
    double* density0 = dens0->getPointer();
    double* density1  = dens1->getPointer();
    double* energy0 = en0->getPointer();
    double* energy1 = en1->getPointer();
    double* pressure = pres->getPointer();
    double* viscosity = visc->getPointer();
    double* xvel0 = v0->getPointer(0);
    double* xvel1 = v1->getPointer(0);
    double* yvel0 = v0->getPointer(1);
    double* yvel1 = v1->getPointer(1);

    pdat::CellData<double> v_volchange(patch.getBox(), 1, hier::IntVector(d_dim, 0));
    double* volume_change = v_volchange.getPointer();

    if (predict) {

        for (int  k = ifirst(1); k <= ilast(1); k++) {
            for (int j = ifirst(0); j <= ilast(0); j++) {

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

#ifdef DEBUG
                tbox::pout << "energy change = " << energy_change << ", volume change = " << volume_change(j,k) << std::endl;
                tbox::pout << "[" << j << "][" << k << "]" << std::endl;
#endif

                energy1(j,k)=energy0(j,k)-energy_change;

                density1(j,k)=density0(j,k)*volume_change(j,k);
            }
        }

    } else {

        for (int  k = ifirst(1); k <= ilast(1); k++) {
            for (int j = ifirst(0); j <= ilast(0); j++) {

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

#ifdef DEBUG
                tbox::pout << "energy change = " << energy_change << ", volume change = " << volume_change(j,k) << std::endl;
                tbox::pout << "[" << j << "][" << k << "]" << std::endl;
#endif

                energy1(j,k)=energy0(j,k)-energy_change;

                density1(j,k)=density0(j,k)*volume_change(j,k);
            }
        }

    }
}

void Cleverleaf::flux_calc_knl(
        hier::Patch& patch,
        double dt)
{

    tbox::Pointer<pdat::NodeData<double> > v_vel0 = patch.getPatchData(d_velocity, getCurrentDataContext());
    tbox::Pointer<pdat::NodeData<double> > v_vel1 = patch.getPatchData(d_velocity, getNewDataContext());

    tbox::Pointer<pdat::EdgeData<double> > v_volflux = patch.getPatchData(d_volflux, getCurrentDataContext());

    tbox::Pointer<pdat::CellData<double> > v_celldeltas = patch.getPatchData(d_celldeltas, getCurrentDataContext());

    hier::IntVector ghosts = v_celldeltas->getGhostCellWidth();

    const hier::Index ifirst = patch.getBox().lower();
    const hier::Index ilast = patch.getBox().upper();

    int xmin = ifirst(0) - ghosts(0); 
    int xmax = ilast(0) + ghosts(0); 
    int ymin = ifirst(1) - ghosts(1); 
    int ymax = ilast(1) + ghosts(1); 

    /**
     * nx needs to account for the number of ghosts or things go very wrong!
     */
    int nx = (xmax) - (xmin) + 1;

    double* xvel0 = v_vel0->getPointer(0);
    double* xvel1 = v_vel1->getPointer(0);
    double* xarea = v_celldeltas->getPointer(1);
    double* vol_flux_x = v_volflux->getPointer(1);

    double* yvel0 = v_vel0->getPointer(1);
    double* yvel1 = v_vel1->getPointer(1);
    double* yarea = v_celldeltas->getPointer(0);
    double* vol_flux_y = v_volflux->getPointer(0);

    for (int k = ifirst(1); k <= ilast(1); k++) {
        for (int j = ifirst(0); j <= ilast(0)+1; j++) {
            vol_flux_x(j,k)=0.25*dt*xarea(j,k)
                *(xvel0(j,k)+xvel0(j,k+1)+xvel1(j,k)+xvel1(j,k+1));

#ifdef DEBUG
            tbox::pout << "vol_flux_x(" << j << "," << k << ") = " << vol_flux_x(j,k) << std::endl;
            tbox::pout << "\tdt\txarea(" << j << "," << k << ")\txvel0(" << j << "," << k << ")\txvel0(" << j << "," << k+1 << ")\txvel1(" << j << "," << k << ")\txvel1(" << j << "," << k+1 << ")" << std::endl;
            tbox::pout << "\t" << dt << "\t" << xarea(j,k) << "\t" << xvel0(j,k) << "\t" << xvel0(j,k+1) << "\t" << xvel1(j,k) << "\t" << xvel1(j,k+1) << std::endl;
#endif
        }
    }

    for (int k = ifirst(1); k <= ilast(1)+1; k++) {
        for (int j = ifirst(0); j <= ilast(0); j++) {
            vol_flux_y(j,k)=0.25*dt*yarea(j,k)
                *(yvel0(j,k)+yvel0(j+1,k)+yvel1(j,k)+yvel1(j+1,k));

#ifdef DEBUG
            tbox::pout << "vol_flux_y(" << j << "," << k << ") = " << vol_flux_y(j,k) << std::endl;
#endif
        }
    }
}

void Cleverleaf::advec_cell(hier::Patch& patch,
        int sweep_number,
        ADVEC_DIR dir)
{
    tbox::Pointer<pdat::CellData<double> > v_density1 = patch.getPatchData(d_density, getNewDataContext());
    tbox::Pointer<pdat::CellData<double> > v_energy1 = patch.getPatchData(d_energy, getNewDataContext());
    tbox::Pointer<pdat::CellData<double> > v_volume = patch.getPatchData(d_volume, getCurrentDataContext());

    tbox::Pointer<pdat::EdgeData<double> > v_volflux = patch.getPatchData(d_volflux, getCurrentDataContext());
    tbox::Pointer<pdat::EdgeData<double> > v_massflux = patch.getPatchData(d_massflux, getCurrentDataContext());

    tbox::Pointer<pdat::NodeData<double> > vertexdeltas = patch.getPatchData(d_vertexdeltas, getCurrentDataContext());

    hier::IntVector ghosts = v_density1->getGhostCellWidth();

    const hier::Index ifirst = patch.getBox().lower();
    const hier::Index ilast = patch.getBox().upper();

    int xmin = ifirst(0) - ghosts(0); 
    int xmax = ilast(0) + ghosts(0); 
    int ymin = ifirst(1) - ghosts(1); 
    int ymax = ilast(1) + ghosts(1); 

    /**
     * nx needs to account for the number of ghosts or things go very wrong!
     */
    int nx = (xmax) - (xmin) + 1;

    double* volume = v_volume->getPointer();
    double* density1 = v_density1->getPointer();
    double* energy1 = v_energy1->getPointer();
    double* vol_flux_x = v_volflux->getPointer(1);
    double* vol_flux_y = v_volflux->getPointer(0);
    double* mass_flux_x = v_massflux->getPointer(1);
    double* mass_flux_y = v_massflux->getPointer(0);

    double* vertexdx = vertexdeltas->getPointer(0);
    double* vertexdy = vertexdeltas->getPointer(1);

    int upwind,donor,downwind,dif;

    double sigma,sigmat,sigmav,sigmam,sigma3,sigma4;
    double diffuw,diffdw,limiter;
    double one_by_six;

    pdat::CellData<double> v_prevol(patch.getBox(), 1, hier::IntVector(d_dim, 2));
    pdat::CellData<double> v_postvol(patch.getBox(), 1, hier::IntVector(d_dim, 2));
    pdat::CellData<double> v_premass(patch.getBox(), 1, hier::IntVector(d_dim, 2));
    pdat::CellData<double> v_postmass(patch.getBox(), 1, hier::IntVector(d_dim, 2));
    pdat::CellData<double> v_advecvol(patch.getBox(), 1, hier::IntVector(d_dim, 2));
    pdat::CellData<double> v_postener(patch.getBox(), 1, hier::IntVector(d_dim, 2));
    pdat::CellData<double> v_enerflux(patch.getBox(), 1, hier::IntVector(d_dim, 2));

    double* pre_vol = v_prevol.getPointer();
    double* post_vol = v_postvol.getPointer();
    double* pre_mass = v_premass.getPointer();
    double* post_mass = v_postmass.getPointer();
    double* advec_vol = v_advecvol.getPointer();
    double* post_ener = v_postener.getPointer();
    double* ener_flux = v_enerflux.getPointer();

    one_by_six=1.0/6.0;


    if(dir == X) {
        if (sweep_number == 1){
            for(int k= ymin; k <= ymax; k++) {
                for(int j=xmin; j <= xmax; j++) {
                    pre_vol(j,k)=volume(j,k)+(vol_flux_x(j+1,k  )-vol_flux_x(j,k)+vol_flux_y(j  ,k+1)-vol_flux_y(j,k));
                    post_vol(j,k)=pre_vol(j,k)-(vol_flux_x(j+1,k  )-vol_flux_x(j,k));
                }
            } 
        } else {
            for(int k=ymin; k <= ymax; k++) {
                for(int j = xmin; j <= xmax; j++) {
                    pre_vol(j,k)=volume(j,k)+vol_flux_x(j+1,k)-vol_flux_x(j,k);
                    post_vol(j,k)=volume(j,k);
                }
            } 
        }

        for(int k = ifirst(1); k <= ilast(1); k++) {
            for(int j = ifirst(0); j <= xmax; j++) {

                if(vol_flux_x(j,k)>0.0){
                    upwind   =j-2;
                    donor    =j-1;
                    downwind =j;
                    dif      =donor;
                } else {
                    upwind   =min(j+1,xmax);
                    donor    =j;
                    downwind =j-1;
                    dif      =upwind;
                }

                sigmat=abs(vol_flux_x(j,k))/pre_vol(donor,k);
                sigma3=(1.0+sigmat)*(vertexdx(j)/vertexdx(dif));
                sigma4=2.0-sigmat;

                sigma=sigmat;
                sigmav=sigmat;

                diffuw=density1(donor,k)-density1(upwind,k);
                diffdw=density1(downwind,k)-density1(donor,k);

                if(diffuw*diffdw>0.0){
                    limiter=(1.0-sigmav)*copysign(1.0,diffdw)*min(abs(diffuw),min(abs(diffdw),one_by_six*(sigma3*abs(diffuw)+sigma4*abs(diffdw))));
                } else {
                    limiter=0.0;
                }

                mass_flux_x(j,k)=vol_flux_x(j,k)*(density1(donor,k)+limiter);

                sigmam=abs(mass_flux_x(j,k))/(density1(donor,k)*pre_vol(donor,k));
                diffuw=energy1(donor,k)-energy1(upwind,k);
                diffdw=energy1(downwind,k)-energy1(donor,k);

                if(diffuw*diffdw>0.0){
                    limiter=(1.0-sigmam)*copysign(1.0,diffdw)*min(abs(diffuw),min(abs(diffdw),one_by_six*(sigma3*abs(diffuw)+sigma4*abs(diffdw))));
                } else {
                    limiter=0.0;
                }

                ener_flux(j,k)=mass_flux_x(j,k)*(energy1(donor,k)+limiter);
            }
        }

        for(int k=ifirst(1); k <= ilast(1); k++) {
            for(int j=ifirst(0); j <= ilast(0); j++) {
                pre_mass(j,k)=density1(j,k)*pre_vol(j,k);
                post_mass(j,k)=pre_mass(j,k)+mass_flux_x(j,k)-mass_flux_x(j+1,k);
                post_ener(j,k)=(energy1(j,k)*pre_mass(j,k)+ener_flux(j,k)-ener_flux(j+1,k))/post_mass(j,k);
                advec_vol(j,k)=pre_vol(j,k)+vol_flux_x(j,k)-vol_flux_x(j+1,k);
                density1(j,k)=post_mass(j,k)/advec_vol(j,k);
                energy1(j,k)=post_ener(j,k);
#ifdef DEBUG
                tbox::pout << "density(" << j << "," << k << ") = " << density1(j,k) << std::endl;
                tbox::pout << "energy(" << j << "," << k << ") = " << energy1(j,k) << std::endl;
                tbox::pout << "advec_vol(" << j << "," << k << ") = " << advec_vol(j,k) << std::endl;
#endif

            }
        }

    } else if (dir == Y) {

    if(sweep_number==1){

      for(int k=ymin; k <= ymax; k++) {
        for(int j=xmin; j <= xmax; j++) {
          pre_vol(j,k)=volume(j,k)+(vol_flux_y(j  ,k+1)-vol_flux_y(j,k)+vol_flux_x(j+1,k  )-vol_flux_x(j,k));
          post_vol(j,k)=pre_vol(j,k)-(vol_flux_y(j  ,k+1)-vol_flux_y(j,k));
        }
      }

    } else {

      for(int k=ymin; k <= ymax; k++) {
        for(int j=xmin; j <= xmax; j++) {
          pre_vol(j,k)=volume(j,k)+vol_flux_y(j  ,k+1)-vol_flux_y(j,k);
          post_vol(j,k)=volume(j,k);
        }
      }

    }

    for(int k=ifirst(1); k <= ymax; k++) {
      for(int j=ifirst(0); j <= ilast(0); j++) {

        if(vol_flux_y(j,k)>0.0){
          upwind   =k-2;
          donor    =k-1;
          downwind =k;
          dif      =donor;
        } else {
          upwind   =min(k+1,ymax);
          donor    =k;
          downwind =k-1;
          dif      =upwind;
        }

        sigmat=abs(vol_flux_y(j,k))/pre_vol(j,donor);
        sigma3=(1.0+sigmat)*(vertexdy(k)/vertexdy(dif));
        sigma4=2.0-sigmat;

        sigma=sigmat;
        sigmav=sigmat;

        diffuw=density1(j,donor)-density1(j,upwind);
        diffdw=density1(j,downwind)-density1(j,donor);

        if(diffuw*diffdw>0.0){
          limiter=(1.0-sigmav)*copysign(1.0,diffdw)*min(abs(diffuw),min(abs(diffdw),one_by_six*(sigma3*abs(diffuw)+sigma4*abs(diffdw))));
        } else {
          limiter=0.0;
        }
        mass_flux_y(j,k)=vol_flux_y(j,k)*(density1(j,donor)+limiter);

        sigmam=abs(mass_flux_y(j,k))/(density1(j,donor)*pre_vol(j,donor));
        diffuw=energy1(j,donor)-energy1(j,upwind);
        diffdw=energy1(j,downwind)-energy1(j,donor);

        if(diffuw*diffdw>0.0){
          limiter=(1.0-sigmam)*copysign(1.0,diffdw)*min(abs(diffuw),min(abs(diffdw),one_by_six*(sigma3*abs(diffuw)+sigma4*abs(diffdw))));
        } else {
          limiter=0.0;
        }
        ener_flux(j,k)=mass_flux_y(j,k)*(energy1(j,donor)+limiter);
      }
    }

    for( int k=ifirst(1); k <= ilast(1); k++) {
      for( int j=ifirst(0); j <= ilast(0); j++) {
        pre_mass(j,k)=density1(j,k)*pre_vol(j,k);
        post_mass(j,k)=pre_mass(j,k)+mass_flux_y(j,k)-mass_flux_y(j,k+1);
        post_ener(j,k)=(energy1(j,k)*pre_mass(j,k)+ener_flux(j,k)-ener_flux(j,k+1))/post_mass(j,k);
        advec_vol(j,k)=pre_vol(j,k)+vol_flux_y(j,k)-vol_flux_y(j,k+1);
        density1(j,k)=post_mass(j,k)/advec_vol(j,k);
        energy1(j,k)=post_ener(j,k);

#ifdef DEBUG
        tbox::pout << "advec_vol(" << j << "," << k << ") = " << advec_vol(j,k) << std::endl;
#endif
      }
    }
  }
}

void Cleverleaf::advec_mom(hier::Patch& patch,
        int sweep_number,
        ADVEC_DIR direction,
        ADVEC_DIR which_vel)
{
    tbox::Pointer<pdat::CellData<double> > v_density1 = patch.getPatchData(d_density, getNewDataContext());
    tbox::Pointer<pdat::CellData<double> > v_volume = patch.getPatchData(d_volume, getCurrentDataContext());
    tbox::Pointer<pdat::NodeData<double> > v_vel1 = patch.getPatchData(d_velocity, getNewDataContext());

    tbox::Pointer<pdat::EdgeData<double> > v_volflux = patch.getPatchData(d_volflux, getCurrentDataContext());
    tbox::Pointer<pdat::EdgeData<double> > v_massflux = patch.getPatchData(d_massflux, getCurrentDataContext());

    tbox::Pointer<pdat::CellData<double> > v_celldeltas = patch.getPatchData(d_celldeltas, getCurrentDataContext());

    hier::IntVector ghosts = v_density1->getGhostCellWidth();

    const hier::Index ifirst = patch.getBox().lower();
    const hier::Index ilast = patch.getBox().upper();

    int xmin = ifirst(0) - ghosts(0); 
    int xmax = ilast(0) + ghosts(0); 
    int ymin = ifirst(1) - ghosts(1); 
    int ymax = ilast(1) + ghosts(1); 

    /**
     * nx needs to account for the number of ghosts or things go very wrong!
     */
    int nx = (xmax) - (xmin) + 1;

    double* volume = v_volume->getPointer();
    double* density1 = v_density1->getPointer();
    double* vol_flux_x = v_volflux->getPointer(1);
    double* vol_flux_y = v_volflux->getPointer(0);
    double* mass_flux_x = v_massflux->getPointer(1);
    double* mass_flux_y = v_massflux->getPointer(0);

    double* vel1;

    double* celldx = v_celldeltas->getPointer(0);
    double* celldy = v_celldeltas->getPointer(1);

    int upwind,donor,downwind,dif;

    double sigma,width,wind;
    double vdiffuw,vdiffdw,auw,adw,limiter;

    pdat::NodeData<double> v_nodeflux(patch.getBox(), 1, hier::IntVector(d_dim, 2));
    pdat::NodeData<double> v_nodemasspost(patch.getBox(), 1, hier::IntVector(d_dim, 2));
    pdat::NodeData<double> v_nodemasspre(patch.getBox(), 1, hier::IntVector(d_dim, 2));
    pdat::NodeData<double> v_advecvel(patch.getBox(), 1, hier::IntVector(d_dim, 2));
    pdat::NodeData<double> v_momflux(patch.getBox(), 1, hier::IntVector(d_dim, 2));

    pdat::CellData<double> v_prevol(patch.getBox(), 1, hier::IntVector(d_dim, 2));
    pdat::CellData<double> v_postvol(patch.getBox(), 1, hier::IntVector(d_dim, 2));

    double* node_flux = v_nodeflux.getPointer();
    double* node_mass_post = v_nodemasspost.getPointer();
    double* node_mass_pre = v_nodemasspre.getPointer();
    double* advec_vel = v_advecvel.getPointer();
    double* mom_flux = v_momflux.getPointer();
    double* pre_vol = v_prevol.getPointer();
    double* post_vol = v_postvol.getPointer();

    if (which_vel == X) {
        vel1 = v_vel1->getPointer(0);
    } else {
        vel1 = v_vel1->getPointer(1);
    }

    int mom_sweep=direction+2*(sweep_number-1);

    tbox::perr << mom_sweep << std::endl;


    if (mom_sweep == 1) { //! x 1
        for(int k=ifirst(1)-2; k <= ilast(1)+2; k++) {
            for(int j=ifirst(0)-2; j <= ilast(0)+2; j++) {
                post_vol(j,k)= volume(j,k)+vol_flux_y(j  ,k+1)-vol_flux_y(j,k);
                pre_vol(j,k)=post_vol(j,k)+vol_flux_x(j+1,k  )-vol_flux_x(j,k);
            }
        }
    } else if (mom_sweep == 2) { //! y 1
        for(int k=ifirst(1)-2; k <= ilast(1)+2; k++) {
            for(int j=ifirst(0)-2; j <= ilast(0)+2; j++) {
                post_vol(j,k)= volume(j,k)+vol_flux_x(j+1,k  )-vol_flux_x(j,k);
                pre_vol(j,k)=post_vol(j,k)+vol_flux_y(j  ,k+1)-vol_flux_y(j,k);
            }
        }
    } else if (mom_sweep == 3) { //! x 2
        for(int k=ifirst(1)-2; k <= ilast(1)+2; k++) {
            for(int j=ifirst(0)-2; j <= ilast(0)+2; j++) {
                post_vol(j,k)=volume(j,k);
                pre_vol(j,k)=post_vol(j,k)+vol_flux_y(j  ,k+1)-vol_flux_y(j,k);
            }
        }
    } else if (mom_sweep == 4) { //! y 2
        for(int k=ifirst(1)-2; k <= ilast(1)+2; k++) {
            for(int j=ifirst(0)-2; j <= ilast(0)+2; j++) {
                post_vol(j,k)=volume(j,k);
                pre_vol(j,k)=post_vol(j,k)+vol_flux_x(j+1,k  )-vol_flux_x(j,k);
            }
        }
    } 

    if (direction == 1) {
        for(int k=ifirst(1); k <= ilast(1)+1; k++) {
            //! Find staggered mesh mass fluxes, nodal masses and volumes.
            for(int j=ifirst(0)-2; j <= ilast(0)+2; j++) {
                node_flux(j,k)=0.25*(mass_flux_x(j,k-1  )+mass_flux_x(j  ,k)+mass_flux_x(j+1,k-1)+mass_flux_x(j+1,k)); //! Mass Flux
            }
            for(int j=ifirst(0)-1; j <= ilast(0)+2; j++) {
                //! Staggered cell mass post advection
                node_mass_post(j,k)=0.25*(density1(j  ,k-1)*post_vol(j  ,k-1)                   
                        +density1(j  ,k  )*post_vol(j  ,k  )                   
                        +density1(j-1,k-1)*post_vol(j-1,k-1)                   
                        +density1(j-1,k  )*post_vol(j-1,k  ));
#ifdef DEBUG
                tbox::pout << "post_vol(" << j << ", " << k-1 << ") = " << post_vol(j,k-1) << std::endl;
                tbox::pout << "post_vol(" << j << ", " << k << ") = " << post_vol(j,k) << std::endl;
                tbox::pout << "post_vol(" << j-1 << ", " << k-1 << ") = " << post_vol(j-1,k-1) << std::endl;
                tbox::pout << "post_vol(" << j-1 << ", " << k << ") = " << post_vol(j-1,k) << std::endl;
                tbox::pout << "density1(" << j << ", " << k-1 << ") = " << density1(j,k-1) << std::endl;
                tbox::pout << "density1(" << j << ", " << k << ") = " << density1(j,k) << std::endl;
                tbox::pout << "density1(" << j-1 << ", " << k-1 << ") = " << density1(j-1,k-1) << std::endl;
                tbox::pout << "density1(" << j-1 << ", " << k << ") = " << density1(j-1,k) << std::endl;
#endif
            }
            //! Stagered cell mass pre advection
            for(int j=ifirst(0)-1; j <= ilast(0)+2; j++) {
                node_mass_pre(j,k)=node_mass_post(j,k)-node_flux(j-1,k)+node_flux(j,k);
            }
            for(int j=ifirst(0)-1; j <= ilast(0)+1; j++) {
                if (node_flux(j,k) < 0.0) {
                    upwind=j+2;
                    donor=j+1;
                    downwind=j;
                    dif=donor;
                } else {
                    upwind=j-1;
                    donor=j;
                    downwind=j+1;
                    dif=upwind;
                } 
                sigma=abs(node_flux(j,k))/(node_mass_pre(donor,k));
                width=celldx(j,k);
                vdiffuw=vel1(donor,k)-vel1(upwind,k);
                vdiffdw=vel1(downwind,k)-vel1(donor,k);
                limiter=0.0;
                if (vdiffuw*vdiffdw > 0.0) {
                    auw=abs(vdiffuw);
                    adw=abs(vdiffdw);
                    wind=1.0;
                    if (vdiffdw <= 0.0) wind=-1.0;
                    limiter=wind*min(width*((2.0-sigma)*adw/width+(1.0+sigma)*auw/celldx(dif,k))/6.0,min(auw,adw));
                } 
                //! Set advection velocity and mometum flux
                advec_vel(j,k)=vel1(donor,k)+(1.0-sigma)*limiter;
                mom_flux(j,k)=advec_vel(j,k)*node_flux(j,k);
            }
            for(int j=ifirst(0); j <= ilast(0)+1; j++) {
                double tv1 = (vel1(j,k)*node_mass_pre(j,k)+mom_flux(j-1,k)-mom_flux(j,k))/node_mass_post(j,k);
                double nmp = node_mass_pre(j,k);
                double mf = mom_flux(j-1,k);
                double mf1 = mom_flux(j,k);
                double nmpo = node_mass_post(j,k);

                vel1(j,k)=(vel1(j,k)*node_mass_pre(j,k)+mom_flux(j-1,k)-mom_flux(j,k))/node_mass_post(j,k);
#ifdef DEBUG
                tbox::pout << "node_mass_pre(" << j << "," << k << ") = " << node_mass_pre(j,k) << std::endl;
                tbox::pout << "mom_flux(" << j-1 << "," << k << ") = " << mom_flux(j-1,k) << std::endl;
                tbox::pout << "mom_flux(" << j << "," << k << ") = " << mom_flux(j,k) << std::endl;
                tbox::pout << "node_mass_post(" << j << "," << k << ") = " << node_mass_post(j,k) << std::endl;
                tbox::pout << "vel1(" << j << "," << k << ") = " << vel1(j,k) << std::endl;
#endif
            }
        }
    } else if (direction == 2) {
        for(int j=ifirst(0); j <= ilast(0)+1; j++) {
            //! Find staggered mesh mass fluxes and nodal masses and volumes.
            for(int k=ifirst(1)-2; k <= ilast(1)+2; k++) {
                node_flux(j,k)=0.25*(mass_flux_y(j-1,k  )+mass_flux_y(j  ,k  )+mass_flux_y(j-1,k+1)+mass_flux_y(j  ,k+1));
            }
            for(int k=ifirst(1)-1; k <= ilast(1)+2; k++) {
                node_mass_post(j,k)=0.25*(density1(j  ,k-1)*post_vol(j  ,k-1)                     
                        +density1(j  ,k  )*post_vol(j  ,k  )                     
                        +density1(j-1,k-1)*post_vol(j-1,k-1)                     
                        +density1(j-1,k  )*post_vol(j-1,k  ));
            }
            for(int k=ifirst(1)-1; k <= ilast(1)+2; k++) {
                node_mass_pre(j,k)=node_mass_post(j,k)-node_flux(j,k-1)+node_flux(j,k);
            }
            for(int k=ifirst(1)-1; k <= ilast(1)+1; k++) {
                if (node_flux(j,k) < 0.0) {
                    upwind=k+2;
                    donor=k+1;
                    downwind=k;
                    dif=donor;
                } else {
                    upwind=k-1;
                    donor=k;
                    downwind=k+1;
                    dif=upwind;
                } 

                sigma=abs(node_flux(j,k))/(node_mass_pre(j,donor));
                width=celldy(j,k);
                vdiffuw=vel1(j,donor)-vel1(j,upwind);
                vdiffdw=vel1(j,downwind)-vel1(j,donor);
                limiter=0.0;
                if (vdiffuw*vdiffdw > 0.0) {
                    auw=abs(vdiffuw);
                    adw=abs(vdiffdw);
                    wind=1.0;
                    if (vdiffdw <= 0.0) wind=-1.0;
                    limiter=wind*min(width*((2.0-sigma)*adw/width+(1.0+sigma)*auw/celldy(j,dif))/6.0,min(auw,adw));
                } 
                advec_vel(j,k)=vel1(j,donor)+(1.0-sigma)*limiter;
                mom_flux(j,k)=advec_vel(j,k)*node_flux(j,k);
            }

            for(int k=ifirst(1); k <= ilast(1)+1; k++) {
                double tv1 = (vel1(j,k)*node_mass_pre(j,k)+mom_flux(j,k-1)-mom_flux(j,k))/node_mass_post(j,k);
                double nmp = node_mass_pre(j,k);
                double mf = mom_flux(j,k-1);
                double mf1 = mom_flux(j,k);
                double nmpo = node_mass_post(j,k);

                vel1(j,k)=(vel1(j,k)*node_mass_pre(j,k)+mom_flux(j,k-1)-mom_flux(j,k))/node_mass_post(j,k);
#ifdef DEBUG
                tbox::pout << "node_mass_pre(" << j << "," << k << ") = " << node_mass_pre(j,k) << std::endl;
                tbox::pout << "mom_flux(" << j << "," << k-1 << ") = " << mom_flux(j-1,k) << std::endl;
                tbox::pout << "mom_flux(" << j << "," << k << ") = " << mom_flux(j,k) << std::endl;
                tbox::pout << "node_mass_post(" << j << "," << k << ") = " << node_mass_post(j,k) << std::endl;
                tbox::pout << "vel1(" << j << "," << k << ") = " << vel1(j,k) << std::endl;
#endif
            }
        }
    } 
}


void Cleverleaf::setPhysicalBoundaryConditions(
        hier::Patch& patch,
        const double fill_time,
        const hier::IntVector& ghost_width_to_fill)
{

    //tbox::pout << "In Cleverleaf::setPhysicalBoundaryConditions..." << std::endl;

    tbox::Pointer<pdat::CellData<double> > v_pressure =
        patch.getPatchData(d_pressure, getScratchDataContext());

    tbox::Pointer<pdat::CellData<double> > v_density0 =
        patch.getPatchData(d_density, getScratchDataContext());

    tbox::Pointer<pdat::CellData<double> > v_density1 =
        patch.getPatchData(d_density, getNewDataContext());

    tbox::Pointer<pdat::CellData<double> > v_energy0 =
        patch.getPatchData(d_energy, getScratchDataContext());
    
    tbox::Pointer<pdat::CellData<double> > v_energy1 =
        patch.getPatchData(d_energy, getNewDataContext());

    tbox::Pointer<pdat::CellData<double> > v_viscosity =
        patch.getPatchData(d_viscosity, getScratchDataContext());

    tbox::Pointer<pdat::NodeData<double> > v_vel0 =
        patch.getPatchData(d_velocity, getScratchDataContext());

    tbox::Pointer<pdat::NodeData<double> > v_vel1 =
        patch.getPatchData(d_velocity, getNewDataContext());

    tbox::Pointer<pdat::EdgeData<double> > v_massflux = 
        patch.getPatchData(d_massflux, getScratchDataContext());

    tbox::Pointer<pdat::EdgeData<double> > v_volflux = 
        patch.getPatchData(d_volflux, getScratchDataContext());

    hier::IntVector ghosts = v_pressure->getGhostCellWidth();

    const hier::Index ifirst = patch.getBox().lower();
    const hier::Index ilast = patch.getBox().upper();

    int xmin = ifirst(0) - ghosts(0); 
    int xmax = ilast(0) + ghosts(0); 
    int ymin = ifirst(1) - ghosts(1); 
    int ymax = ilast(1) + ghosts(1); 

    int nx = xmax - xmin +1;

    double* pressure = v_pressure->getPointer();

    double* density0 = v_density0->getPointer();
    double* density1 = v_density1->getPointer();

    double* energy0 = v_energy0->getPointer();
    double* energy1 = v_energy1->getPointer();

    double* viscosity = v_viscosity->getPointer();

    double* xvel0 = v_vel0->getPointer(0);
    double* xvel1 = v_vel1->getPointer(0);

    double* yvel0 = v_vel0->getPointer(1);
    double* yvel1 = v_vel1->getPointer(1);

    double* mass_flux_x = v_massflux->getPointer(1);
    double* mass_flux_y = v_massflux->getPointer(0);

    double* vol_flux_x = v_volflux->getPointer(1);
    double* vol_flux_y = v_volflux->getPointer(0);

    int depth = ghost_width_to_fill[0];

    const tbox::Pointer<geom::CartesianPatchGeometry> pgeom = 
        patch.getPatchGeometry();

    const tbox::Array<hier::BoundaryBox>& edge_bdry = pgeom->getCodimensionBoundaries(Bdry::EDGE2D);

    for(int i = 0; i < edge_bdry.getSize(); i++) {
        switch(edge_bdry[i].getLocationIndex()) {
            case (BdryLoc::YLO) :

                /*
                 * update pressure boundary...
                 */
                reflectPhysicalBoundary(
                        pressure,
                        BdryLoc::YLO,
                        depth,
                        ifirst, ilast,
                        xmin, xmax, ymin, ymax,
                        nx);

                reflectPhysicalBoundary(
                        density0,
                        BdryLoc::YLO,
                        depth,
                        ifirst, ilast,
                        xmin, xmax, ymin, ymax,
                        nx);

                reflectPhysicalBoundary(
                        density1,
                        BdryLoc::YLO,
                        depth,
                        ifirst, ilast,
                        xmin, xmax, ymin, ymax, nx);

                reflectPhysicalBoundary(
                        energy0,
                        BdryLoc::YLO,
                        depth,
                        ifirst, ilast,
                        xmin, xmax, ymin, ymax, nx);

                reflectPhysicalBoundary(
                        energy1,
                        BdryLoc::YLO,
                        depth,
                        ifirst, ilast,
                        xmin, xmax, ymin, ymax, nx);

                reflectPhysicalBoundary(
                        viscosity,
                        BdryLoc::YLO,
                        depth,
                        ifirst, ilast,
                        xmin, xmax, ymin, ymax, nx);

                reflectXNodeBoundary(
                        xvel0,
                        BdryLoc::YLO,
                        depth,
                        ifirst, ilast,
                        xmin, xmax+1, ymin, ymax+1, nx+1);

                reflectXNodeBoundary(
                        xvel1,
                        BdryLoc::YLO,
                        depth,
                        ifirst, ilast,
                        xmin, xmax+1, ymin, ymax+1, nx+1);

                reflectYNodeBoundary(
                        yvel0,
                        BdryLoc::YLO,
                        depth,
                        ifirst, ilast,
                        xmin, xmax+1, ymin, ymax+1, nx+1);

                reflectYNodeBoundary(
                        yvel1,
                        BdryLoc::YLO,
                        depth,
                        ifirst, ilast,
                        xmin, xmax+1, ymin, ymax+1, nx+1);

                reflectXEdgeBoundary(
                        vol_flux_x,
                        BdryLoc::YLO,
                        depth,
                        ifirst, ilast,
                        xmin, xmax+1, ymin, ymax, nx+1);

                reflectYEdgeBoundary(
                        vol_flux_y,
                        BdryLoc::YLO,
                        depth,
                        ifirst, ilast,
                        xmin, xmax, ymin, ymax+1, nx);

                reflectXEdgeBoundary(
                        mass_flux_x,
                        BdryLoc::YLO,
                        depth,
                        ifirst, ilast,
                        xmin, xmax+1, ymin, ymax, nx+1);

                reflectYEdgeBoundary(
                        mass_flux_y,
                        BdryLoc::YLO,
                        depth,
                        ifirst, ilast,
                        xmin, xmax, ymin, ymax+1, nx);
                break;


            case (BdryLoc::YHI) :
                
                /*
                 * Update pressure boundary...
                 */
                reflectPhysicalBoundary(
                        pressure,
                        BdryLoc::YHI,
                        depth,
                        ifirst, ilast,
                        xmin, xmax, ymin, ymax,
                        nx);

                reflectPhysicalBoundary(
                        density0,
                        BdryLoc::YHI,
                        depth,
                        ifirst, ilast,
                        xmin, xmax, ymin, ymax,
                        nx);

                reflectPhysicalBoundary(
                        density1,
                        BdryLoc::YHI,
                        depth,
                        ifirst, ilast,
                        xmin, xmax, ymin, ymax,
                        nx);

                reflectPhysicalBoundary(
                        energy0,
                        BdryLoc::YHI,
                        depth,
                        ifirst, ilast,
                        xmin, xmax, ymin, ymax,
                        nx);

                reflectPhysicalBoundary(
                        energy1,
                        BdryLoc::YHI,
                        depth,
                        ifirst, ilast,
                        xmin, xmax, ymin, ymax,
                        nx);

                reflectPhysicalBoundary(
                        viscosity,
                        BdryLoc::YHI,
                        depth,
                        ifirst, ilast,
                        xmin, xmax, ymin, ymax,
                        nx);

                reflectXNodeBoundary(
                        xvel0,
                        BdryLoc::YHI,
                        depth,
                        ifirst, ilast,
                        xmin, xmax+1, ymin, ymax+1, nx+1);

                reflectXNodeBoundary(
                        xvel1,
                        BdryLoc::YHI,
                        depth,
                        ifirst, ilast,
                        xmin, xmax+1, ymin, ymax+1, nx+1);

                reflectYNodeBoundary(
                        yvel0,
                        BdryLoc::YHI,
                        depth,
                        ifirst, ilast,
                        xmin, xmax+1, ymin, ymax+1, nx+1);

                reflectYNodeBoundary(
                        yvel1,
                        BdryLoc::YHI,
                        depth,
                        ifirst, ilast,
                        xmin, xmax+1, ymin, ymax+1, nx+1);

                reflectXEdgeBoundary(
                        vol_flux_x,
                        BdryLoc::YHI,
                        depth,
                        ifirst, ilast,
                        xmin, xmax+1, ymin, ymax, nx+1);

                reflectYEdgeBoundary(
                        vol_flux_y,
                        BdryLoc::YHI,
                        depth,
                        ifirst, ilast,
                        xmin, xmax, ymin, ymax+1, nx);

                reflectXEdgeBoundary(
                        mass_flux_x,
                        BdryLoc::YHI,
                        depth,
                        ifirst, ilast,
                        xmin, xmax+1, ymin, ymax, nx+1);

                reflectYEdgeBoundary(
                        mass_flux_y,
                        BdryLoc::YHI,
                        depth,
                        ifirst, ilast,
                        xmin, xmax, ymin, ymax+1, nx);
                break;


            case (BdryLoc::XLO) :

                reflectPhysicalBoundary(
                        pressure,
                        BdryLoc::XLO,
                        depth,
                        ifirst, ilast,
                        xmin, xmax, ymin, ymax,
                        nx);

                reflectPhysicalBoundary(
                        density0,
                        BdryLoc::XLO,
                        depth,
                        ifirst, ilast,
                        xmin, xmax, ymin, ymax, nx);

                reflectPhysicalBoundary(
                        density1,
                        BdryLoc::XLO,
                        depth,
                        ifirst, ilast,
                        xmin, xmax, ymin, ymax,
                        nx);

                reflectPhysicalBoundary(
                        energy0,
                        BdryLoc::XLO,
                        depth,
                        ifirst, ilast,
                        xmin, xmax, ymin, ymax,
                        nx);

                reflectPhysicalBoundary(
                        energy1,
                        BdryLoc::XLO,
                        depth,
                        ifirst, ilast,
                        xmin, xmax, ymin, ymax,
                        nx);

                reflectPhysicalBoundary(
                        viscosity,
                        BdryLoc::XLO,
                        depth,
                        ifirst, ilast,
                        xmin, xmax, ymin, ymax,
                        nx);

                reflectXNodeBoundary(
                        xvel0,
                        BdryLoc::XLO,
                        depth,
                        ifirst, ilast,
                        xmin, xmax+1, ymin, ymax+1, nx+1);

                reflectXNodeBoundary(
                        xvel1,
                        BdryLoc::XLO,
                        depth,
                        ifirst, ilast,
                        xmin, xmax+1, ymin, ymax+1, nx+1);

                reflectYNodeBoundary(
                        yvel0,
                        BdryLoc::XLO,
                        depth,
                        ifirst, ilast,
                        xmin, xmax+1, ymin, ymax+1, nx+1);

                reflectYNodeBoundary(
                        yvel1,
                        BdryLoc::XLO,
                        depth,
                        ifirst, ilast,
                        xmin, xmax+1, ymin, ymax+1, nx+1);

                reflectXEdgeBoundary(
                        vol_flux_x,
                        BdryLoc::XLO,
                        depth,
                        ifirst, ilast,
                        xmin, xmax+1, ymin, ymax, nx+1);

                reflectYEdgeBoundary(
                        vol_flux_y,
                        BdryLoc::XLO,
                        depth,
                        ifirst, ilast,
                        xmin, xmax, ymin, ymax+1, nx);

                reflectXEdgeBoundary(
                        mass_flux_x,
                        BdryLoc::XLO,
                        depth,
                        ifirst, ilast,
                        xmin, xmax+1, ymin, ymax, nx+1);

                reflectYEdgeBoundary(
                        mass_flux_y,
                        BdryLoc::XLO,
                        depth,
                        ifirst, ilast,
                        xmin, xmax, ymin, ymax+1, nx);

                //tbox::pout << "XLO" << std::endl;
                break;

            case (BdryLoc::XHI) :

                reflectPhysicalBoundary(
                        pressure,
                        BdryLoc::XHI,
                        depth,
                        ifirst, ilast,
                        xmin, xmax, ymin, ymax,
                        nx);

                reflectPhysicalBoundary(
                        density0,
                        BdryLoc::XHI,
                        depth,
                        ifirst, ilast,
                        xmin, xmax, ymin, ymax,
                        nx);

                reflectPhysicalBoundary(
                        density1,
                        BdryLoc::XHI,
                        depth,
                        ifirst, ilast,
                        xmin, xmax, ymin, ymax,
                        nx);

                reflectPhysicalBoundary(
                        energy0,
                        BdryLoc::XHI,
                        depth,
                        ifirst, ilast,
                        xmin, xmax, ymin, ymax,
                        nx);

                reflectPhysicalBoundary(
                        energy1,
                        BdryLoc::XHI,
                        depth,
                        ifirst, ilast,
                        xmin, xmax, ymin, ymax,
                        nx);

                reflectPhysicalBoundary(
                        viscosity,
                        BdryLoc::XHI,
                        depth,
                        ifirst, ilast,
                        xmin, xmax, ymin, ymax,
                        nx);

                reflectXNodeBoundary(
                        xvel0,
                        BdryLoc::XHI,
                        depth,
                        ifirst, ilast,
                        xmin, xmax+1, ymin, ymax+1, nx+1);

                reflectXNodeBoundary(
                        xvel1,
                        BdryLoc::XHI,
                        depth,
                        ifirst, ilast,
                        xmin, xmax+1, ymin, ymax+1, nx+1);

                reflectYNodeBoundary(
                        yvel0,
                        BdryLoc::XHI,
                        depth,
                        ifirst, ilast,
                        xmin, xmax+1, ymin, ymax+1, nx+1);

                reflectYNodeBoundary(
                        yvel1,
                        BdryLoc::XHI,
                        depth,
                        ifirst, ilast,
                        xmin, xmax+1, ymin, ymax+1, nx+1);

                reflectXEdgeBoundary(
                        vol_flux_x,
                        BdryLoc::XHI,
                        depth,
                        ifirst, ilast,
                        xmin, xmax+1, ymin, ymax, nx+1);

                reflectYEdgeBoundary(
                        vol_flux_y,
                        BdryLoc::XHI,
                        depth,
                        ifirst, ilast,
                        xmin, xmax, ymin, ymax+1, nx);

                reflectXEdgeBoundary(
                        mass_flux_x,
                        BdryLoc::XHI,
                        depth,
                        ifirst, ilast,
                        xmin, xmax+1, ymin, ymax, nx+1);

                reflectYEdgeBoundary(
                        mass_flux_y,
                        BdryLoc::XHI,
                        depth,
                        ifirst, ilast,
                        xmin, xmax, ymin, ymax+1, nx);

                //tbox::pout << "XHI" << std::endl;
                break;

            default : tbox::perr << "[ERROR] Unknown edge location in setPhysicalBoundaryConditions... " << std::endl;
                      exit(-1);
        }
    }

    //tbox::pout << "Leaving Cleverleaf::setPhysicalBoundaryConditions..." << std::endl;

}

void Cleverleaf::reflectPhysicalBoundary(
        double* data,
        BdryLoc::Type boundary,
        int depth,
        hier::Index ifirst,
        hier::Index ilast,
        int xmin,
        int xmax,
        int ymin,
        int ymax,
        int nx)
{

        switch(boundary) {
            case (BdryLoc::YLO) :
                /*
                 * Reflect bottom edge...
                 */
                for (int k=1; k <= depth; k++) {
                    for (int j=xmin; j <= xmax; j++) {
                        data(j, ifirst(1)-k) = data(j, (ifirst(1)+(k-1)));
                    }
                }
                break;

            case (BdryLoc::YHI) :
                
                /*
                 * Reflect top edge...
                 */
                for (int k=1; k <= depth; k++) {
                    for (int j=xmin; j <= xmax; j++) {
                        data(j,ilast(1)+k)=data(j,ilast(1)-(k-1));
                    }
                }
                break;

            case (BdryLoc::XLO) :

                for (int k=ymin; k <= ymax; k++) {
                    for (int j=1; j <= depth; j++) {
                        data(ifirst(0)-j,k)=data(ifirst(0)+(j-1),k);
                    }
                }
                break;

            case (BdryLoc::XHI) :
                for (int k=ymin; k <= ymax; k++) {
                    for (int j=1; j <= depth; j++) {
                        data(ilast(0)+j,k)=data(ilast(0)-(j-1),k);
                    }
                }
                break;

            default : tbox::perr << "[ERROR] Unknown edge location in reflectPhysicalBoundary... " << std::endl;
                      exit(-1);
        }
}


void Cleverleaf::reflectXNodeBoundary(
        double* data,
        BdryLoc::Type boundary,
        int depth,
        hier::Index ifirst,
        hier::Index ilast,
        int xmin,
        int xmax,
        int ymin,
        int ymax,
        int nx)
{
        switch(boundary) {
            case (BdryLoc::YLO) :
                /*
                 * Reflect bottom edge...
                 */
                for (int k=1; k <= depth; k++) {
                    for (int j= ifirst(0)-depth; j <= ilast(0)+depth; j++) {
                        data(j, ifirst(1)-k) = data(j, (ifirst(1)+(k)));
                    }
                }
                break;

            case (BdryLoc::YHI) :
                
                /*
                 * Reflect top edge...
                 */
                for (int k=1; k <= depth; k++) {
                    for (int j=ifirst(0)-depth; j <= ilast(0)+1+depth; j++) {
                        data(j,ilast(1)+1+k) = data(j,ilast(1)+1-k);
                    }
                }
                break;

            case (BdryLoc::XLO) :

                for (int k=ifirst(1)-depth; k <= ilast(1)+1+depth; k++) {
                    for (int j=1; j <= depth; j++) {
                        data(ifirst(0)-j,k)= -data(ifirst(0)+j,k);

                        //std::cerr << "q(" << j << "," << ifirst(1)-k << ") = " << data(j,ifirst(1)-k) << ", setting q(" << j << "," << ifirst(1)+(k-1) << ") = " << -data(j, (ifirst(1)+(k-1))) << std::endl;
                    }
                }
                break;

            case (BdryLoc::XHI) :
                for (int k=ifirst(1)-depth; k <= ilast(1)+1+depth; k++) {
                    for (int j=1; j <= depth; j++) {
                        data(ilast(0)+1+j,k)= -data(ilast(0)+1-j,k);
                    }
                }
                break;

            default : tbox::perr << "[ERROR] Unknown edge location in reflectXQuantBoundary... " << std::endl;
                      exit(-1);
        }
}

void Cleverleaf::reflectYNodeBoundary(
        double* data,
        BdryLoc::Type boundary,
        int depth,
        hier::Index ifirst,
        hier::Index ilast,
        int xmin,
        int xmax,
        int ymin,
        int ymax,
        int nx)
{
        switch(boundary) {
            case (BdryLoc::YLO) :
                /*
                 * Reflect bottom edge...
                 */
                for (int k=1; k <= depth; k++) {
                    for (int j=ifirst(0)-depth; j <= ilast(0)+1+depth; j++) {
                        data(j, ifirst(1)-k) = -data(j, ifirst(1)+k);
                    }
                }
                break;

            case (BdryLoc::YHI) :
                
                /*
                 * Reflect top edge...
                 */
                for (int k=1; k <= depth; k++) {
                    for (int j=ifirst(0)-depth; j <= ilast(0)+1+depth; j++) {
                        data(j,ilast(1)+1+k) = -data(j,ilast(1)+1-k);
                    }
                }
                break;

            case (BdryLoc::XLO) :

                for (int k=ifirst(1)-depth; k <= ilast(1)+1+depth; k++) {
                    for (int j=1; j <= depth; j++) {
                        data(ifirst(0)-j,k) = data(ifirst(0)+j,k);
                    }
                }
                break;

            case (BdryLoc::XHI) :
                for (int k=ifirst(1)-depth; k <= ilast(1)+1+depth; k++) {
                    for (int j=1; j <= depth; j++) {
                        data(ilast(0)+1+j,k) = data(ilast(0)+1-j,k);
                    }
                }
                break;

            default : tbox::perr << "[ERROR] Unknown edge location in reflectYQuantBoundary... " << std::endl;
                      exit(-1);
        }
}

void Cleverleaf::reflectXEdgeBoundary(
        double* data,
        BdryLoc::Type boundary,
        int depth,
        hier::Index ifirst,
        hier::Index ilast,
        int xmin,
        int xmax,
        int ymin,
        int ymax,
        int nx)
{
        switch(boundary) {
            case (BdryLoc::YLO) :
                for (int k=1; k <= depth; k++) {
                    for (int j= ifirst(0)-depth; j <= ilast(0)+1+depth; j++) {
                        data(j, ifirst(1)-k) = data(j, (ifirst(1)+(k)));
                    }
                }
                break;

            case (BdryLoc::YHI) :
                for (int k=1; k <= depth; k++) {
                    for (int j=ifirst(0)-depth; j <= ilast(0)+1+depth; j++) {
                        data(j,ilast(1)+k) = data(j,ilast(1)-k);
                    }
                }
                break;

            case (BdryLoc::XLO) :
                for (int k=ifirst(1)-depth; k <= ilast(1)+depth; k++) {
                    for (int j=1; j <= depth; j++) {
                        data(ifirst(0)-j,k)= -data(ifirst(0)+j,k);
                    }
                }
                break;

            case (BdryLoc::XHI) :
                for (int k=ifirst(1)-depth; k <= ilast(1)+depth; k++) {
                    for (int j=1; j <= depth; j++) {
                        data(ilast(0)+1+j,k)= -data(ilast(0)+1-j,k);
                    }
                }
                break;

            default : tbox::perr << "[ERROR] Unknown edge location in reflectXQuantBoundary... " << std::endl;
                      exit(-1);
        }
}

void Cleverleaf::reflectYEdgeBoundary(
        double* data,
        BdryLoc::Type boundary,
        int depth,
        hier::Index ifirst,
        hier::Index ilast,
        int xmin,
        int xmax,
        int ymin,
        int ymax,
        int nx)
{
        switch(boundary) {
            case (BdryLoc::YLO) :
                /*
                 * Reflect bottom edge...
                 */
                for (int k=1; k <= depth; k++) {
                    for (int j=ifirst(0)-depth; j <= ilast(0)+depth; j++) {
                        data(j, ifirst(1)-k) = -data(j, ifirst(1)+k);
                    }
                }
                break;

            case (BdryLoc::YHI) :
                
                /*
                 * Reflect top edge...
                 */
                for (int k=1; k <= depth; k++) {
                    for (int j=ifirst(0)-depth; j <= ilast(0)+depth; j++) {
                        data(j,ilast(1)+1+k) = -data(j,ilast(1)+1-k);
                    }
                }
                break;

            case (BdryLoc::XLO) :

                for (int k=ifirst(1)-depth; k <= ilast(1)+1+depth; k++) {
                    for (int j=1; j <= depth; j++) {
                        data(ifirst(0)-j,k) = data(ifirst(0)+j,k);
                    }
                }
                break;

            case (BdryLoc::XHI) :
                for (int k=ifirst(1)-depth; k <= ilast(1)+1+depth; k++) {
                    for (int j=1; j <= depth; j++) {
                        data(ilast(0)+j,k) = data(ilast(0)-j,k);
                    }
                }
                break;

            default : tbox::perr << "[ERROR] Unknown edge location in reflectYQuantBoundary... " << std::endl;
                      exit(-1);
        }
}

void Cleverleaf::field_summary(
        hier::Patch& patch,
        double* vol,
        double* mass,
        double* press,
        double* ie,
        double* ke)
{
    tbox::Pointer<pdat::CellData<double> > v_volume = patch.getPatchData(d_volume, getCurrentDataContext());
    tbox::Pointer<pdat::CellData<double> > v_density0 = patch.getPatchData(d_density, getCurrentDataContext());
    tbox::Pointer<pdat::CellData<double> > v_energy0 = patch.getPatchData(d_energy, getCurrentDataContext());
    tbox::Pointer<pdat::CellData<double> > v_pressure = patch.getPatchData(d_pressure, getCurrentDataContext());
    tbox::Pointer<pdat::NodeData<double> > v_vel0 = patch.getPatchData(d_velocity, getCurrentDataContext());

    const hier::Index ifirst = patch.getBox().lower();
    const hier::Index ilast = patch.getBox().upper();

    hier::IntVector ghosts = v_density0->getGhostCellWidth();

    int xmin = ifirst(0) - ghosts(0);
    int xmax = ilast(0) + ghosts(0);
    int ymin = ifirst(1) - ghosts(1);
    int ymax = ilast(1) + ghosts(1);

    int nx = xmax - xmin + 1;
    int ny = ymax - ymin + 1;

    double* volume = v_volume->getPointer();
    double* density0 = v_density0->getPointer();
    double* energy0 = v_energy0->getPointer();
    double* pressure = v_pressure->getPointer();
    double* xvel0 = v_vel0->getPointer(0);
    double* yvel0 = v_vel0->getPointer(1);

    double vsqrd;
    double cell_vol;
    double cell_mass;

    *vol=0.0;
    *mass=0.0;
    *ie=0.0;
    *ke=0.0;
    *press=0.0;

    for(int k = ifirst(1); k <= ilast(1); k++) {
        for(int j = ifirst(0); j <= ilast(0); j++) {
            vsqrd=0.0;

            vsqrd=vsqrd+0.25*(pow(xvel0(j,k),2)+pow(yvel0(j,k),2));
            vsqrd=vsqrd+0.25*(pow(xvel0(j+1,k),2)+pow(yvel0(j+1,k),2));
            vsqrd=vsqrd+0.25*(pow(xvel0(j,k+1),2)+pow(yvel0(j,k+1),2));
            vsqrd=vsqrd+0.25*(pow(xvel0(j+1,k+1),2)+pow(yvel0(j+1,k+1),2));

            cell_vol=volume(j,k);
            cell_mass=cell_vol*density0(j,k);

            *vol=*vol+cell_vol;
            *mass=*mass+cell_mass;
            *ie=*ie+cell_mass*energy0(j,k);
            *ke=*ke+cell_mass*0.5*vsqrd;
            *press=*press+cell_vol*pressure(j,k);
        }
    }
}
