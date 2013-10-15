#include "Cleverleaf.h"
#include "CartesianCellDoubleVolumeWeightedAverage.h"
#include "CartesianCellDoubleMassWeightedAverage.h"
#include "CartesianCellIntConstantCoarsen.h"

#include "SAMRAI/hier/VariableDatabase.h"
#include "SAMRAI/geom/CartesianPatchGeometry.h"

#include "SAMRAI/pdat/CellData.h"
#include "SAMRAI/pdat/NodeData.h"
#include "SAMRAI/pdat/EdgeData.h"

#include "SAMRAI/appu/CartesianBoundaryUtilities2.h"
#include "SAMRAI/appu/CartesianBoundaryDefines.h"

#include "SAMRAI/pdat/NodeDoubleInjection.h"
#include "SAMRAI/geom/CartesianNodeDoubleLinearRefine.h"
#include "SAMRAI/geom/CartesianEdgeDoubleConservativeLinearRefine.h"
#include "SAMRAI/geom/CartesianCellDoubleConservativeLinearRefine.h"

#include <iostream>
#include <cmath>

#define F90_FUNC(name,NAME) name ## _

extern "C" {
    void F90_FUNC(ideal_gas_kernel,IDEAL_GAS_KERNEL)
        (int*,int*,int*,int*,double*,double*,double*,double*);

    void F90_FUNC(accelerate_kernel, ACCELERATE_KERNEL)
        (int*,int*,int*,int*,double*,double*,double*,double*,double*,double*,
         double*,double*,double*,double*,double*,double*);

    void F90_FUNC(viscosity_kernel, VISCOSITY_KERNEL)
        (int*,int*,int*,int*,double*,double*,double*,double*,double*,double*,double*);

    void F90_FUNC(pdv_kernel, PDV_KERNEL)
        (int*,int*,int*,int*,int*,double*,double*,double*,double*,double*,double*,
         double*,double*,double*,double*,double*,double*,double*,double*,double*);

    void F90_FUNC(advec_cell_kernel, ADVEC_CELL_KERNEL)
        (int*,int*,int*,int*,int*,int*,double*,double*,double*,double*,double*,double*,
         double*,double*,double*,double*,double*,double*,double*,double*,double*,double*);

    void F90_FUNC(advec_mom_kernel, ADVEC_MOM_KERNEL)
        (int*,int*,int*,int*,double*,double*,double*,double*,double*,double*,double*,double*,
         double*,double*,double*,double*,double*,double*,double*,double*,double*,int*,int*,int*);

    void F90_FUNC(calc_dt_kernel, CALC_DT_KERNEL)
        (int*,int*,int*,int*,double*,double*,double*,double*,double*,double*,double*,double*,
         double*,double*,double*,double*,double*,double*,double*,double*,double*,double*,
         double*,double*,double*,double*,double*,int*,double*,double*,int*,int*,int*);

    void F90_FUNC(flux_calc_kernel, FLUX_CALC_KERNEL)
        (int*,int*,int*,int*,double*,double*,double*,double*,double*,double*,
            double*,double*,double*);

    void F90_FUNC(update_halo_kernel_top, UPDATE_HALO_KERNEL_TOP)
        (int*,int*,int*,int*,double*,double*,double*,double*,double*,double*,
         double*,double*,double*,double*,double*,double*,double*,double*,double*,int*,int*);

    void F90_FUNC(update_halo_kernel_bottom, UPDATE_HALO_KERNEL_BOTTOM)
        (int*,int*,int*,int*,double*,double*,double*,double*,double*,double*,
         double*,double*,double*,double*,double*,double*,double*,double*,double*,int*,int*);

    void F90_FUNC(update_halo_kernel_left, UPDATE_HALO_KERNEL_LEFT)
        (int*,int*,int*,int*,double*,double*,double*,double*,double*,double*,
         double*,double*,double*,double*,double*,double*,double*,double*,double*,int*,int*);

    void F90_FUNC(update_halo_kernel_right, UPDATE_HALO_KERNEL_RIGHT)
        (int*,int*,int*,int*,double*,double*,double*,double*,double*,double*,
         double*,double*,double*,double*,double*,double*,double*,double*,double*,int*,int*);

    void F90_FUNC(tag_q_kernel,TAG_Q_KERNEL)
        (int*,int*,int*,int*,double*,int*);

    void F90_FUNC(tag_energy_kernel,TAG_ENERGY_KERNEL)
        (int*,int*,int*,int*,double*,int*);

    void F90_FUNC(tag_density_kernel,TAG_DENSITY_KERNEL)
        (int*,int*,int*,int*,double*,int*);

    void F90_FUNC(tag_all_kernel,TAG_ALL_KERNEL)
        (int*,int*,int*,int*,int*);

    void F90_FUNC(debug_kernel, DEBUG_KERNEL)
        (int*, int*, int*, int*, double*, double*, double*, double*);

    void F90_FUNC(field_summary_kernel, FIELD_SUMMARY_KERNEL)
        (int*, int*, int*, int*, double*, double*, double*, double*, double*, double*, int*,
         double*, double*, double*, double*, double*, int*);
}


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

#define data(i,j) data[(((i)-xmin)) + (((j)-ymin)*nx)]

Cleverleaf::Cleverleaf(
        boost::shared_ptr<tbox::Database> input_database,
        boost::shared_ptr<hier::PatchHierarchy> hierarchy,
        const tbox::Dimension& dim,
        boost::shared_ptr<geom::CartesianGridGeometry> grid_geometry):
    LagrangianEulerianPatchStrategy(dim),
    d_dim(dim),
    d_nghosts(dim),
    input_db(input_database),
    state_prefix("state"),
    d_velocity(new pdat::NodeVariable<double>(d_dim, "velocity", d_dim.getValue())),
    d_massflux(new pdat::EdgeVariable<double>(d_dim, "massflux", 2)),
    d_volflux(new pdat::EdgeVariable<double>(d_dim, "volflux", 2)),
    d_pressure(new pdat::CellVariable<double>(d_dim, "pressure", 1)),
    d_viscosity(new pdat::CellVariable<double>(d_dim, "viscosity", 1)),
    d_soundspeed(new pdat::CellVariable<double>(d_dim, "soundspeed", 1)),
    d_density(new pdat::CellVariable<double>(d_dim, "density", 1)),
    d_energy(new pdat::CellVariable<double>(d_dim, "energy", 1)),
    d_volume(new pdat::CellVariable<double>(d_dim, "volume", 1)),
    d_celldeltas(new pdat::CellVariable<double>(d_dim, "celldelta", d_dim.getValue())),
    d_cellcoords(new pdat::CellVariable<double>(d_dim, "cellcoords", d_dim.getValue())),
    d_vertexdeltas(new pdat::NodeVariable<double>(d_dim, "vertexdeltas", d_dim.getValue())),
    d_vertexcoords(new pdat::NodeVariable<double>(d_dim, "vertexcoords", d_dim.getValue())),
    d_level_indicator(new pdat::CellVariable<int>(d_dim, "level_indicator", d_dim.getValue())),
    d_exchange_fields(new int[15])
{
    d_hierarchy = hierarchy;
    d_grid_geometry = grid_geometry;

    d_nghosts = hier::IntVector(d_dim, 2);

    /*
     * Add our coarsen operators to the registry.
     */
    boost::shared_ptr<hier::CoarsenOperator> vol_weighted_avg(new CartesianCellDoubleVolumeWeightedAverage(dim));
    boost::shared_ptr<hier::CoarsenOperator> mass_weighted_avg(new CartesianCellDoubleMassWeightedAverage(dim));
    boost::shared_ptr<hier::CoarsenOperator> constant_cell_coarsen(new CartesianCellIntConstantCoarsen(dim));

    boost::shared_ptr<hier::CoarsenOperator> ndi(new pdat::NodeDoubleInjection());
    boost::shared_ptr<hier::RefineOperator> cndlr(new geom::CartesianNodeDoubleLinearRefine());
    boost::shared_ptr<hier::RefineOperator> cedclr(new geom::CartesianEdgeDoubleConservativeLinearRefine());
    boost::shared_ptr<hier::RefineOperator> ccdclr(new geom::CartesianCellDoubleConservativeLinearRefine());

    d_grid_geometry->addCoarsenOperator(typeid(pdat::CellVariable<double>).name(), vol_weighted_avg);
    d_grid_geometry->addCoarsenOperator(typeid(pdat::CellVariable<double>).name(), mass_weighted_avg);
    d_grid_geometry->addCoarsenOperator(typeid(pdat::CellVariable<int>).name(), constant_cell_coarsen);
    d_grid_geometry->addCoarsenOperator(typeid(pdat::NodeVariable<double>).name(), ndi);
    d_grid_geometry->addRefineOperator(typeid(pdat::NodeVariable<double>).name(), cndlr);
    d_grid_geometry->addRefineOperator(typeid(pdat::EdgeVariable<double>).name(), cedclr);
    d_grid_geometry->addRefineOperator(typeid(pdat::CellVariable<double>).name(), ccdclr);

    d_tag_all = input_database->getBoolWithDefault("tag_all", false);
}

void Cleverleaf::registerModelVariables(
        LagrangianEulerianLevelIntegrator* integrator) 
{
    integrator->registerVariable(
            d_velocity,
            LagrangianEulerianLevelIntegrator::FIELD,
            LagrangianEulerianLevelIntegrator::PRIME_CELLS_EXCH |
                LagrangianEulerianLevelIntegrator::PRE_LAGRANGE_EXCH |
                LagrangianEulerianLevelIntegrator::PRE_SWEEP_1_MOM_EXCH |
                LagrangianEulerianLevelIntegrator::PRE_SWEEP_2_MOM_EXCH,
            d_nghosts,
            d_grid_geometry);

    integrator->registerVariable(
            d_massflux,
            LagrangianEulerianLevelIntegrator::NORMAL,
            LagrangianEulerianLevelIntegrator::PRE_SWEEP_1_MOM_EXCH |
                LagrangianEulerianLevelIntegrator::PRE_SWEEP_2_MOM_EXCH,
            d_nghosts,
            d_grid_geometry);

    integrator->registerVariable(
            d_volflux,
            LagrangianEulerianLevelIntegrator::NORMAL,
            LagrangianEulerianLevelIntegrator::PRE_SWEEP_1_CELL_EXCH,
            d_nghosts,
            d_grid_geometry);

    integrator->registerVariable(
            d_pressure,
            LagrangianEulerianLevelIntegrator::NORMAL,
            LagrangianEulerianLevelIntegrator::PRIME_CELLS_EXCH |
                LagrangianEulerianLevelIntegrator::PRE_LAGRANGE_EXCH |
                LagrangianEulerianLevelIntegrator::HALF_STEP_EXCH,
            d_nghosts,
            d_grid_geometry);

    integrator->registerVariable(
            d_viscosity,
            LagrangianEulerianLevelIntegrator::NORMAL,
            LagrangianEulerianLevelIntegrator::PRIME_CELLS_EXCH |
                LagrangianEulerianLevelIntegrator::POST_VISCOSITY_EXCH,
            d_nghosts,
            d_grid_geometry);

    integrator->registerVariable(
            d_soundspeed,
            LagrangianEulerianLevelIntegrator::NORMAL,
            LagrangianEulerianLevelIntegrator::NO_EXCH,
            d_nghosts,
            d_grid_geometry);

    integrator->registerVariable(
            d_density,
            LagrangianEulerianLevelIntegrator::FIELD | LagrangianEulerianLevelIntegrator::REVERT,
            LagrangianEulerianLevelIntegrator::PRIME_CELLS_EXCH |
                LagrangianEulerianLevelIntegrator::PRE_LAGRANGE_EXCH |
                LagrangianEulerianLevelIntegrator::PRE_SWEEP_1_CELL_EXCH |
                LagrangianEulerianLevelIntegrator::PRE_SWEEP_1_MOM_EXCH |
                LagrangianEulerianLevelIntegrator::PRE_SWEEP_2_MOM_EXCH,
            d_nghosts, 
            d_grid_geometry);

    integrator->registerVariable(
            d_energy,
            LagrangianEulerianLevelIntegrator::FIELD | LagrangianEulerianLevelIntegrator::REVERT,
            LagrangianEulerianLevelIntegrator::PRIME_CELLS_EXCH |
                LagrangianEulerianLevelIntegrator::PRE_LAGRANGE_EXCH |
                LagrangianEulerianLevelIntegrator::PRE_SWEEP_1_CELL_EXCH |
                LagrangianEulerianLevelIntegrator::PRE_SWEEP_1_MOM_EXCH |
                LagrangianEulerianLevelIntegrator::PRE_SWEEP_2_MOM_EXCH,
            d_nghosts,
            d_grid_geometry);

    integrator->registerVariable(
            d_volume,
            LagrangianEulerianLevelIntegrator::NORMAL,
            LagrangianEulerianLevelIntegrator::NO_EXCH,
            d_nghosts,
            d_grid_geometry);

    integrator->registerVariable(
            d_celldeltas,
            LagrangianEulerianLevelIntegrator::NORMAL,
            LagrangianEulerianLevelIntegrator::NO_EXCH,
            d_nghosts,
            d_grid_geometry);
    integrator->registerVariable(
            d_cellcoords,
            LagrangianEulerianLevelIntegrator::NORMAL,
            LagrangianEulerianLevelIntegrator::NO_EXCH,
            d_nghosts,
            d_grid_geometry);

    integrator->registerVariable(
            d_vertexdeltas,
            LagrangianEulerianLevelIntegrator::NORMAL,
            LagrangianEulerianLevelIntegrator::NO_EXCH,
            d_nghosts,
            d_grid_geometry);

    integrator->registerVariable(
            d_vertexcoords,
            LagrangianEulerianLevelIntegrator::NORMAL,
            LagrangianEulerianLevelIntegrator::NO_EXCH,
            d_nghosts,
            d_grid_geometry);

    integrator->registerVariable(
            d_level_indicator,
            LagrangianEulerianLevelIntegrator::INDICATOR,
            LagrangianEulerianLevelIntegrator::NO_EXCH,
            d_nghosts,
            d_grid_geometry);

    hier::VariableDatabase* vardb = hier::VariableDatabase::getDatabase();

    d_plot_context = integrator->getPlotContext();

    if (d_visit_writer) {
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

        d_visit_writer->registerPlotQuantity("Level Indicator",
                "SCALAR",
                vardb->mapVariableAndContextToIndex(
                    d_level_indicator, d_plot_context));

        /*
         * Register vectors with the VisIt writer.
         */
        d_visit_writer->registerPlotQuantity("Velocity",
                "VECTOR",
                vardb->mapVariableAndContextToIndex(
                    d_velocity, d_plot_context));

//        d_visit_writer->registerPlotQuantity("Massflux",
//                "VECTOR",
//                vardb->mapVariableAndContextToIndex(
//                    d_massflux, d_plot_context));
//
//        d_visit_writer->registerPlotQuantity("Volflux",
//                "VECTOR",
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

void Cleverleaf::registerVisItDataWriter(boost::shared_ptr<appu::VisItDataWriter> writer) {
    d_visit_writer = writer;
}

void Cleverleaf::initializeDataOnPatch(
        hier::Patch& patch,
        double init_data_time,
        bool initial_time)
{
    if (initial_time) {
        boost::shared_ptr<pdat::NodeData<double> > velocity(
                patch.getPatchData(d_velocity, getCurrentDataContext()),
                boost::detail::dynamic_cast_tag());

        boost::shared_ptr<pdat::EdgeData<double> > massflux(
                patch.getPatchData(d_massflux, getCurrentDataContext()),
                boost::detail::dynamic_cast_tag());

        boost::shared_ptr<pdat::EdgeData<double> > volflux(
                patch.getPatchData(d_volflux, getCurrentDataContext()),
                boost::detail::dynamic_cast_tag());

        boost::shared_ptr<pdat::CellData<double> > v_pressure(
                patch.getPatchData(d_pressure, getCurrentDataContext()),
                boost::detail::dynamic_cast_tag());

        boost::shared_ptr<pdat::CellData<double> > v_viscosity(
                patch.getPatchData(d_viscosity, getCurrentDataContext()),
                boost::detail::dynamic_cast_tag());

        boost::shared_ptr<pdat::CellData<double> > v_soundspeed(
                patch.getPatchData(d_soundspeed, getCurrentDataContext()),
                boost::detail::dynamic_cast_tag());

        boost::shared_ptr<pdat::CellData<double> > v_density(
                patch.getPatchData(d_density, getCurrentDataContext()),
                boost::detail::dynamic_cast_tag());

        boost::shared_ptr<pdat::CellData<double> > v_energy(
                patch.getPatchData(d_energy, getCurrentDataContext()),
                boost::detail::dynamic_cast_tag());

        boost::shared_ptr<pdat::CellData<double> > v_volume(
                patch.getPatchData(d_volume, getCurrentDataContext()),
                boost::detail::dynamic_cast_tag());

        boost::shared_ptr<pdat::CellData<double> > celldeltas(
                patch.getPatchData(d_celldeltas,getCurrentDataContext()),
                boost::detail::dynamic_cast_tag());

        boost::shared_ptr<pdat::CellData<double> > cellcoords(
                patch.getPatchData(d_cellcoords, getCurrentDataContext()),
                boost::detail::dynamic_cast_tag());

        boost::shared_ptr<pdat::NodeData<double> > vertexdeltas(
                patch.getPatchData(d_vertexdeltas, getCurrentDataContext()),
                boost::detail::dynamic_cast_tag());

        boost::shared_ptr<pdat::NodeData<double> > vertexcoords(
                patch.getPatchData( d_vertexcoords, getCurrentDataContext()),
                boost::detail::dynamic_cast_tag());

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

        const boost::shared_ptr<geom::CartesianPatchGeometry> pgeom(patch.getPatchGeometry(), boost::detail::dynamic_cast_tag());

        const double* dxs = pgeom->getDx();
        const double* coords = pgeom->getXLower();
        const double* ucoords = pgeom->getXUpper();

        double rxmin = coords[0];
        double rymin = coords[1];
        double rxmax = ucoords[0];
        double rymax = ucoords[1];

        double dx = dxs[0];
        double dy = dxs[1];
        double vol = dx*dy;

        //std::cout << "density (pre-fill) ==" << v_density->getPointer()[POLY2(1, 1, -2, -2, nx)] << std::endl;


        /*
         * Use the fillAll() methods to initialise other variables for now...
         */
        velocity->fillAll(0.0);
        massflux->fillAll(0.0);
        volflux->fillAll(0.0);
        v_viscosity->fillAll(0.0);
        v_soundspeed->fillAll(0.0);
        v_volume->fillAll(vol);
        v_pressure->fillAll(0.0);

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
                vertexx[ind] = (rxmin-2.0*dx) + dx*(xcount);
                vertexy[ind] = (rymin-2.0*dy) + dy*(ycount);

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
        v_energy->fillAll(0.0);
        v_density->fillAll(0.0);
        double* density = v_density->getPointer();
        double* energy = v_energy->getPointer();
        double* xvel0 = velocity->getPointer(0);
        double* yvel0 = velocity->getPointer(1);

        boost::shared_ptr<tbox::Database> states_db = input_db->getDatabase("states");
        int nstates = states_db->getInteger("num_states");

        for(int state = 0; state < nstates; state++) {
            std::ostringstream state_stream;
            state_stream << state_prefix << state;

            boost::shared_ptr<tbox::Database> current_state = states_db->getDatabase(state_stream.str());

            double state_xmin = rxmin;
            double state_xmax = rxmax;
            double state_ymin = rymin;
            double state_ymax = rymax;
            std::string state_geometry;

            state_geometry = current_state->getStringWithDefault("geometry", "RECTANGLE");

            double state_density = current_state->getDouble("density");
            double state_energy = current_state->getDouble("energy");
            double state_xvel = current_state->getDoubleWithDefault("xvel", 0.0);
            double state_yvel = current_state->getDoubleWithDefault("yvel", 0.0);

            if (state_geometry.compare("RECTANGLE") == 0) {

                if (state != 0) {
                    double* state_min = new double[2];
                    double* state_max = new double[2];

                    current_state->getDoubleArray("min", state_min, 2);
                    current_state->getDoubleArray("max", state_max, 2);

                    state_xmin = state_min[0];
                    state_ymin = state_min[1];

                    state_xmax = state_max[0];
                    state_ymax = state_max[1];
                }

                for(int j = ymin; j <= ymax; j++) {
                    for(int i = xmin; i <= xmax; i++) {
                        int n1 = POLY2(i,j,xmin,ymin,nx);
                        int v1 = POLY2(i,j,vimin,vjmin,vnx);

                        if (((float)vertexx[v1] >= (float)state_xmin && (float)vertexx[v1] < (float)state_xmax)) { 
                            if ((float)vertexy[v1] >= (float)state_ymin && (float)vertexy[v1] < (float)state_ymax) {
                                density(i,j) = state_density;
                                energy(i,j) = state_energy;

                                xvel0(i,j) = state_xvel;
                                xvel0(i,j+1) = state_xvel;
                                xvel0(i+1,j) = state_xvel;
                                xvel0(i+1,j+1) = state_xvel;

                                yvel0(i,j) = state_yvel;
                                yvel0(i,j+1) = state_yvel;
                                yvel0(i+1,j) = state_yvel;
                                yvel0(i+1,j+1) = state_yvel;
                            }
                        }
                    }
                }
            } else if (state_geometry.compare("CIRCLE") == 0) {
                double* center = new double[2];
                current_state->getDoubleArray("center", center, 2);

                double x_center = center[0];
                double y_center = center[1];

                double state_radius = current_state->getDouble("radius");

                for(int j = ymin; j <= ymax; j++) {
                    for(int i = xmin; i <= xmax; i++) {
                        int n1 = POLY2(i,j,xmin,ymin,nx);
                        int v1 = POLY2(i,j,vimin,vjmin,vnx);

                        double cell_radius =
                            sqrt((cellx[n1]-x_center)*(cellx[n1]-x_center)
                                    +(celly[n1]-y_center)*(celly[n1]-y_center));

                        if((float)cell_radius <= (float)state_radius) {


                            density(i,j) = state_density;
                            energy(i,j) = state_energy;

                            xvel0(i,j) = state_xvel;
                            xvel0(i,j+1) = state_xvel;
                            xvel0(i+1,j) = state_xvel;
                            xvel0(i+1,j+1) = state_xvel;

                            yvel0(i,j) = state_yvel;
                            yvel0(i,j+1) = state_yvel;
                            yvel0(i+1,j) = state_yvel;
                            yvel0(i+1,j+1) = state_yvel;
                        }
                    }
                }
            } else if (state_geometry.compare("POINT") == 0) {
                double* center = new double[2];
                current_state->getDoubleArray("center", center, 2);

                double x_center = center[0];
                double y_center = center[1];

                for(int j = ymin; j <= ymax; j++) {
                    for(int i = xmin; i <= xmax; i++) {
                        int n1 = POLY2(i,j,xmin,ymin,nx);
                        int v1 = POLY2(i,j,vimin,vjmin,vnx);

                        if ((float)vertexx[v1] == (float)x_center) { 
                            if ((float)vertexy[v1] == (float)y_center) {
                                density(i,j) = state_density;
                                energy(i,j) = state_energy;
                            }
                        }
                    }
                }
            }
        }
    } else {
        boost::shared_ptr<pdat::CellData<double> > celldeltas(patch.getPatchData(
                    d_celldeltas,
                    getCurrentDataContext()),
                boost::detail::dynamic_cast_tag());

        boost::shared_ptr<pdat::CellData<double> > cellcoords(
                patch.getPatchData(d_cellcoords, getCurrentDataContext()),
                boost::detail::dynamic_cast_tag());

        boost::shared_ptr<pdat::NodeData<double> > vertexdeltas(
                patch.getPatchData(d_vertexdeltas, getCurrentDataContext()),
                boost::detail::dynamic_cast_tag());

        boost::shared_ptr<pdat::NodeData<double> > vertexcoords(
                patch.getPatchData(d_vertexcoords, getCurrentDataContext()),
                boost::detail::dynamic_cast_tag());

        boost::shared_ptr<pdat::CellData<double> > v_volume(
                patch.getPatchData(d_volume, getCurrentDataContext()),
                boost::detail::dynamic_cast_tag());

        boost::shared_ptr<pdat::EdgeData<double> > v_volflux(
                patch.getPatchData(d_volflux, getCurrentDataContext()),
                boost::detail::dynamic_cast_tag());

        // Fill volflux with 0.0 for now
        v_volflux->fillAll(0.0);


        /*
         * Fill in the volume array...
         */
        const hier::Index ifirst = patch.getBox().lower();
        const hier::Index ilast = patch.getBox().upper();

        hier::IntVector ghosts = celldeltas->getGhostCellWidth();

        int xmin = ifirst(0) - ghosts(0);
        int xmax = ilast(0) + ghosts(0);
        int ymin = ifirst(1) - ghosts(1);
        int ymax = ilast(1) + ghosts(1);

        int xminng = ifirst(0);
        int xmaxng = ilast(0);
        int yminng = ifirst(1);
        int ymaxng = ilast(1);

        int nx = xmax - xmin + 1;
        int ny = ymax - ymin + 1;

        const boost::shared_ptr<geom::CartesianPatchGeometry> pgeom(
                patch.getPatchGeometry(),
                boost::detail::dynamic_cast_tag());

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
        v_volume->fillAll(vol);

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
    }

    int level_number = patch.getPatchLevelNumber();

    boost::shared_ptr<pdat::CellData<int> > level_indicator(patch.getPatchData(
                d_level_indicator,
                getCurrentDataContext()),
            boost::detail::dynamic_cast_tag());

    level_indicator->fillAll(level_number);


    ideal_gas_knl(patch, false);
}

void Cleverleaf::accelerate(
        hier::Patch& patch,
        double dt)
{
    double nodal_mass;

    boost::shared_ptr<pdat::CellData<double> > v_density(
            patch.getPatchData(d_density, getCurrentDataContext()),
            boost::detail::dynamic_cast_tag());

    boost::shared_ptr<pdat::CellData<double> > v_volume(
            patch.getPatchData(d_volume, getCurrentDataContext()),
            boost::detail::dynamic_cast_tag());

    boost::shared_ptr<pdat::CellData<double> > v_celldeltas(
            patch.getPatchData(d_celldeltas, getCurrentDataContext()),
            boost::detail::dynamic_cast_tag());

    boost::shared_ptr<pdat::CellData<double> > v_pressure(
            patch.getPatchData(d_pressure, getCurrentDataContext()),
            boost::detail::dynamic_cast_tag());

    boost::shared_ptr<pdat::CellData<double> > v_viscosity(
            patch.getPatchData(d_viscosity, getCurrentDataContext()),
            boost::detail::dynamic_cast_tag());

    boost::shared_ptr<pdat::NodeData<double> > v_vel0(
            patch.getPatchData(d_velocity, getCurrentDataContext()),
            boost::detail::dynamic_cast_tag());

    boost::shared_ptr<pdat::NodeData<double> > v_vel1(
            patch.getPatchData(d_velocity, getNewDataContext()),
            boost::detail::dynamic_cast_tag());
        
    hier::IntVector ghosts = v_pressure->getGhostCellWidth();

    const hier::Index ifirst = patch.getBox().lower();
    const hier::Index ilast = patch.getBox().upper();

    int xmin = ifirst(0);// - ghosts(0); 
    int xmax = ilast(0);// + ghosts(0); 
    int ymin = ifirst(1);// - ghosts(1); 
    int ymax = ilast(1);// + ghosts(1); 

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

    pdat::NodeData<double> v_stepbymass(patch.getBox(), 1, hier::IntVector(d_dim, 2));
    double* stepbymass = v_stepbymass.getPointer();

    F90_FUNC(accelerate_kernel,ACCELERATE_KERNEL)
        (&xmin,&xmax,&ymin,&ymax,&dt,
         xarea,
         yarea,
         volume,
         density0,
         pressure,
         viscosity,
         xvel0,
         yvel0,
         xvel1,
         yvel1,
         stepbymass);
}

void Cleverleaf::ideal_gas_knl(
        hier::Patch& patch,
        bool predict)
{
    boost::shared_ptr<pdat::CellData<double> > v_density;
    boost::shared_ptr<pdat::CellData<double> > v_energy;

    if (predict) {
        boost::shared_ptr<pdat::CellData<double> > d_tmp(
                patch.getPatchData(d_density, getNewDataContext()),
                boost::detail::dynamic_cast_tag());

        v_density = d_tmp;
    } else {
        boost::shared_ptr<pdat::CellData<double> > d_tmp(
                patch.getPatchData(d_density, getCurrentDataContext()),
                boost::detail::dynamic_cast_tag());

        v_density = d_tmp;
    }

    if (predict) {
        boost::shared_ptr<pdat::CellData<double> > e_tmp(
                patch.getPatchData(d_energy, getNewDataContext()),
                boost::detail::dynamic_cast_tag());

        v_energy = e_tmp;
    } else {
        boost::shared_ptr<pdat::CellData<double> > e_tmp(
                patch.getPatchData(d_energy, getCurrentDataContext()),
                boost::detail::dynamic_cast_tag());

        v_energy = e_tmp;
    }

    boost::shared_ptr<pdat::CellData<double> > v_pressure(
            patch.getPatchData(d_pressure, getCurrentDataContext()),
            boost::detail::dynamic_cast_tag());

    boost::shared_ptr<pdat::CellData<double> > v_soundspeed(
            patch.getPatchData(d_soundspeed, getCurrentDataContext()),
            boost::detail::dynamic_cast_tag());


    hier::IntVector ghosts = v_pressure->getGhostCellWidth();

    const hier::Index ifirst = patch.getBox().lower();
    const hier::Index ilast = patch.getBox().upper();

    int xmin = ifirst(0);// - ghosts(0); 
    int xmax = ilast(0);// + ghosts(0); 
    int ymin = ifirst(1);// - ghosts(1); 
    int ymax = ilast(1);// + ghosts(1); 

    /**
     * nx needs to account for the number of ghosts or things go very wrong!
     */
    int nx = (xmax) - (xmin) + 1;

    //tbox::pout << "xmin: " << nx << std::endl;

    double* density = v_density->getPointer();
    double* energy = v_energy->getPointer();
    double* pressure = v_pressure->getPointer();
    double* soundspeed = v_soundspeed->getPointer();

    F90_FUNC(ideal_gas_kernel,IDEAL_GAS_KERNEL)
        (&xmin, &xmax, &ymin, &ymax,
         density,
         energy,
         pressure,
         soundspeed);
}

void Cleverleaf::viscosity_knl(
        hier::Patch& patch)
{
  boost::shared_ptr<pdat::CellData<double> > v_density0(
          patch.getPatchData(d_density, getCurrentDataContext()),
          boost::detail::dynamic_cast_tag());

  boost::shared_ptr<pdat::CellData<double> > v_pressure(
          patch.getPatchData(d_pressure, getCurrentDataContext()),
          boost::detail::dynamic_cast_tag());

  boost::shared_ptr<pdat::CellData<double> > v_viscosity(
          patch.getPatchData(d_viscosity, getCurrentDataContext()),
          boost::detail::dynamic_cast_tag());

  boost::shared_ptr<pdat::NodeData<double> > v_vel0(
          patch.getPatchData(d_velocity, getCurrentDataContext()),
          boost::detail::dynamic_cast_tag());

  boost::shared_ptr<pdat::CellData<double> > v_celldeltas(
          patch.getPatchData(d_celldeltas, getCurrentDataContext()),
          boost::detail::dynamic_cast_tag());


    hier::IntVector ghosts = v_pressure->getGhostCellWidth();

    const hier::Index ifirst = patch.getBox().lower();
    const hier::Index ilast = patch.getBox().upper();

    int xmin = ifirst(0); 
    int xmax = ilast(0); 
    int ymin = ifirst(1); 
    int ymax = ilast(1); 

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

    F90_FUNC(viscosity_kernel,VISCOSITY_KERNEL)
        (&xmin, &xmax, &ymin, &ymax,
         celldx,
         celldy,
         density,
         pressure,
         viscosity,
         xvel0,
         yvel0); 
}

double Cleverleaf::calc_dt_knl(
        hier::Patch& patch)
{
    boost::shared_ptr<pdat::CellData<double> > v_density1(
            patch.getPatchData(d_density, getNewDataContext()),
            boost::detail::dynamic_cast_tag());
    boost::shared_ptr<pdat::CellData<double> > v_density0(
            patch.getPatchData(d_density, getCurrentDataContext()),
            boost::detail::dynamic_cast_tag());

    boost::shared_ptr<pdat::CellData<double> > v_energy0(
            patch.getPatchData(d_energy, getCurrentDataContext()),
            boost::detail::dynamic_cast_tag());
    boost::shared_ptr<pdat::CellData<double> > v_energy1(
            patch.getPatchData(d_energy, getNewDataContext()),
            boost::detail::dynamic_cast_tag());

    boost::shared_ptr<pdat::CellData<double> > v_pressure(
            patch.getPatchData(d_pressure, getCurrentDataContext()),
            boost::detail::dynamic_cast_tag());

    boost::shared_ptr<pdat::CellData<double> > v_soundspeed(
            patch.getPatchData(d_soundspeed, getCurrentDataContext()),
            boost::detail::dynamic_cast_tag());

    boost::shared_ptr<pdat::CellData<double> > v_viscosity(
            patch.getPatchData(d_viscosity, getCurrentDataContext()),
            boost::detail::dynamic_cast_tag());

    boost::shared_ptr<pdat::NodeData<double> > v_vel0(
            patch.getPatchData(d_velocity, getCurrentDataContext()),
            boost::detail::dynamic_cast_tag());

    boost::shared_ptr<pdat::CellData<double> > v_celldeltas(
            patch.getPatchData(d_celldeltas, getCurrentDataContext()),
            boost::detail::dynamic_cast_tag());

    boost::shared_ptr<pdat::CellData<double> > v_volume(
            patch.getPatchData(d_volume, getCurrentDataContext()),
            boost::detail::dynamic_cast_tag());

    boost::shared_ptr<pdat::CellData<double> > v_cellcoords(
            patch.getPatchData(d_cellcoords, getCurrentDataContext()),
            boost::detail::dynamic_cast_tag());

    boost::shared_ptr<pdat::CellData<double> > v_energy(
            patch.getPatchData(d_energy, getCurrentDataContext()),
            boost::detail::dynamic_cast_tag());

    hier::IntVector ghosts = v_pressure->getGhostCellWidth();

    const hier::Index ifirst = patch.getBox().lower();
    const hier::Index ilast = patch.getBox().upper();

    int xmin = ifirst(0);
    int xmax = ilast(0);
    int ymin = ifirst(1);
    int ymax = ilast(1);

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
    double* energy0 = v_energy->getPointer();
    double* xarea = v_celldeltas->getPointer(1);
    double* yarea = v_celldeltas->getPointer(0);
    double* volume = v_volume->getPointer();
    double* cellx = v_cellcoords->getPointer(0);
    double* celly = v_cellcoords->getPointer(1);

    double dt_min_val = 1.0e+21;
    int small=0;
    double dtc_safe = 0.7;
    int dtl_control;
    double dtu_safe = 0.5;
    double dtv_safe = 0.5;
    double dtdiv_safe = 0.7;
    double xl_pos,
           yl_pos;
    double dt_min = 0.0000001;

    int kldt;
    int jldt;

    double g_small = 1.0e-16;
    double g_big = 1.0e+21;

    pdat::NodeData<double> v_dtmin(patch.getBox(), 1, hier::IntVector(d_dim, 2));
    double* dtmin = v_dtmin.getPointer();

    F90_FUNC(calc_dt_kernel, CALC_DT_KERNEL)
        (&xmin, &xmax, &ymin, &ymax,
         &g_small,
         &g_big,
         &dt_min,
         &dtc_safe,
         &dtu_safe,
         &dtv_safe,
         &dtdiv_safe,
         xarea,
         yarea,
         cellx,
         celly,
         celldx,
         celldy,
         volume,
         density0,
         energy0,
         pressure,
         viscosity,
         soundspeed,
         xvel0,
         yvel0,
         dtmin,
         &dt_min_val,
         &dtl_control,
         &xl_pos,
         &yl_pos,
         &jldt,
         &kldt,
         &small);

    return dt_min_val;
}

void Cleverleaf::pdv_knl(
        hier::Patch& patch,
        double dt,
        bool predict)
{
    /*
     * Get necessary variables
     */
    boost::shared_ptr<pdat::CellData<double> > area(
            patch.getPatchData(d_celldeltas, getCurrentDataContext()),
            boost::detail::dynamic_cast_tag());
    boost::shared_ptr<pdat::CellData<double> > vol(
            patch.getPatchData(d_volume, getCurrentDataContext()),
            boost::detail::dynamic_cast_tag());
    boost::shared_ptr<pdat::CellData<double> > dens0(
            patch.getPatchData(d_density, getCurrentDataContext()),
            boost::detail::dynamic_cast_tag());
    boost::shared_ptr<pdat::CellData<double> > dens1(
            patch.getPatchData(d_density, getNewDataContext()),
            boost::detail::dynamic_cast_tag());
    boost::shared_ptr<pdat::CellData<double> > en0(
            patch.getPatchData(d_energy, getCurrentDataContext()),
            boost::detail::dynamic_cast_tag());
    boost::shared_ptr<pdat::CellData<double> > en1(
            patch.getPatchData(d_energy, getNewDataContext()),
            boost::detail::dynamic_cast_tag());
    boost::shared_ptr<pdat::CellData<double> > pres(
            patch.getPatchData(d_pressure, getCurrentDataContext()),
            boost::detail::dynamic_cast_tag());
    boost::shared_ptr<pdat::CellData<double> > visc(
            patch.getPatchData(d_viscosity, getCurrentDataContext()),
            boost::detail::dynamic_cast_tag());
    boost::shared_ptr<pdat::NodeData<double> > v0(
            patch.getPatchData(d_velocity, getCurrentDataContext()),
            boost::detail::dynamic_cast_tag());
    boost::shared_ptr<pdat::NodeData<double> > v1(
            patch.getPatchData(d_velocity, getNewDataContext()),
            boost::detail::dynamic_cast_tag());

    hier::IntVector ghosts = pres->getGhostCellWidth();

    const hier::Index ifirst = patch.getBox().lower();
    const hier::Index ilast = patch.getBox().upper();

    int xmin = ifirst(0);
    int xmax = ilast(0);
    int ymin = ifirst(1);
    int ymax = ilast(1);

    int prdct = -1;

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

    pdat::NodeData<double> v_volchange(patch.getBox(), 1, hier::IntVector(d_dim, 2));
    double* volume_change = v_volchange.getPointer();

    if (predict)
        prdct = 1;
    else
        prdct = 0;

    F90_FUNC(pdv_kernel, PDV_KERNEL)
        (&prdct,&xmin,&xmax,&ymin,&ymax,&dt,
         xarea,
         yarea,
         volume,
         density0,
         density1,
         energy0,
         energy1,
         pressure,
         viscosity,
         xvel0,
         xvel1,
         yvel0,
         yvel1,
         volume_change);
}

void Cleverleaf::flux_calc_knl(
        hier::Patch& patch,
        double dt)
{

    boost::shared_ptr<pdat::NodeData<double> > v_vel0(
            patch.getPatchData(d_velocity, getCurrentDataContext()),
            boost::detail::dynamic_cast_tag());
    boost::shared_ptr<pdat::NodeData<double> > v_vel1(
            patch.getPatchData(d_velocity, getNewDataContext()),
            boost::detail::dynamic_cast_tag());

    boost::shared_ptr<pdat::EdgeData<double> > v_volflux(
            patch.getPatchData(d_volflux, getCurrentDataContext()),
            boost::detail::dynamic_cast_tag());

    boost::shared_ptr<pdat::CellData<double> > v_celldeltas(
            patch.getPatchData(d_celldeltas, getCurrentDataContext()),
            boost::detail::dynamic_cast_tag());

    hier::IntVector ghosts = v_celldeltas->getGhostCellWidth();

    const hier::Index ifirst = patch.getBox().lower();
    const hier::Index ilast = patch.getBox().upper();

    int xmin = ifirst(0); 
    int xmax = ilast(0); 
    int ymin = ifirst(1); 
    int ymax = ilast(1); 

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

    F90_FUNC(flux_calc_kernel, FLUX_CALC_KERNEL)
        (&xmin, &xmax, &ymin, &ymax,
         &dt,
         xarea,
         yarea,
         xvel0,
         yvel0,
         xvel1,
         yvel1,
         vol_flux_x,
         vol_flux_y);
}

void Cleverleaf::advec_cell(hier::Patch& patch,
        int sweep_number,
        ADVEC_DIR dir)
{
    boost::shared_ptr<pdat::CellData<double> > v_density1(
            patch.getPatchData(d_density, getNewDataContext()),
            boost::detail::dynamic_cast_tag());

    boost::shared_ptr<pdat::CellData<double> > v_energy1(
            patch.getPatchData(d_energy, getNewDataContext()),
            boost::detail::dynamic_cast_tag());

    boost::shared_ptr<pdat::CellData<double> > v_volume(
            patch.getPatchData(d_volume, getCurrentDataContext()),
            boost::detail::dynamic_cast_tag());

    boost::shared_ptr<pdat::EdgeData<double> > v_volflux(
            patch.getPatchData(d_volflux, getCurrentDataContext()),
            boost::detail::dynamic_cast_tag());

    boost::shared_ptr<pdat::EdgeData<double> > v_massflux(
            patch.getPatchData(d_massflux, getCurrentDataContext()),
            boost::detail::dynamic_cast_tag());

    boost::shared_ptr<pdat::NodeData<double> > vertexdeltas(
            patch.getPatchData(d_vertexdeltas, getCurrentDataContext()),
            boost::detail::dynamic_cast_tag());

    hier::IntVector ghosts = v_density1->getGhostCellWidth();

    const hier::Index ifirst = patch.getBox().lower();
    const hier::Index ilast = patch.getBox().upper();

    int xmin = ifirst(0); 
    int xmax = ilast(0); 
    int ymin = ifirst(1); 
    int ymax = ilast(1); 

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

    pdat::NodeData<double> v_prevol(patch.getBox(), 1, hier::IntVector(d_dim, 2));
    pdat::NodeData<double> v_postvol(patch.getBox(), 1, hier::IntVector(d_dim, 2));
    pdat::NodeData<double> v_premass(patch.getBox(), 1, hier::IntVector(d_dim, 2));
    pdat::NodeData<double> v_postmass(patch.getBox(), 1, hier::IntVector(d_dim, 2));
    pdat::NodeData<double> v_advecvol(patch.getBox(), 1, hier::IntVector(d_dim, 2));
    pdat::NodeData<double> v_postener(patch.getBox(), 1, hier::IntVector(d_dim, 2));
    pdat::NodeData<double> v_enerflux(patch.getBox(), 1, hier::IntVector(d_dim, 2));

    double* pre_vol = v_prevol.getPointer();
    double* post_vol = v_postvol.getPointer();
    double* pre_mass = v_premass.getPointer();
    double* post_mass = v_postmass.getPointer();
    double* advec_vol = v_advecvol.getPointer();
    double* post_ener = v_postener.getPointer();
    double* ener_flux = v_enerflux.getPointer();

    int idir;

    if(dir == X)
        idir = 1;
    else idir = 2;

    F90_FUNC(advec_cell_kernel, ADVEC_CELL_KERNEL)
        (&xmin, &xmax, &ymin, &ymax, &idir, &sweep_number,
         vertexdx,
         vertexdy,
         volume,
         density1,
         energy1,
         mass_flux_x,
         vol_flux_x,
         mass_flux_y,
         vol_flux_y,
         pre_vol,
         post_vol,
         pre_mass,
         post_mass,
         advec_vol,
         post_ener,
         ener_flux);
}

void Cleverleaf::advec_mom(hier::Patch& patch,
        int sweep_number,
        ADVEC_DIR direction,
        ADVEC_DIR which_vel)
{
    boost::shared_ptr<pdat::CellData<double> > v_density1(
            patch.getPatchData(d_density, getNewDataContext()),
            boost::detail::dynamic_cast_tag());

    boost::shared_ptr<pdat::CellData<double> > v_volume(
            patch.getPatchData(d_volume, getCurrentDataContext()),
            boost::detail::dynamic_cast_tag());

    boost::shared_ptr<pdat::NodeData<double> > v_vel1(
            patch.getPatchData(d_velocity, getNewDataContext()),
            boost::detail::dynamic_cast_tag());

    boost::shared_ptr<pdat::EdgeData<double> > v_volflux(
            patch.getPatchData(d_volflux, getCurrentDataContext()),
            boost::detail::dynamic_cast_tag());

    boost::shared_ptr<pdat::EdgeData<double> > v_massflux(
            patch.getPatchData(d_massflux, getCurrentDataContext()),
            boost::detail::dynamic_cast_tag());

    boost::shared_ptr<pdat::CellData<double> > v_celldeltas(
            patch.getPatchData(d_celldeltas, getCurrentDataContext()),
            boost::detail::dynamic_cast_tag());

    hier::IntVector ghosts = v_density1->getGhostCellWidth();

    const hier::Index ifirst = patch.getBox().lower();
    const hier::Index ilast = patch.getBox().upper();

    int xmin = ifirst(0); 
    int xmax = ilast(0); 
    int ymin = ifirst(1); 
    int ymax = ilast(1); 

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

    double* xvel1 = v_vel1->getPointer(0);
    double* yvel1 = v_vel1->getPointer(1);

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

    pdat::NodeData<double> v_prevol(patch.getBox(), 1, hier::IntVector(d_dim, 2));
    pdat::NodeData<double> v_postvol(patch.getBox(), 1, hier::IntVector(d_dim, 2));

    double* node_flux = v_nodeflux.getPointer();
    double* node_mass_post = v_nodemasspost.getPointer();
    double* node_mass_pre = v_nodemasspre.getPointer();
    double* advec_vel = v_advecvel.getPointer();
    double* mom_flux = v_momflux.getPointer();
    double* pre_vol = v_prevol.getPointer();
    double* post_vol = v_postvol.getPointer();

    int iwhich, idir;

    if (which_vel == X)
        iwhich = 1;
    else iwhich = 2;

    if (direction == X)
        idir = 1;
    else idir = 2;

    F90_FUNC(advec_mom_kernel, ADVEC_MOM_KERNEL)
        (&xmin, &xmax, &ymin, &ymax,
         xvel1,
         yvel1,
         mass_flux_x,
         vol_flux_x,
         mass_flux_y,
         vol_flux_y,
         volume,
         density1,
         node_flux,
         node_mass_post,
         node_mass_pre,
         advec_vel,
         mom_flux,
         pre_vol,
         post_vol,
         celldx,
         celldy,
         &iwhich,
         &sweep_number,
         &idir);
}


void Cleverleaf::setPhysicalBoundaryConditions(
        hier::Patch& patch,
        const double fill_time,
        const hier::IntVector& ghost_width_to_fill)
{

    boost::shared_ptr<pdat::CellData<double> > v_pressure(
            patch.getPatchData(d_pressure, getScratchDataContext()),
            boost::detail::dynamic_cast_tag());

    boost::shared_ptr<pdat::CellData<double> > v_density0(
            patch.getPatchData(d_density, getScratchDataContext()),
            boost::detail::dynamic_cast_tag());

    boost::shared_ptr<pdat::CellData<double> > v_density1(
            patch.getPatchData(d_density, getScratchNewDataContext()),
            boost::detail::dynamic_cast_tag());

    boost::shared_ptr<pdat::CellData<double> > v_energy0(
            patch.getPatchData(d_energy, getScratchDataContext()),
            boost::detail::dynamic_cast_tag());

    boost::shared_ptr<pdat::CellData<double> > v_energy1(
            patch.getPatchData(d_energy, getScratchNewDataContext()),
            boost::detail::dynamic_cast_tag());

    boost::shared_ptr<pdat::CellData<double> > v_viscosity(
            patch.getPatchData(d_viscosity, getScratchDataContext()),
            boost::detail::dynamic_cast_tag());

    boost::shared_ptr<pdat::NodeData<double> > v_vel0(
            patch.getPatchData(d_velocity, getScratchDataContext()),
            boost::detail::dynamic_cast_tag());

    boost::shared_ptr<pdat::NodeData<double> > v_vel1(
            patch.getPatchData(d_velocity, getScratchNewDataContext()),
            boost::detail::dynamic_cast_tag());

    boost::shared_ptr<pdat::EdgeData<double> > v_massflux( 
            patch.getPatchData(d_massflux, getScratchDataContext()),
            boost::detail::dynamic_cast_tag());

    boost::shared_ptr<pdat::EdgeData<double> > v_volflux( 
            patch.getPatchData(d_volflux, getScratchDataContext()),
            boost::detail::dynamic_cast_tag());

    boost::shared_ptr<pdat::CellData<double> > v_soundspeed( 
            patch.getPatchData(d_soundspeed, getScratchDataContext()),
            boost::detail::dynamic_cast_tag());

    const hier::Index ifirst = patch.getBox().lower();
    const hier::Index ilast = patch.getBox().upper();

    int xmin = ifirst(0); 
    int xmax = ilast(0); 
    int ymin = ifirst(1); 
    int ymax = ilast(1); 

    int nx = xmax - xmin + 1;

    double* pressure;

    double* density0;
    double* density1;

    double* energy0;
    double* energy1;

    double* viscosity;

    double* xvel0;
    double* xvel1;

    double* yvel0;
    double* yvel1;

    double* mass_flux_x;
    double* mass_flux_y;

    double* vol_flux_x;
    double* vol_flux_y;

    double* soundspeed;

    int depth = ghost_width_to_fill[0];

    const boost::shared_ptr<geom::CartesianPatchGeometry> pgeom(
            patch.getPatchGeometry(),
            boost::detail::dynamic_cast_tag());

    const std::vector<hier::BoundaryBox>& edge_bdry = pgeom->getCodimensionBoundaries(Bdry::EDGE2D);

    for(int i = 0; i < 15; i++) {
        d_exchange_fields[i] = 0;
    }

    switch(d_which_exchange) {
        case LagrangianEulerianLevelIntegrator::FIELD_EXCH:
            {
                d_exchange_fields[FIELD_DENSITY0] = 1;
                d_exchange_fields[FIELD_ENERGY0] = 1;
                d_exchange_fields[FIELD_XVEL0] = 1;
                d_exchange_fields[FIELD_YVEL0] = 1;

                density0 = v_density0->getPointer();
                energy0 = v_energy0->getPointer();
                xvel0 = v_vel0->getPointer(0);
                yvel0 = v_vel0->getPointer(1);
            } break;
        case LagrangianEulerianLevelIntegrator::PRIME_CELLS_EXCH: 
            {
                d_exchange_fields[FIELD_DENSITY0] = 1;
                d_exchange_fields[FIELD_DENSITY1] = 1;
                d_exchange_fields[FIELD_ENERGY0] = 1;
                d_exchange_fields[FIELD_ENERGY1] = 1;
                d_exchange_fields[FIELD_PRESSURE] = 1;
                d_exchange_fields[FIELD_VISCOSITY] = 1;
                d_exchange_fields[FIELD_XVEL0] = 1;
                d_exchange_fields[FIELD_XVEL1] = 1;
                d_exchange_fields[FIELD_YVEL0] = 1;
                d_exchange_fields[FIELD_YVEL1] = 1;

                pressure = v_pressure->getPointer();

                density0 = v_density0->getPointer();
                density1 = v_density1->getPointer();

                energy0 = v_energy0->getPointer();
                energy1 = v_energy1->getPointer();

                viscosity = v_viscosity->getPointer();

                xvel0 = v_vel0->getPointer(0);
                xvel1 = v_vel1->getPointer(0);

                yvel0 = v_vel0->getPointer(1);
                yvel1 = v_vel1->getPointer(1);
            } break;
        case LagrangianEulerianLevelIntegrator::PRE_LAGRANGE_EXCH:
            {
                d_exchange_fields[FIELD_XVEL0] = 1;
                d_exchange_fields[FIELD_YVEL0] = 1;
                d_exchange_fields[FIELD_PRESSURE] = 1;
                d_exchange_fields[FIELD_DENSITY0] = 1;
                d_exchange_fields[FIELD_ENERGY0] = 1;

                pressure = v_pressure->getPointer();
                density0 = v_density0->getPointer();
                energy0 = v_energy0->getPointer();
                xvel0 = v_vel0->getPointer(0);
                yvel0 = v_vel0->getPointer(1);
            } break;
        case LagrangianEulerianLevelIntegrator::POST_VISCOSITY_EXCH:
            {
                d_exchange_fields[FIELD_VISCOSITY] = 1;

                viscosity = v_viscosity->getPointer();
            } break;
        case LagrangianEulerianLevelIntegrator::HALF_STEP_EXCH:
            {
                d_exchange_fields[FIELD_PRESSURE] = 1;

                pressure = v_pressure->getPointer();
            } break;
        case LagrangianEulerianLevelIntegrator::PRE_SWEEP_1_CELL_EXCH:
            {
                d_exchange_fields[FIELD_VOL_FLUX_X] = 1;
                d_exchange_fields[FIELD_VOL_FLUX_Y] = 1;
                d_exchange_fields[FIELD_DENSITY1] = 1;
                d_exchange_fields[FIELD_ENERGY1] = 1;

                density1 = v_density1->getPointer();
                energy1 = v_energy1->getPointer();
                vol_flux_x = v_volflux->getPointer(1);
                vol_flux_y = v_volflux->getPointer(0);
            } break;
        case LagrangianEulerianLevelIntegrator::PRE_SWEEP_1_MOM_EXCH:
            {
                d_exchange_fields[FIELD_XVEL1] = 1;
                d_exchange_fields[FIELD_YVEL1] = 1;
                d_exchange_fields[FIELD_DENSITY1] = 1;
                d_exchange_fields[FIELD_ENERGY1] = 1;
                d_exchange_fields[FIELD_MASS_FLUX_X] = 1;
                d_exchange_fields[FIELD_MASS_FLUX_Y] = 1;

                density1 = v_density1->getPointer();
                energy1 = v_energy1->getPointer();

                xvel1 = v_vel1->getPointer(0);
                yvel1 = v_vel1->getPointer(1);

                mass_flux_x = v_massflux->getPointer(1);
                mass_flux_y = v_massflux->getPointer(0);
            } break;
        case LagrangianEulerianLevelIntegrator::PRE_SWEEP_2_MOM_EXCH:
            {
                d_exchange_fields[FIELD_XVEL1] = 1;
                d_exchange_fields[FIELD_YVEL1] = 1;
                d_exchange_fields[FIELD_DENSITY1] = 1;
                d_exchange_fields[FIELD_ENERGY1] = 1;
                d_exchange_fields[FIELD_MASS_FLUX_X] = 1;
                d_exchange_fields[FIELD_MASS_FLUX_Y] = 1;

                density1 = v_density1->getPointer();
                energy1 = v_energy1->getPointer();

                xvel1 = v_vel1->getPointer(0);
                yvel1 = v_vel1->getPointer(1);

                mass_flux_x = v_massflux->getPointer(1);
                mass_flux_y = v_massflux->getPointer(0);
            } break;
        default : tbox::perr << "[ERROR] Unknown exchange id in setPhysicalBoundaryConditions... " << std::endl;
                  exit(-1);
    }

    for(int i = 0; i < edge_bdry.size(); i++) {
        switch(edge_bdry[i].getLocationIndex()) {
            case (BdryLoc::YLO) :
                F90_FUNC(update_halo_kernel_bottom, UPDATE_HALO_KERNEL_BOTTOM)
                    (&xmin, &xmax, &ymin, &ymax,
                     density0,
                     energy0,
                     pressure,
                     viscosity,
                     soundspeed,
                     density1,
                     energy1,
                     xvel0,
                     yvel0,
                     xvel1,
                     yvel1,
                     vol_flux_x,
                     vol_flux_y,
                     mass_flux_x,
                     mass_flux_y,
                     d_exchange_fields,
                     &depth);
                break;


            case (BdryLoc::YHI) :
                F90_FUNC(update_halo_kernel_top, UPDATE_HALO_KERNEL_TOP)
                    (&xmin, &xmax, &ymin, &ymax,
                     density0,
                     energy0,
                     pressure,
                     viscosity,
                     soundspeed,
                     density1,
                     energy1,
                     xvel0,
                     yvel0,
                     xvel1,
                     yvel1,
                     vol_flux_x,
                     vol_flux_y,
                     mass_flux_x,
                     mass_flux_y,
                     d_exchange_fields,
                     &depth);
                break;


            case (BdryLoc::XLO) :
                F90_FUNC(update_halo_kernel_left, UPDATE_HALO_KERNEL_LEFT)
                    (&xmin, &xmax, &ymin, &ymax,
                     density0,
                     energy0,
                     pressure,
                     viscosity,
                     soundspeed,
                     density1,
                     energy1,
                     xvel0,
                     yvel0,
                     xvel1,
                     yvel1,
                     vol_flux_x,
                     vol_flux_y,
                     mass_flux_x,
                     mass_flux_y,
                     d_exchange_fields,
                     &depth);
                break;

            case (BdryLoc::XHI) :
                F90_FUNC(update_halo_kernel_right, UPDATE_HALO_KERNEL_RIGHT)
                    (&xmin, &xmax, &ymin, &ymax,
                     density0,
                     energy0,
                     pressure,
                     viscosity,
                     soundspeed,
                     density1,
                     energy1,
                     xvel0,
                     yvel0,
                     xvel1,
                     yvel1,
                     vol_flux_x,
                     vol_flux_y,
                     mass_flux_x,
                     mass_flux_y,
                     d_exchange_fields,
                     &depth);
                break;
            default : tbox::perr << "[ERROR] Unknown edge location in setPhysicalBoundaryConditions... " << std::endl;
                      exit(-1);
        }
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
    boost::shared_ptr<pdat::CellData<double> > v_volume(
            patch.getPatchData(d_volume, getCurrentDataContext()),
            boost::detail::dynamic_cast_tag());

    boost::shared_ptr<pdat::CellData<double> > v_density0(
            patch.getPatchData(d_density, getCurrentDataContext()),
            boost::detail::dynamic_cast_tag());

    boost::shared_ptr<pdat::CellData<double> > v_energy0(
            patch.getPatchData(d_energy, getCurrentDataContext()),
            boost::detail::dynamic_cast_tag());

    boost::shared_ptr<pdat::CellData<double> > v_pressure(
            patch.getPatchData(d_pressure, getCurrentDataContext()),
            boost::detail::dynamic_cast_tag());

    boost::shared_ptr<pdat::NodeData<double> > v_vel0(
            patch.getPatchData(d_velocity, getCurrentDataContext()),
            boost::detail::dynamic_cast_tag());

    boost::shared_ptr<pdat::CellData<int> > v_level_indicator(
            patch.getPatchData(d_level_indicator, getCurrentDataContext()),
            boost::detail::dynamic_cast_tag());

    const hier::Index ifirst = patch.getBox().lower();
    const hier::Index ilast = patch.getBox().upper();

    int xmin = ifirst(0);
    int xmax = ilast(0);
    int ymin = ifirst(1);
    int ymax = ilast(1);

    double* density0 = v_density0->getPointer();
    double* energy0 = v_energy0->getPointer();
    double* pressure = v_pressure->getPointer();
    double* volume = v_volume->getPointer();
    double* xvel0 = v_vel0->getPointer(0);
    double* yvel0 = v_vel0->getPointer(1);

    int* level_indicator = v_level_indicator->getPointer();

    int level_number = patch.getPatchLevelNumber();

    *vol=0.0;
    *mass=0.0;
    *ie=0.0;
    *ke=0.0;
    *press=0.0;

    F90_FUNC(field_summary_kernel, FIELD_SUMMARY_KERNEL)
        (&xmin, &xmax, &ymin, &ymax,
         volume,
         density0,
         energy0,
         pressure,
         xvel0,
         yvel0,
         level_indicator,
         vol,
         mass,
         ie,
         ke,
         press,
         &level_number);

    std::cout << "values: " << *vol << ", " << *mass << std::endl;
}

void Cleverleaf::tagGradientDetectorCells(
        hier::Patch& patch,
        const double regrid_time,
        const bool initial_error,
        const int tag_index)
{

    const int error_level_number = patch.getPatchLevelNumber();

    const boost::shared_ptr<geom::CartesianPatchGeometry> patch_geom(
            patch.getPatchGeometry(),
            boost::detail::dynamic_cast_tag());

    const double* dx = patch_geom->getDx();

    boost::shared_ptr<pdat::CellData<int> > tags(
            patch.getPatchData(tag_index),
            boost::detail::dynamic_cast_tag());

    boost::shared_ptr<pdat::CellData<double> > v_density(
            patch.getPatchData(d_density, getCurrentDataContext()),
            boost::detail::dynamic_cast_tag());

    boost::shared_ptr<pdat::CellData<double> > v_energy(
            patch.getPatchData(d_energy, getCurrentDataContext()),
            boost::detail::dynamic_cast_tag());

    boost::shared_ptr<pdat::CellData<double> > v_viscosity(
            patch.getPatchData(d_viscosity, getCurrentDataContext()),
            boost::detail::dynamic_cast_tag());

    hier::Box pbox = patch.getBox();
    hier::BoxContainer domain_boxes;
    d_grid_geometry->computePhysicalDomain(domain_boxes, patch_geom->getRatio(), hier::BlockId::zero());

    // IF xvel > 0.1 tag a cell!

    /*
     * Construct domain bounding box
     */
    hier::Box domain(d_dim);
    for (hier::BoxContainer::iterator ib = domain_boxes.begin(); ib != domain_boxes.end(); ++ib) {
        domain += *ib;
    }

    hier::Index domfirst = domain.lower();
    hier::Index domlast = domain.upper();
    hier::Index ifirst = patch.getBox().lower();
    hier::Index ilast = patch.getBox().upper();

    hier::IntVector nghosts = v_density->getGhostCellWidth();
    hier::IntVector tnghosts = tags->getGhostCellWidth();

    int xmin = ifirst(0) - nghosts[0];
    int ymin = ifirst(1) - nghosts[1];
    int xmax = ilast(0) + nghosts[0];
    int ymax = ilast(1) + nghosts[1];

    int tx_min = ifirst(0) - tnghosts[0];
    int ty_min = ifirst(1) - tnghosts[1];
    int tx_max = ilast(0) + tnghosts[0];
    int ty_max = ilast(1) + tnghosts[1];

    int nx = (xmax - xmin) + 1;
    int tnx = (tx_max - tx_min) + 1;

    /*
     * Create a set of temporary tags and set to untagged value.
     */

    boost::shared_ptr<pdat::CellData<int> > temp_tags(new pdat::CellData<int>(
                pbox,
                1,
                nghosts));

    temp_tags->fillAll(0);
    tags->fillAll(0);

    double* density0 = v_density->getPointer();
    double* energy0 = v_energy->getPointer();
    double* viscosity = v_viscosity->getPointer();
    int* temp_tags_array = temp_tags->getPointer();

    F90_FUNC(tag_q_kernel,TAG_Q_KERNEL)
        (&ifirst(0), &ilast(0), &ifirst(1), &ilast(1),
         viscosity,
         temp_tags_array);

    F90_FUNC(tag_density_kernel,TAG_DENSITY_KERNEL)
        (&ifirst(0), &ilast(0), &ifirst(1), &ilast(1),
         density0,
         temp_tags_array);

    F90_FUNC(tag_energy_kernel,TAG_ENERGY_KERNEL)
        (&ifirst(0), &ilast(0), &ifirst(1), &ilast(1),
         energy0,
         temp_tags_array);

    if(d_tag_all) {
        F90_FUNC(tag_all_kernel,TAG_ALL_KERNEL)
            (&ifirst(0), &ilast(0), &ifirst(1), &ilast(1),
             temp_tags_array);
    }

    /*
     * Update tags
     */
   pdat::CellIterator icend(pdat::CellGeometry::end(pbox));

   for (pdat::CellIterator ic(pdat::CellGeometry::begin(pbox)); ic != icend; ++ic) {
      (*tags)(*ic, 0) = (*temp_tags)(*ic, 0);
   }
}

void Cleverleaf::debug_knl(hier::Patch& patch)
{
    boost::shared_ptr<pdat::NodeData<double> > v_vel0(
            patch.getPatchData(d_velocity, getCurrentDataContext()),
            boost::detail::dynamic_cast_tag());

    boost::shared_ptr<pdat::NodeData<double> > v_vel1(
            patch.getPatchData(d_velocity, getNewDataContext()),
            boost::detail::dynamic_cast_tag());

    const hier::Index ifirst = patch.getBox().lower();
    const hier::Index ilast = patch.getBox().upper();

    int xmin = ifirst(0); 
    int xmax = ilast(0); 
    int ymin = ifirst(1); 
    int ymax = ilast(1); 

    double* xvel0 = v_vel0->getPointer(0);
    double* xvel1 = v_vel1->getPointer(0);

    double* yvel0 = v_vel0->getPointer(1);
    double* yvel1 = v_vel1->getPointer(1);

    F90_FUNC(debug_kernel, DEBUG_KERNEL)
        (&xmin, &xmax, &ymin, &ymax,
         xvel0,
         yvel0,
         xvel1,
         yvel1);
}
