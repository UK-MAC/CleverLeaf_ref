//////////////////////////////////////////////////////////////////////////////
// Crown Copyright 2014 AWE, Copyright 2014 David Beckingsale.
//
// This file is part of CleverLeaf.
//
// CleverLeaf is free software: you can redistribute it and/or modify it under
// the terms of the GNU General Public License as published by the Free Software
// Foundation, either version 3 of the License, or (at your option) any later
// version.
//
// CleverLeaf is distributed in the hope that it will be useful, but WITHOUT ANY
// WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR
// A PARTICULAR PURPOSE. See the GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License along with
// CleverLeaf. If not, see http://www.gnu.org/licenses/.
//////////////////////////////////////////////////////////////////////////////
#include "Cleverleaf.h"

#include <iostream>
#include <cmath>

#include "SAMRAI/hier/VariableDatabase.h"
#include "SAMRAI/geom/CartesianPatchGeometry.h"

#include "SAMRAI/pdat/SideDoubleConstantRefine.h"

//#include "SAMRAI/appu/VisItDataWriter.h"

#include "geom/CartesianCleverCellDoubleVolumeWeightedAverage.h"
#include "geom/CartesianCleverCellDoubleMassWeightedAverage.h"
#include "geom/CartesianCleverCellIntConstantCoarsen.h"
#include "geom/CartesianCleverSideDoubleFirstOrderRefine.h"

#include "geom/CartesianCleverCellDoubleConservativeLinearRefine.h"
#include "geom/CartesianCleverNodeDoubleLinearRefine.h"

#include "pdat/CleverNodeInjection.h"

#include "macros.h"

#define F90_FUNC(name,NAME) name ## _

extern "C" {
  void F90_FUNC(ideal_gas_kernel,IDEAL_GAS_KERNEL)
    (int*,int*,int*,int*,double*,double*,double*,double*);

  void F90_FUNC(accelerate_kernel, ACCELERATE_KERNEL)
    (int*,int*,int*,int*,double*,int*,double*,double*,double*,double*,double*,
     double*,double*,double*,double*,double*,double*);

  void F90_FUNC(viscosity_kernel, VISCOSITY_KERNEL)
    (int*,int*,int*,int*,double*,double*,double*,double*,double*,double*,
     double*);

  void F90_FUNC(pdv_kernel, PDV_KERNEL)
    (int*,int*,int*,int*,int*,double*,double*,double*,double*,double*,double*,
     double*,double*,double*,double*,double*,double*,double*,double*,double*);

  void F90_FUNC(advec_cell_kernel, ADVEC_CELL_KERNEL)
    (int*,int*,int*,int*,int*,int*,double*,double*,double*,double*,double*,
     double*,double*,double*,double*,double*,double*,double*,double*,double*,
     double*,double*);

  void F90_FUNC(advec_mom_kernel, ADVEC_MOM_KERNEL)
    (int*,int*,int*,int*,double*,double*,double*,double*,double*,double*,
     double*,double*,double*,double*,double*,double*,double*,double*,double*,
     double*,double*,int*,int*,int*);

  void F90_FUNC(calc_dt_kernel, CALC_DT_KERNEL)
    (int*,int*,int*,int*,const double&,const double&,double*,double*,double*,
     double*,double*,double*,double*,double*,double*,double*,double*,double*,
     double*,double*,double*,double*,double*,double*,double*,double*,double*,
     int*,double*,double*,int*,int*,int*);

  void F90_FUNC(flux_calc_kernel, FLUX_CALC_KERNEL)
    (int*,int*,int*,int*,double*,double*,double*,double*,double*,double*,
     double*,double*,double*);

  void F90_FUNC(update_halo_kernel_top, UPDATE_HALO_KERNEL_TOP)
    (int*,int*,int*,int*,double*,double*,double*,double*,double*,double*,
     double*,double*,double*,double*,double*,double*,double*,double*,double*,
     int*,int*);

  void F90_FUNC(update_halo_kernel_bottom, UPDATE_HALO_KERNEL_BOTTOM)
    (int*,int*,int*,int*,double*,double*,double*,double*,double*,double*,
     double*,double*,double*,double*,double*,double*,double*,double*,double*,
     int*,int*);

  void F90_FUNC(update_halo_kernel_left, UPDATE_HALO_KERNEL_LEFT)
    (int*,int*,int*,int*,double*,double*,double*,double*,double*,double*,
     double*,double*,double*,double*,double*,double*,double*,double*,double*,
     int*,int*);

  void F90_FUNC(update_halo_kernel_right, UPDATE_HALO_KERNEL_RIGHT)
    (int*,int*,int*,int*,double*,double*,double*,double*,double*,double*,
     double*,double*,double*,double*,double*,double*,double*,double*,double*,
     int*,int*);

  void F90_FUNC(tag_q_kernel,TAG_Q_KERNEL)
    (int*,int*,int*,int*,double*,double*,int*);

  void F90_FUNC(tag_energy_kernel,TAG_ENERGY_KERNEL)
    (int*,int*,int*,int*,double*,double*,int*);

  void F90_FUNC(tag_density_kernel,TAG_DENSITY_KERNEL)
    (int*,int*,int*,int*,double*,double*,int*);

  void F90_FUNC(tag_pressure_kernel,TAG_PRESSURE_KERNEL)
    (int*,int*,int*,int*,double*,double*,int*);

  void F90_FUNC(tag_all_kernel,TAG_ALL_KERNEL)
    (int*,int*,int*,int*,int*);

  void F90_FUNC(debug_kernel, DEBUG_KERNEL)
    (int*, int*, int*, int*, double*, double*, double*, double*, double*, double*, 
     double*, double*, double*, double*, double*, double*, double*, double*, double*,
     double*, double*, double*, double*, double*, double*, double*);

  void F90_FUNC(field_summary_kernel, FIELD_SUMMARY_KERNEL)
    (int*,int*,int*,int*,double*,double*,double*,double*,double*,double*,int*,
     double*,double*,double*,double*,double*,int*,int*);

  void F90_FUNC(initialise_chunk_kernel, INITIALISE_CHUNK_KERNEL)
    (int*,int*,int*,int*,double*,double*,double*,double*,double*,double*,
     double*,double*,double*,double*,double*,double*,double*);

  void F90_FUNC(generate_chunk_kernel, GENERATE_CHUNK_KERNEL)
    (int*,int*,int*,int*,double*,double*,double*,double*,double*,double*,
     double*,double*,int*,double*,double*,double*,double*,double*,double*,
     double*,double*,double*,int*,const int*,const int*,const int*);
}

#define POLY2(i, j, imin, jmin, nx) ((i - imin) + (j-jmin) * nx)

const int Cleverleaf::g_rectangle;
const int Cleverleaf::g_circle;
const int Cleverleaf::g_point;

const double Cleverleaf::g_small = 1.0e-16;
const double Cleverleaf::g_big = 1.0e+21;

tbox::StartupShutdownManager::Handler
Cleverleaf::s_initialize_handler(
    Cleverleaf::initializeCallback,
    0,
    0,
    Cleverleaf::finalizeCallback,
    tbox::StartupShutdownManager::priorityTimers);

std::shared_ptr<tbox::Timer> Cleverleaf::t_fill_boundary;

Cleverleaf::Cleverleaf(
    std::shared_ptr<tbox::Database> input_database,
    std::shared_ptr<hier::PatchHierarchy> hierarchy,
    const tbox::Dimension& dim,
    std::shared_ptr<geom::CartesianGridGeometry> grid_geometry):
  LagrangianEulerianPatchStrategy(dim),
  d_dim(dim),
  d_nghosts(dim),
  input_db(input_database),
  state_prefix("state"),
  d_velocity(new clever::pdat::CleverNodeVariable<double>(
        d_dim, "velocity", d_dim.getValue())),
  d_massflux(new clever::pdat::CleverSideVariable<double>(d_dim, "massflux", 2)),
  d_volflux(new clever::pdat::CleverSideVariable<double>(d_dim, "volflux", 2)),
  d_pressure(new clever::pdat::CleverCellVariable<double>(d_dim, "pressure", 1)),
  d_viscosity(new clever::pdat::CleverCellVariable<double>(d_dim, "viscosity", 1)),
  d_soundspeed(new clever::pdat::CleverCellVariable<double>(d_dim, "soundspeed", 1)),
  d_density(new clever::pdat::CleverCellVariable<double>(d_dim, "density", 1)),
  d_energy(new clever::pdat::CleverCellVariable<double>(d_dim, "energy", 1)),
  d_volume(new clever::pdat::CleverCellVariable<double>(d_dim, "volume", 1)),
  d_celldeltas(new clever::pdat::CleverCellVariable<double>(d_dim, "celldelta", dim.getValue())),
  d_cellcoords(new clever::pdat::CleverCellVariable<double>(d_dim, "cellcoords", dim.getValue())),
  d_vertexdeltas(new clever::pdat::CleverNodeVariable<double>(d_dim, "vertexdeltas",dim.getValue())),
  d_vertexcoords(new clever::pdat::CleverNodeVariable<double>(d_dim, "vertexcoords", dim.getValue())),
  d_workarray1(new clever::pdat::CleverNodeVariable<double>(d_dim, "workarray 1", 1)),
  d_workarray2(new clever::pdat::CleverNodeVariable<double>(d_dim, "workarray 2", 1)),
  d_workarray3(new clever::pdat::CleverNodeVariable<double>(d_dim, "workarray 3", 1)),
  d_workarray4(new clever::pdat::CleverNodeVariable<double>(d_dim, "workarray 4", 1)),
  d_workarray5(new clever::pdat::CleverNodeVariable<double>(d_dim, "workarray 5", 1)),
  d_workarray6(new clever::pdat::CleverNodeVariable<double>(d_dim, "workarray 6", 1)),
  d_workarray7(new clever::pdat::CleverNodeVariable<double>(d_dim, "workarray 7", 1)),
  d_level_indicator(new clever::pdat::CleverCellVariable<int>(d_dim, "level_indicator", 1)),
  d_exchange_fields(new int[15])
{
  d_hierarchy = hierarchy;
  d_grid_geometry = grid_geometry;

  d_nghosts = hier::IntVector(d_dim, 2);

  /*
   * Add our coarsen operators to the registry.
   */
  std::shared_ptr<hier::CoarsenOperator> vol_weighted_avg(
      new clever::geom::CartesianCleverCellDoubleVolumeWeightedAverage(dim));
  std::shared_ptr<hier::CoarsenOperator> mass_weighted_avg(
      new clever::geom::CartesianCleverCellDoubleMassWeightedAverage(dim));
  std::shared_ptr<hier::CoarsenOperator> constant_cell_coarsen(
      new clever::geom::CartesianCleverCellIntConstantCoarsen(dim));

  std::shared_ptr<hier::CoarsenOperator> ndi(
      new clever::pdat::CleverNodeInjection<double>());
  std::shared_ptr<hier::RefineOperator> cndlr(
      new clever::geom::CartesianCleverNodeDoubleLinearRefine());
  std::shared_ptr<hier::RefineOperator> cedclr(
      new clever::geom::CartesianCleverSideDoubleFirstOrderRefine());
  std::shared_ptr<hier::RefineOperator> ccdclr(
      new clever::geom::CartesianCleverCellDoubleConservativeLinearRefine());

  d_grid_geometry->addCoarsenOperator(
      typeid(clever::pdat::CleverCellVariable<double>).name(), vol_weighted_avg);
  d_grid_geometry->addCoarsenOperator(
      typeid(clever::pdat::CleverCellVariable<double>).name(), mass_weighted_avg);
  d_grid_geometry->addCoarsenOperator(
      typeid(clever::pdat::CleverCellVariable<int>).name(), constant_cell_coarsen);
  d_grid_geometry->addCoarsenOperator(
      typeid(clever::pdat::CleverNodeVariable<double>).name(), ndi);
  d_grid_geometry->addRefineOperator(
      typeid(clever::pdat::CleverNodeVariable<double>).name(), cndlr);
  d_grid_geometry->addRefineOperator(
      typeid(clever::pdat::CleverSideVariable<double>).name(), cedclr);
  d_grid_geometry->addRefineOperator(
      typeid(clever::pdat::CleverCellVariable<double>).name(), ccdclr);

  d_tag_all = input_database->getBoolWithDefault("tag_all", false);
  d_tag_q = input_database->getBoolWithDefault("tag_q", true);
  d_tag_density = input_database->getBoolWithDefault("tag_density", true);
  d_tag_energy = input_database->getBoolWithDefault("tag_energy", true);
  d_tag_pressure = input_database->getBoolWithDefault("tag_pressure", true);

  d_tag_q_threshold = input_database->getDoubleWithDefault(
      "tag_q_threshold", 0.001);
  d_tag_density_gradient = input_database->getDoubleWithDefault(
      "tag_density_threshold", 0.1);
  d_tag_energy_gradient = input_database->getDoubleWithDefault(
      "tag_energy_threshold", 0.1);
  d_tag_pressure_gradient = input_database->getDoubleWithDefault(
      "tag_pressure_threshold", 0.1);

  d_pdv_weight = input_database->getIntegerWithDefault(
      "physics_weight", 1);

  d_gravity = input_database->getBoolWithDefault("gravity", false);
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
      LagrangianEulerianLevelIntegrator::FIELD |
      LagrangianEulerianLevelIntegrator::REVERT,
      LagrangianEulerianLevelIntegrator::PRIME_CELLS_EXCH |
      LagrangianEulerianLevelIntegrator::PRE_LAGRANGE_EXCH |
      LagrangianEulerianLevelIntegrator::PRE_SWEEP_1_CELL_EXCH |
      LagrangianEulerianLevelIntegrator::PRE_SWEEP_1_MOM_EXCH |
      LagrangianEulerianLevelIntegrator::PRE_SWEEP_2_MOM_EXCH,
      d_nghosts, 
      d_grid_geometry);

  integrator->registerVariable(
      d_energy,
      LagrangianEulerianLevelIntegrator::FIELD |
      LagrangianEulerianLevelIntegrator::REVERT,
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
      d_workarray1,
      LagrangianEulerianLevelIntegrator::NORMAL,
      LagrangianEulerianLevelIntegrator::NO_EXCH,
      d_nghosts,
      d_grid_geometry);

  integrator->registerVariable(
      d_workarray2,
      LagrangianEulerianLevelIntegrator::NORMAL,
      LagrangianEulerianLevelIntegrator::NO_EXCH,
      d_nghosts,
      d_grid_geometry);

  integrator->registerVariable(
      d_workarray3,
      LagrangianEulerianLevelIntegrator::NORMAL,
      LagrangianEulerianLevelIntegrator::NO_EXCH,
      d_nghosts,
      d_grid_geometry);

  integrator->registerVariable(
      d_workarray4,
      LagrangianEulerianLevelIntegrator::NORMAL,
      LagrangianEulerianLevelIntegrator::NO_EXCH,
      d_nghosts,
      d_grid_geometry);

  integrator->registerVariable(
      d_workarray5,
      LagrangianEulerianLevelIntegrator::NORMAL,
      LagrangianEulerianLevelIntegrator::NO_EXCH,
      d_nghosts,
      d_grid_geometry);

  integrator->registerVariable(
      d_workarray6,
      LagrangianEulerianLevelIntegrator::NORMAL,
      LagrangianEulerianLevelIntegrator::NO_EXCH,
      d_nghosts,
      d_grid_geometry);

  integrator->registerVariable(
      d_workarray7,
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
    d_visit_writer->registerDerivedPlotQuantity(
        "Pressure",
        "SCALAR",
        d_pressure.get());

    d_visit_writer->registerDerivedPlotQuantity(
        "Viscosity",
        "SCALAR",
        d_viscosity.get());

    d_visit_writer->registerDerivedPlotQuantity(
        "Soundspeed",
        "SCALAR",
        d_soundspeed.get());

    d_visit_writer->registerDerivedPlotQuantity(
        "Density",
        "SCALAR",
        d_density.get());

    d_visit_writer->registerDerivedPlotQuantity(
        "Energy",
        "SCALAR",
        d_energy.get());

    d_visit_writer->registerDerivedPlotQuantity(
        "Volume",
        "SCALAR",
        d_volume.get());

    d_visit_writer->registerDerivedPlotQuantity(
        "Level Indicator",
        "SCALAR",
        d_level_indicator.get());

    d_visit_writer->registerDerivedPlotQuantity(
        "Velocity",
        "VECTOR",
        d_velocity.get(),
        1.0,
        "NODE");

    d_visit_writer->registerDerivedPlotQuantity(
        "Vertexcoords",
        "VECTOR",
        d_vertexcoords.get(),
        1.0,
        "NODE");

    d_visit_writer->registerDerivedPlotQuantity(
        "Cellcoords",
        "VECTOR",
        d_cellcoords.get());
  }
}

void Cleverleaf::registerVisItDataWriter(
    std::shared_ptr<appu::VisItDataWriter> writer)
{
  d_visit_writer = writer;
}

void Cleverleaf::initializeDataOnPatch(
    hier::Patch& patch,
    double init_data_time,
    bool initial_time)
{
  std::shared_ptr<clever::pdat::CleverCellData<double> > cell_deltas(
      SHARED_PTR_CAST(clever::pdat::CleverCellData<double>,
        patch.getPatchData(d_celldeltas,getCurrentDataContext())));

  std::shared_ptr<clever::pdat::CleverCellData<double> > cell_coordinates(
      SHARED_PTR_CAST(clever::pdat::CleverCellData<double> ,
        patch.getPatchData(d_cellcoords, getCurrentDataContext())));

  std::shared_ptr<clever::pdat::CleverNodeData<double> > vertex_deltas(
      SHARED_PTR_CAST(clever::pdat::CleverNodeData<double> ,
        patch.getPatchData(d_vertexdeltas, getCurrentDataContext())));

  std::shared_ptr<clever::pdat::CleverNodeData<double> > vertex_coordinates(
      SHARED_PTR_CAST(clever::pdat::CleverNodeData<double> ,
        patch.getPatchData( d_vertexcoords, getCurrentDataContext())));

  std::shared_ptr<clever::pdat::CleverCellData<double> > volume(
      SHARED_PTR_CAST(clever::pdat::CleverCellData<double> ,
        patch.getPatchData(d_volume, getCurrentDataContext())));

  const hier::Index ifirst = patch.getBox().lower();
  const hier::Index ilast = patch.getBox().upper();

  int x_min = ifirst(0);
  int x_max = ilast(0);
  int y_min = ifirst(1);
  int y_max = ilast(1);

  /*
   * Get the patch geometry - this stores the coordinates of the patch corner
   * in our Cartesian geometry, and the as well as the cell sizes.
   */
  const std::shared_ptr<geom::CartesianPatchGeometry> pgeom(
      SHARED_PTR_CAST(geom::CartesianPatchGeometry,
        patch.getPatchGeometry()));

  /* Array containing dx and dy */
  const double* dxs = pgeom->getDx();

  /* Lower left coordinate of the patch */
  const double* coords = pgeom->getXLower();

  double physical_xmin = coords[0];
  double physical_ymin = coords[1];

  double dx = dxs[0];
  double dy = dxs[1];

  F90_FUNC(initialise_chunk_kernel,INITIALISE_CHUNK_KERNEL)
    (&x_min,
     &x_max,
     &y_min,
     &y_max,
     &physical_xmin,
     &physical_ymin,
     &dx,
     &dy,
     vertex_coordinates->getPointer(0),
     vertex_deltas->getPointer(0),
     vertex_coordinates->getPointer(1),
     vertex_deltas->getPointer(1),
     cell_coordinates->getPointer(0),
     cell_deltas->getPointer(0),
     cell_coordinates->getPointer(1),
     cell_deltas->getPointer(1),
     volume->getPointer());

  if (initial_time) {
    std::shared_ptr<clever::pdat::CleverNodeData<double> > velocity(
        SHARED_PTR_CAST(clever::pdat::CleverNodeData<double> ,
          patch.getPatchData(d_velocity, getCurrentDataContext())));

    std::shared_ptr<clever::pdat::CleverSideData<double> > mass_flux(
        SHARED_PTR_CAST(clever::pdat::CleverSideData<double> ,
          patch.getPatchData(d_massflux, getCurrentDataContext())));

    std::shared_ptr<clever::pdat::CleverSideData<double> > volume_flux(
        SHARED_PTR_CAST(clever::pdat::CleverSideData<double> ,
          patch.getPatchData(d_volflux, getCurrentDataContext())));

    std::shared_ptr<clever::pdat::CleverCellData<double> > pressure(
        SHARED_PTR_CAST(clever::pdat::CleverCellData<double> ,
          patch.getPatchData(d_pressure, getCurrentDataContext())));

    std::shared_ptr<clever::pdat::CleverCellData<double> > viscosity(
        SHARED_PTR_CAST(clever::pdat::CleverCellData<double> ,
          patch.getPatchData(d_viscosity, getCurrentDataContext())));

    std::shared_ptr<clever::pdat::CleverCellData<double> > soundspeed(
        SHARED_PTR_CAST(clever::pdat::CleverCellData<double> ,
          patch.getPatchData(d_soundspeed, getCurrentDataContext())));

    std::shared_ptr<clever::pdat::CleverCellData<double> > density(
        SHARED_PTR_CAST(clever::pdat::CleverCellData<double> ,
          patch.getPatchData(d_density, getCurrentDataContext())));

    std::shared_ptr<clever::pdat::CleverCellData<double> > energy(
        SHARED_PTR_CAST(clever::pdat::CleverCellData<double> ,
          patch.getPatchData(d_energy, getCurrentDataContext())));

    velocity->fillAll(0.0);
    mass_flux->fillAll(0.0);
    volume_flux->fillAll(0.0);
    viscosity->fillAll(0.0);
    soundspeed->fillAll(0.0);
    pressure->fillAll(0.0);
    energy->fillAll(0.0);
    density->fillAll(0.0);

    std::shared_ptr<tbox::Database> states_db = input_db->getDatabase(
        "states");

    int number_of_states = states_db->getInteger("num_states");

    double* state_density = new double[number_of_states];
    double* state_energy = new double[number_of_states];
    double* state_xvel = new double[number_of_states];
    double* state_yvel = new double[number_of_states];
    double* state_xmin = new double[number_of_states];
    double* state_ymin = new double[number_of_states];
    double* state_xmax = new double[number_of_states];
    double* state_ymax = new double[number_of_states];
    double* state_radius = new double[number_of_states];
    int* state_geometry = new int[number_of_states];

    for(int state = 0; state < number_of_states; state++) {
      std::ostringstream state_stream;
      state_stream << state_prefix << state;

      std::shared_ptr<tbox::Database> current_state = states_db->getDatabase(
          state_stream.str());

      std::string state_geometry_string;
      state_geometry_string = current_state->getStringWithDefault(
          "geometry", "RECTANGLE");


      if (state_geometry_string.compare("RECTANGLE") == 0) {
        state_geometry[state] = g_rectangle;
      } else if (state_geometry_string.compare("CIRCLE") == 0) {
        state_geometry[state] = g_circle;
      } else if (state_geometry_string.compare("POINT") == 0) {
        state_geometry[state] = g_point;
      }

      if(state_geometry[state] == g_circle ||
          state_geometry[state] == g_point) {
        double* center = new double[2];
        current_state->getDoubleArray("center", center, 2);

        state_xmin[state] = center[0];
        state_ymin[state] = center[1];

        state_xmax[state] = -1;
        state_ymax[state] = -1;
      } else {
        if(state == 0) {
          state_xmin[state] = -1;
          state_ymin[state] = -1;

          state_xmax[state] = -1;
          state_ymax[state] = -1;
        } else {
          double* state_min = new double[2];
          double* state_max = new double[2];
          current_state->getDoubleArray("min", state_min, 2);
          current_state->getDoubleArray("max", state_max, 2);

          state_xmin[state] = state_min[0];
          state_ymin[state] = state_min[1];

          state_xmax[state] = state_max[0];
          state_ymax[state] = state_max[1];
        }
      }

      state_density[state] = current_state->getDouble("density");
      state_energy[state] = current_state->getDouble("energy");
      state_xvel[state] = current_state->getDoubleWithDefault("xvel", 0.0);
      state_yvel[state] = current_state->getDoubleWithDefault("yvel", 0.0);
      state_radius[state] = current_state->getDoubleWithDefault("radius", -1);
    }

    F90_FUNC(generate_chunk_kernel,GENERATE_CHUNK_KERNEL)
      (&x_min,
       &x_max,
       &y_min,
       &y_max,
       vertex_coordinates->getPointer(0),
       vertex_coordinates->getPointer(1),
       cell_coordinates->getPointer(0),
       cell_coordinates->getPointer(1),
       density->getPointer(),
       energy->getPointer(),
       velocity->getPointer(0),
       velocity->getPointer(1),
       &number_of_states,
       state_density,
       state_energy,
       state_xvel,
       state_yvel,
       state_xmin,
       state_xmax,
       state_ymin,
       state_ymax,
       state_radius,
       state_geometry,
       &g_rectangle,
       &g_circle,
       &g_point);

    delete[] state_density;
    delete[] state_energy;
    delete[] state_xvel;
    delete[] state_yvel;
    delete[] state_xmin;
    delete[] state_ymin;
    delete[] state_xmax;
    delete[] state_ymax;
    delete[] state_radius;
    delete[] state_geometry;
  }

  ideal_gas_knl(patch, false);
}

void Cleverleaf::accelerate(hier::Patch& patch, double dt)
{
  std::shared_ptr<clever::pdat::CleverCellData<double> > density0(
      SHARED_PTR_CAST(clever::pdat::CleverCellData<double> ,
        patch.getPatchData(d_density, getCurrentDataContext())));

  std::shared_ptr<clever::pdat::CleverCellData<double> > volume(
      SHARED_PTR_CAST(clever::pdat::CleverCellData<double> ,
        patch.getPatchData(d_volume, getCurrentDataContext())));

  std::shared_ptr<clever::pdat::CleverCellData<double> > cell_deltas(
      SHARED_PTR_CAST(clever::pdat::CleverCellData<double> ,
        patch.getPatchData(d_celldeltas, getCurrentDataContext())));

  std::shared_ptr<clever::pdat::CleverCellData<double> > pressure(
      SHARED_PTR_CAST(clever::pdat::CleverCellData<double> ,
        patch.getPatchData(d_pressure, getCurrentDataContext())));

  std::shared_ptr<clever::pdat::CleverCellData<double> > viscosity(
      SHARED_PTR_CAST(clever::pdat::CleverCellData<double> ,
        patch.getPatchData(d_viscosity, getCurrentDataContext())));

  std::shared_ptr<clever::pdat::CleverNodeData<double> > velocity0(
      SHARED_PTR_CAST(clever::pdat::CleverNodeData<double> ,
        patch.getPatchData(d_velocity, getCurrentDataContext())));

  std::shared_ptr<clever::pdat::CleverNodeData<double> > velocity1(
      SHARED_PTR_CAST(clever::pdat::CleverNodeData<double> ,
        patch.getPatchData(d_velocity, getNewDataContext())));

  std::shared_ptr<clever::pdat::CleverNodeData<double> > step_by_mass(
      SHARED_PTR_CAST(clever::pdat::CleverNodeData<double> ,
        patch.getPatchData(d_workarray1, getCurrentDataContext())));

  const hier::Index ifirst = patch.getBox().lower();
  const hier::Index ilast = patch.getBox().upper();

  int xmin = ifirst(0);
  int xmax = ilast(0);
  int ymin = ifirst(1);
  int ymax = ilast(1);

  int igravity;

  if (d_gravity) {
    igravity = 1;
  } else {
    igravity = 0;
  }

  F90_FUNC(accelerate_kernel,ACCELERATE_KERNEL)
    (&xmin,
     &xmax,
     &ymin,
     &ymax,
     &dt,
     &igravity,
     cell_deltas->getPointer(1),
     cell_deltas->getPointer(0),
     volume->getPointer(),
     density0->getPointer(),
     pressure->getPointer(),
     viscosity->getPointer(),
     velocity0->getPointer(0),
     velocity0->getPointer(1),
     velocity1->getPointer(0),
     velocity1->getPointer(1),
     step_by_mass->getPointer());
}

void Cleverleaf::ideal_gas_knl(hier::Patch& patch, bool predict)
{
  std::shared_ptr<clever::pdat::CleverCellData<double> > density;
  std::shared_ptr<clever::pdat::CleverCellData<double> > energy;

  if (predict) {
    std::shared_ptr<clever::pdat::CleverCellData<double> > density1(
        SHARED_PTR_CAST(clever::pdat::CleverCellData<double> ,
          patch.getPatchData(d_density, getNewDataContext())));
    density = density1;

    std::shared_ptr<clever::pdat::CleverCellData<double> > energy1(
        SHARED_PTR_CAST(clever::pdat::CleverCellData<double> ,
          patch.getPatchData(d_energy, getNewDataContext())));
    energy = energy1;
  } else {
    std::shared_ptr<clever::pdat::CleverCellData<double> > density0(
        SHARED_PTR_CAST(clever::pdat::CleverCellData<double> ,
          patch.getPatchData(d_density, getCurrentDataContext())));
    density = density0;

    std::shared_ptr<clever::pdat::CleverCellData<double> > energy0(
        SHARED_PTR_CAST(clever::pdat::CleverCellData<double> ,
          patch.getPatchData(d_energy, getCurrentDataContext())));
    energy = energy0;
  }

  std::shared_ptr<clever::pdat::CleverCellData<double> > pressure(
      SHARED_PTR_CAST(clever::pdat::CleverCellData<double> ,
        patch.getPatchData(d_pressure, getCurrentDataContext())));

  std::shared_ptr<clever::pdat::CleverCellData<double> > soundspeed(
      SHARED_PTR_CAST(clever::pdat::CleverCellData<double> ,
        patch.getPatchData(d_soundspeed, getCurrentDataContext())));

  const hier::Index ifirst = patch.getBox().lower();
  const hier::Index ilast = patch.getBox().upper();

  int xmin = ifirst(0);
  int xmax = ilast(0);
  int ymin = ifirst(1);
  int ymax = ilast(1);

  F90_FUNC(ideal_gas_kernel,IDEAL_GAS_KERNEL)
    (&xmin,
     &xmax,
     &ymin,
     &ymax,
     density->getPointer(),
     energy->getPointer(),
     pressure->getPointer(),
     soundspeed->getPointer());
}

void Cleverleaf::viscosity_knl(hier::Patch& patch)
{
  std::shared_ptr<clever::pdat::CleverCellData<double> > density0(
      SHARED_PTR_CAST(clever::pdat::CleverCellData<double> ,
        patch.getPatchData(d_density, getCurrentDataContext())));

  std::shared_ptr<clever::pdat::CleverCellData<double> > pressure(
      SHARED_PTR_CAST(clever::pdat::CleverCellData<double> ,
        patch.getPatchData(d_pressure, getCurrentDataContext())));

  std::shared_ptr<clever::pdat::CleverCellData<double> > viscosity(
      SHARED_PTR_CAST(clever::pdat::CleverCellData<double> ,
        patch.getPatchData(d_viscosity, getCurrentDataContext())));

  std::shared_ptr<clever::pdat::CleverNodeData<double> > velocity0(
      SHARED_PTR_CAST(clever::pdat::CleverNodeData<double> ,
        patch.getPatchData(d_velocity, getCurrentDataContext())));

  std::shared_ptr<clever::pdat::CleverCellData<double> > cell_deltas(
      SHARED_PTR_CAST(clever::pdat::CleverCellData<double> ,
        patch.getPatchData(d_celldeltas, getCurrentDataContext())));

  const hier::Index ifirst = patch.getBox().lower();
  const hier::Index ilast = patch.getBox().upper();

  int xmin = ifirst(0); 
  int xmax = ilast(0); 
  int ymin = ifirst(1); 
  int ymax = ilast(1); 

  F90_FUNC(viscosity_kernel,VISCOSITY_KERNEL)
    (&xmin,
     &xmax,
     &ymin,
     &ymax,
     cell_deltas->getPointer(0),
     cell_deltas->getPointer(1),
     density0->getPointer(),
     pressure->getPointer(),
     viscosity->getPointer(),
     velocity0->getPointer(0),
     velocity0->getPointer(1)); 
}

double Cleverleaf::calc_dt_knl(hier::Patch& patch)
{
  std::shared_ptr<clever::pdat::CleverCellData<double> > density0(
      SHARED_PTR_CAST(clever::pdat::CleverCellData<double> ,
        patch.getPatchData(d_density, getCurrentDataContext())));

  std::shared_ptr<clever::pdat::CleverCellData<double> > energy0(
      SHARED_PTR_CAST(clever::pdat::CleverCellData<double> ,
        patch.getPatchData(d_energy, getCurrentDataContext())));

  std::shared_ptr<clever::pdat::CleverCellData<double> > pressure(
      SHARED_PTR_CAST(clever::pdat::CleverCellData<double> ,
        patch.getPatchData(d_pressure, getCurrentDataContext())));

  std::shared_ptr<clever::pdat::CleverCellData<double> > soundspeed(
      SHARED_PTR_CAST(clever::pdat::CleverCellData<double> ,
        patch.getPatchData(d_soundspeed, getCurrentDataContext())));

  std::shared_ptr<clever::pdat::CleverCellData<double> > viscosity(
      SHARED_PTR_CAST(clever::pdat::CleverCellData<double> ,
        patch.getPatchData(d_viscosity, getCurrentDataContext())));

  std::shared_ptr<clever::pdat::CleverNodeData<double> > velocity0(
      SHARED_PTR_CAST(clever::pdat::CleverNodeData<double> ,
        patch.getPatchData(d_velocity, getCurrentDataContext())));

  std::shared_ptr<clever::pdat::CleverCellData<double> > cell_deltas(
      SHARED_PTR_CAST(clever::pdat::CleverCellData<double> ,
        patch.getPatchData(d_celldeltas, getCurrentDataContext())));

  std::shared_ptr<clever::pdat::CleverCellData<double> > volume(
      SHARED_PTR_CAST(clever::pdat::CleverCellData<double> ,
        patch.getPatchData(d_volume, getCurrentDataContext())));

  std::shared_ptr<clever::pdat::CleverCellData<double> > cell_coordinates(
      SHARED_PTR_CAST(clever::pdat::CleverCellData<double> ,
        patch.getPatchData(d_cellcoords, getCurrentDataContext())));

  std::shared_ptr<clever::pdat::CleverNodeData<double> > dtmin(
      SHARED_PTR_CAST(clever::pdat::CleverNodeData<double> ,
        patch.getPatchData(d_workarray1, getCurrentDataContext())));

  const hier::Index ifirst = patch.getBox().lower();
  const hier::Index ilast = patch.getBox().upper();

  int xmin = ifirst(0);
  int xmax = ilast(0);
  int ymin = ifirst(1);
  int ymax = ilast(1);

  double dt_min_val = 1.0e+21;
  double dtc_safe = 0.7;
  double dtu_safe = 0.5;
  double dtv_safe = 0.5;
  double dtdiv_safe = 0.7;
  double dt_min = 0.0000001;

  int dtl_control;
  int jldt, kldt;
  double xl_pos, yl_pos;
  int small;

  F90_FUNC(calc_dt_kernel, CALC_DT_KERNEL)
    (&xmin,
     &xmax,
     &ymin,
     &ymax,
     g_small,
     g_big,
     &dt_min,
     &dtc_safe,
     &dtu_safe,
     &dtv_safe,
     &dtdiv_safe,
     cell_deltas->getPointer(1),
     cell_deltas->getPointer(0),
     cell_coordinates->getPointer(0),
     cell_coordinates->getPointer(1),
     cell_deltas->getPointer(0),
     cell_deltas->getPointer(1),
     volume->getPointer(),
     density0->getPointer(),
     energy0->getPointer(),
     pressure->getPointer(),
     viscosity->getPointer(),
     soundspeed->getPointer(),
     velocity0->getPointer(0),
     velocity0->getPointer(1),
     dtmin->getPointer(),
     &dt_min_val,
     &dtl_control,
     &xl_pos,
     &yl_pos,
     &jldt,
     &kldt,
     &small);

  return dt_min_val;
}

void Cleverleaf::pdv_knl(hier::Patch& patch, double dt, bool predict)
{
  std::shared_ptr<clever::pdat::CleverCellData<double> > cell_deltas(
      SHARED_PTR_CAST(clever::pdat::CleverCellData<double> ,
        patch.getPatchData(d_celldeltas, getCurrentDataContext())));

  std::shared_ptr<clever::pdat::CleverCellData<double> > volume(
      SHARED_PTR_CAST(clever::pdat::CleverCellData<double> ,
        patch.getPatchData(d_volume, getCurrentDataContext())));

  std::shared_ptr<clever::pdat::CleverCellData<double> > density0(
      SHARED_PTR_CAST(clever::pdat::CleverCellData<double> ,
        patch.getPatchData(d_density, getCurrentDataContext())));

  std::shared_ptr<clever::pdat::CleverCellData<double> > density1(
      SHARED_PTR_CAST(clever::pdat::CleverCellData<double> ,
        patch.getPatchData(d_density, getNewDataContext())));

  std::shared_ptr<clever::pdat::CleverCellData<double> > energy0(
      SHARED_PTR_CAST(clever::pdat::CleverCellData<double> ,
        patch.getPatchData(d_energy, getCurrentDataContext())));

  std::shared_ptr<clever::pdat::CleverCellData<double> > energy1(
      SHARED_PTR_CAST(clever::pdat::CleverCellData<double> ,
        patch.getPatchData(d_energy, getNewDataContext())));

  std::shared_ptr<clever::pdat::CleverCellData<double> > pressure(
      SHARED_PTR_CAST(clever::pdat::CleverCellData<double> ,
        patch.getPatchData(d_pressure, getCurrentDataContext())));

  std::shared_ptr<clever::pdat::CleverCellData<double> > viscosity(
      SHARED_PTR_CAST(clever::pdat::CleverCellData<double> ,
        patch.getPatchData(d_viscosity, getCurrentDataContext())));

  std::shared_ptr<clever::pdat::CleverNodeData<double> > velocity0(
      SHARED_PTR_CAST(clever::pdat::CleverNodeData<double> ,
        patch.getPatchData(d_velocity, getCurrentDataContext())));

  std::shared_ptr<clever::pdat::CleverNodeData<double> > velocity1(
      SHARED_PTR_CAST(clever::pdat::CleverNodeData<double> ,
        patch.getPatchData(d_velocity, getNewDataContext())));

  std::shared_ptr<clever::pdat::CleverNodeData<double> > volume_change(
      SHARED_PTR_CAST(clever::pdat::CleverNodeData<double> ,
        patch.getPatchData(d_workarray1, getCurrentDataContext())));

  const hier::Index ifirst = patch.getBox().lower();
  const hier::Index ilast = patch.getBox().upper();

  int xmin = ifirst(0);
  int xmax = ilast(0);
  int ymin = ifirst(1);
  int ymax = ilast(1);

  int prdct;

  if (predict)
    prdct = 1;
  else
    prdct = 0;

  /*
   * Iterate over PdV kernel to increase cost of physics.
   */
  for (int i = 0; i < d_pdv_weight; i++) {
    F90_FUNC(pdv_kernel, PDV_KERNEL)
      (&prdct,
       &xmin,
       &xmax,
       &ymin,
       &ymax,
       &dt,
       cell_deltas->getPointer(1),
       cell_deltas->getPointer(0),
       volume->getPointer(),
       density0->getPointer(),
       density1->getPointer(),
       energy0->getPointer(),
       energy1->getPointer(),
       pressure->getPointer(),
       viscosity->getPointer(),
       velocity0->getPointer(0),
       velocity1->getPointer(0),
       velocity0->getPointer(1),
       velocity1->getPointer(1),
       volume_change->getPointer());
  }
}

void Cleverleaf::flux_calc_knl(hier::Patch& patch, double dt)
{

  std::shared_ptr<clever::pdat::CleverNodeData<double> > velocity0(
      SHARED_PTR_CAST(clever::pdat::CleverNodeData<double> ,
        patch.getPatchData(d_velocity, getCurrentDataContext())));

  std::shared_ptr<clever::pdat::CleverNodeData<double> > velocity1(
      SHARED_PTR_CAST(clever::pdat::CleverNodeData<double> ,
        patch.getPatchData(d_velocity, getNewDataContext())));

  std::shared_ptr<clever::pdat::CleverSideData<double> > volume_flux(
      SHARED_PTR_CAST(clever::pdat::CleverSideData<double> ,
        patch.getPatchData(d_volflux, getCurrentDataContext())));

  std::shared_ptr<clever::pdat::CleverCellData<double> > cell_deltas(
      SHARED_PTR_CAST(clever::pdat::CleverCellData<double> ,
        patch.getPatchData(d_celldeltas, getCurrentDataContext())));

  const hier::Index ifirst = patch.getBox().lower();
  const hier::Index ilast = patch.getBox().upper();

  int xmin = ifirst(0); 
  int xmax = ilast(0); 
  int ymin = ifirst(1); 
  int ymax = ilast(1); 

  F90_FUNC(flux_calc_kernel, FLUX_CALC_KERNEL)
    (&xmin,
     &xmax,
     &ymin,
     &ymax,
     &dt,
     cell_deltas->getPointer(1),
     cell_deltas->getPointer(0),
     velocity0->getPointer(0),
     velocity0->getPointer(1),
     velocity1->getPointer(0),
     velocity1->getPointer(1),
     volume_flux->getPointer(0),
     volume_flux->getPointer(1));
}

void Cleverleaf::advec_cell(hier::Patch& patch, int sweep_number, ADVEC_DIR dir)
{
  std::shared_ptr<clever::pdat::CleverCellData<double> > density1(
      SHARED_PTR_CAST(clever::pdat::CleverCellData<double> ,
        patch.getPatchData(d_density, getNewDataContext())));

  std::shared_ptr<clever::pdat::CleverCellData<double> > energy1(
      SHARED_PTR_CAST(clever::pdat::CleverCellData<double> ,
        patch.getPatchData(d_energy, getNewDataContext())));

  std::shared_ptr<clever::pdat::CleverCellData<double> > volume(
      SHARED_PTR_CAST(clever::pdat::CleverCellData<double> ,
        patch.getPatchData(d_volume, getCurrentDataContext())));

  std::shared_ptr<clever::pdat::CleverSideData<double> > volume_flux(
      SHARED_PTR_CAST(clever::pdat::CleverSideData<double> ,
        patch.getPatchData(d_volflux, getCurrentDataContext())));

  std::shared_ptr<clever::pdat::CleverSideData<double> > mass_flux(
      SHARED_PTR_CAST(clever::pdat::CleverSideData<double> ,
        patch.getPatchData(d_massflux, getCurrentDataContext())));

  std::shared_ptr<clever::pdat::CleverNodeData<double> > vertex_deltas(
      SHARED_PTR_CAST(clever::pdat::CleverNodeData<double> ,
        patch.getPatchData(d_vertexdeltas, getCurrentDataContext())));

  std::shared_ptr<clever::pdat::CleverNodeData<double> > pre_volume(
      SHARED_PTR_CAST(clever::pdat::CleverNodeData<double> ,
        patch.getPatchData(d_workarray1, getCurrentDataContext())));

  std::shared_ptr<clever::pdat::CleverNodeData<double> > post_volume(
      SHARED_PTR_CAST(clever::pdat::CleverNodeData<double> ,
        patch.getPatchData(d_workarray2, getCurrentDataContext())));

  std::shared_ptr<clever::pdat::CleverNodeData<double> > pre_mass(
      SHARED_PTR_CAST(clever::pdat::CleverNodeData<double> ,
        patch.getPatchData(d_workarray3, getCurrentDataContext())));

  std::shared_ptr<clever::pdat::CleverNodeData<double> > post_mass(
      SHARED_PTR_CAST(clever::pdat::CleverNodeData<double> ,
        patch.getPatchData(d_workarray4, getCurrentDataContext())));

  std::shared_ptr<clever::pdat::CleverNodeData<double> > advected_volume(
      SHARED_PTR_CAST(clever::pdat::CleverNodeData<double> ,
        patch.getPatchData(d_workarray5, getCurrentDataContext())));

  std::shared_ptr<clever::pdat::CleverNodeData<double> > post_energy(
      SHARED_PTR_CAST(clever::pdat::CleverNodeData<double> ,
        patch.getPatchData(d_workarray6, getCurrentDataContext())));

  std::shared_ptr<clever::pdat::CleverNodeData<double> > energy_flux(
      SHARED_PTR_CAST(clever::pdat::CleverNodeData<double> ,
        patch.getPatchData(d_workarray7, getCurrentDataContext())));

  const hier::Index ifirst = patch.getBox().lower();
  const hier::Index ilast = patch.getBox().upper();

  int xmin = ifirst(0); 
  int xmax = ilast(0); 
  int ymin = ifirst(1); 
  int ymax = ilast(1); 


  int idir;
  if(dir == X)
    idir = 1;
  else idir = 2;

  F90_FUNC(advec_cell_kernel, ADVEC_CELL_KERNEL)
    (&xmin,
     &xmax,
     &ymin,
     &ymax,
     &idir,
     &sweep_number,
     vertex_deltas->getPointer(0),
     vertex_deltas->getPointer(1),
     volume->getPointer(),
     density1->getPointer(),
     energy1->getPointer(),
     mass_flux->getPointer(0),
     volume_flux->getPointer(0),
     mass_flux->getPointer(1),
     volume_flux->getPointer(1),
     pre_volume->getPointer(),
     post_volume->getPointer(),
     pre_mass->getPointer(),
     post_mass->getPointer(),
     advected_volume->getPointer(),
     post_energy->getPointer(),
     energy_flux->getPointer());
}

void Cleverleaf::advec_mom(
    hier::Patch& patch,
    int sweep_number,
    ADVEC_DIR direction,
    ADVEC_DIR which_vel)
{
  std::shared_ptr<clever::pdat::CleverCellData<double> > density1(
      SHARED_PTR_CAST(clever::pdat::CleverCellData<double> ,
        patch.getPatchData(d_density, getNewDataContext())));

  std::shared_ptr<clever::pdat::CleverCellData<double> > volume(
      SHARED_PTR_CAST(clever::pdat::CleverCellData<double> ,
        patch.getPatchData(d_volume, getCurrentDataContext())));

  std::shared_ptr<clever::pdat::CleverNodeData<double> > velocity1(
      SHARED_PTR_CAST(clever::pdat::CleverNodeData<double> ,
        patch.getPatchData(d_velocity, getNewDataContext())));

  std::shared_ptr<clever::pdat::CleverSideData<double> > volume_flux(
      SHARED_PTR_CAST(clever::pdat::CleverSideData<double> ,
        patch.getPatchData(d_volflux, getCurrentDataContext())));

  std::shared_ptr<clever::pdat::CleverSideData<double> > mass_flux(
      SHARED_PTR_CAST(clever::pdat::CleverSideData<double> ,
        patch.getPatchData(d_massflux, getCurrentDataContext())));

  std::shared_ptr<clever::pdat::CleverCellData<double> > cell_deltas(
      SHARED_PTR_CAST(clever::pdat::CleverCellData<double> ,
        patch.getPatchData(d_celldeltas, getCurrentDataContext())));

  std::shared_ptr<clever::pdat::CleverNodeData<double> > node_flux(
      SHARED_PTR_CAST(clever::pdat::CleverNodeData<double> ,
        patch.getPatchData(d_workarray1, getCurrentDataContext())));

  std::shared_ptr<clever::pdat::CleverNodeData<double> > node_mass_post(
      SHARED_PTR_CAST(clever::pdat::CleverNodeData<double> ,
        patch.getPatchData(d_workarray2, getCurrentDataContext())));

  std::shared_ptr<clever::pdat::CleverNodeData<double> > node_mass_pre(
      SHARED_PTR_CAST(clever::pdat::CleverNodeData<double> ,
        patch.getPatchData(d_workarray3, getCurrentDataContext())));

  std::shared_ptr<clever::pdat::CleverNodeData<double> > advection_velocity(
      SHARED_PTR_CAST(clever::pdat::CleverNodeData<double> ,
        patch.getPatchData(d_workarray4, getCurrentDataContext())));

  std::shared_ptr<clever::pdat::CleverNodeData<double> > momentum_flux(
      SHARED_PTR_CAST(clever::pdat::CleverNodeData<double> ,
        patch.getPatchData(d_workarray5, getCurrentDataContext())));

  std::shared_ptr<clever::pdat::CleverNodeData<double> > pre_volume(
      SHARED_PTR_CAST(clever::pdat::CleverNodeData<double> ,
        patch.getPatchData(d_workarray6, getCurrentDataContext())));

  std::shared_ptr<clever::pdat::CleverNodeData<double> > post_volume(
      SHARED_PTR_CAST(clever::pdat::CleverNodeData<double> ,
        patch.getPatchData(d_workarray7, getCurrentDataContext())));

  const hier::Index ifirst = patch.getBox().lower();
  const hier::Index ilast = patch.getBox().upper();

  int xmin = ifirst(0); 
  int xmax = ilast(0); 
  int ymin = ifirst(1); 
  int ymax = ilast(1); 

  int iwhich, idir;

  if (which_vel == X)
    iwhich = 1;
  else iwhich = 2;

  if (direction == X)
    idir = 1;
  else idir = 2;

  F90_FUNC(advec_mom_kernel, ADVEC_MOM_KERNEL)
    (&xmin,
     &xmax,
     &ymin,
     &ymax,
     velocity1->getPointer(0),
     velocity1->getPointer(1),
     mass_flux->getPointer(0),
     volume_flux->getPointer(0),
     mass_flux->getPointer(1),
     volume_flux->getPointer(1),
     volume->getPointer(),
     density1->getPointer(),
     node_flux->getPointer(),
     node_mass_post->getPointer(),
     node_mass_pre->getPointer(),
     advection_velocity->getPointer(),
     momentum_flux->getPointer(),
     pre_volume->getPointer(),
     post_volume->getPointer(),
     cell_deltas->getPointer(0),
     cell_deltas->getPointer(1),
     &iwhich,
     &sweep_number,
     &idir);
}

void Cleverleaf::setPhysicalBoundaryConditions(
    hier::Patch& patch,
    const double fill_time,
    const hier::IntVector& ghost_width_to_fill)
{
  t_fill_boundary->start();

  std::shared_ptr<clever::pdat::CleverCellData<double> > v_pressure(
      SHARED_PTR_CAST(clever::pdat::CleverCellData<double> ,
        patch.getPatchData(d_pressure, getCurrentDataContext())));

  std::shared_ptr<clever::pdat::CleverCellData<double> > v_density0(
      SHARED_PTR_CAST(clever::pdat::CleverCellData<double> ,
        patch.getPatchData(d_density, getCurrentDataContext())));

  std::shared_ptr<clever::pdat::CleverCellData<double> > v_density1(
      SHARED_PTR_CAST(clever::pdat::CleverCellData<double> ,
        patch.getPatchData(d_density, getNewDataContext())));

  std::shared_ptr<clever::pdat::CleverCellData<double> > v_energy0(
      SHARED_PTR_CAST(clever::pdat::CleverCellData<double> ,
        patch.getPatchData(d_energy, getCurrentDataContext())));

  std::shared_ptr<clever::pdat::CleverCellData<double> > v_energy1(
      SHARED_PTR_CAST(clever::pdat::CleverCellData<double> ,
        patch.getPatchData(d_energy, getNewDataContext())));

  std::shared_ptr<clever::pdat::CleverCellData<double> > v_viscosity(
      SHARED_PTR_CAST(clever::pdat::CleverCellData<double> ,
        patch.getPatchData(d_viscosity, getCurrentDataContext())));

  std::shared_ptr<clever::pdat::CleverNodeData<double> > v_vel0(
      SHARED_PTR_CAST(clever::pdat::CleverNodeData<double> ,
        patch.getPatchData(d_velocity, getCurrentDataContext())));

  std::shared_ptr<clever::pdat::CleverNodeData<double> > v_vel1(
      SHARED_PTR_CAST(clever::pdat::CleverNodeData<double> ,
        patch.getPatchData(d_velocity, getNewDataContext())));

  std::shared_ptr<clever::pdat::CleverSideData<double> > v_massflux( 
      SHARED_PTR_CAST(clever::pdat::CleverSideData<double> ,
        patch.getPatchData(d_massflux, getCurrentDataContext())));

  std::shared_ptr<clever::pdat::CleverSideData<double> > v_volflux( 
      SHARED_PTR_CAST(clever::pdat::CleverSideData<double> ,
        patch.getPatchData(d_volflux, getCurrentDataContext())));

  const hier::Index ifirst = patch.getBox().lower();
  const hier::Index ilast = patch.getBox().upper();

  int xmin = ifirst(0); 
  int xmax = ilast(0); 
  int ymin = ifirst(1); 
  int ymax = ilast(1); 

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

  const std::shared_ptr<geom::CartesianPatchGeometry> pgeom(
      SHARED_PTR_CAST(geom::CartesianPatchGeometry,
        patch.getPatchGeometry()));

  const std::vector<hier::BoundaryBox>& edge_bdry 
    = pgeom->getCodimensionBoundaries(Bdry::EDGE2D);

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
        vol_flux_x = v_volflux->getPointer(0);
        vol_flux_y = v_volflux->getPointer(1);
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

        mass_flux_x = v_massflux->getPointer(0);
        mass_flux_y = v_massflux->getPointer(1);
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

        mass_flux_x = v_massflux->getPointer(0);
        mass_flux_y = v_massflux->getPointer(1);
      } break;
    default : tbox::perr << "[ERROR] Unknown exchange id in setPhysicalBoundaryConditions... " 
              << std::endl;
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
      default : tbox::perr << "[ERROR] Unknown edge location in setPhysicalBoundaryConditions... "
                << std::endl;
                exit(-1);
    }
  }

  t_fill_boundary->stop();
}

void Cleverleaf::field_summary(
    hier::Patch& patch,
    double* total_volume,
    double* total_mass,
    double* total_pressure,
    double* total_internal_energy,
    double* total_kinetic_energy,
    int* total_effective_cells)
{
  std::shared_ptr<clever::pdat::CleverCellData<double> > volume(
      SHARED_PTR_CAST(clever::pdat::CleverCellData<double> ,
        patch.getPatchData(d_volume, getCurrentDataContext())));

  std::shared_ptr<clever::pdat::CleverCellData<double> > density0(
      SHARED_PTR_CAST(clever::pdat::CleverCellData<double> ,
        patch.getPatchData(d_density, getCurrentDataContext())));

  std::shared_ptr<clever::pdat::CleverCellData<double> > energy0(
      SHARED_PTR_CAST(clever::pdat::CleverCellData<double> ,
        patch.getPatchData(d_energy, getCurrentDataContext())));

  std::shared_ptr<clever::pdat::CleverCellData<double> > pressure(
      SHARED_PTR_CAST(clever::pdat::CleverCellData<double> ,
        patch.getPatchData(d_pressure, getCurrentDataContext())));

  std::shared_ptr<clever::pdat::CleverNodeData<double> > velocity0(
      SHARED_PTR_CAST(clever::pdat::CleverNodeData<double> ,
        patch.getPatchData(d_velocity, getCurrentDataContext())));

  std::shared_ptr<clever::pdat::CleverCellData<int> > level_indicator(
      SHARED_PTR_CAST(clever::pdat::CleverCellData<int> ,
        patch.getPatchData(d_level_indicator, getCurrentDataContext())));

  const hier::Index ifirst = patch.getBox().lower();
  const hier::Index ilast = patch.getBox().upper();

  int xmin = ifirst(0);
  int xmax = ilast(0);
  int ymin = ifirst(1);
  int ymax = ilast(1);

  int level_number = patch.getPatchLevelNumber();

  F90_FUNC(field_summary_kernel, FIELD_SUMMARY_KERNEL)
    (&xmin,
     &xmax,
     &ymin,
     &ymax,
     volume->getPointer(),
     density0->getPointer(),
     energy0->getPointer(),
     pressure->getPointer(),
     velocity0->getPointer(0),
     velocity0->getPointer(1),
     level_indicator->getPointer(),
     total_volume,
     total_mass,
     total_internal_energy,
     total_kinetic_energy,
     total_pressure,
     total_effective_cells,
     &level_number);
}

void Cleverleaf::tagGradientDetectorCells(
    hier::Patch& patch,
    const double regrid_time,
    const bool initial_error,
    const int tag_index)
{
  std::shared_ptr<SAMRAI::pdat::CellData<int> > tags(
      SHARED_PTR_CAST(SAMRAI::pdat::CellData<int> ,
        patch.getPatchData(tag_index)));

  std::shared_ptr<clever::pdat::CleverCellData<double> > density0(
      SHARED_PTR_CAST(clever::pdat::CleverCellData<double> ,
        patch.getPatchData(d_density, getCurrentDataContext())));

  std::shared_ptr<clever::pdat::CleverCellData<double> > energy0(
      SHARED_PTR_CAST(clever::pdat::CleverCellData<double> ,
        patch.getPatchData(d_energy, getCurrentDataContext())));

  std::shared_ptr<clever::pdat::CleverCellData<double> > pressure(
      SHARED_PTR_CAST(clever::pdat::CleverCellData<double> ,
        patch.getPatchData(d_pressure, getCurrentDataContext())));

  std::shared_ptr<clever::pdat::CleverCellData<double> > viscosity(
      SHARED_PTR_CAST(clever::pdat::CleverCellData<double> ,
        patch.getPatchData(d_viscosity, getCurrentDataContext())));

  hier::Index ifirst = patch.getBox().lower();
  hier::Index ilast = patch.getBox().upper();

  int xmin = ifirst(0);
  int ymin = ifirst(1);
  int xmax = ilast(0);
  int ymax = ilast(1);

  hier::Box pbox = patch.getBox();
  hier::IntVector nghosts = density0->getGhostCellWidth();

  std::shared_ptr<clever::pdat::CleverCellData<int> > temporary_tags(
      new clever::pdat::CleverCellData<int>(pbox, 1, nghosts));

  temporary_tags->fillAll(0);
  tags->fillAll(0);

  if (d_tag_q) {
    F90_FUNC(tag_q_kernel,TAG_Q_KERNEL)
      (&xmin,
       &xmax,
       &ymin,
       &ymax,
       &d_tag_q_threshold,
       viscosity->getPointer(),
       temporary_tags->getPointer());
  }

  if (d_tag_density) {
    F90_FUNC(tag_density_kernel,TAG_DENSITY_KERNEL)
      (&xmin,
       &xmax,
       &ymin,
       &ymax,
       &d_tag_density_gradient,
       density0->getPointer(),
       temporary_tags->getPointer());
  }

  if (d_tag_energy) {
    F90_FUNC(tag_energy_kernel,TAG_ENERGY_KERNEL)
      (&xmin,
       &xmax,
       &ymin,
       &ymax,
       &d_tag_energy_gradient,
       energy0->getPointer(),
       temporary_tags->getPointer());
  }

  if (d_tag_pressure) {
    F90_FUNC(tag_pressure_kernel,TAG_PRESSURE_KERNEL)
      (&xmin,
       &xmax,
       &ymin,
       &ymax,
       &d_tag_pressure_gradient,
       pressure->getPointer(),
       temporary_tags->getPointer());
  }

  if(d_tag_all) {
    F90_FUNC(tag_all_kernel,TAG_ALL_KERNEL)
      (&xmin,
       &xmax,
       &ymin,
       &ymax,
       temporary_tags->getPointer());
  }

  int data_offset = density0->getGhostBox().offset(ifirst);
  int tag_data_offset = tags->getGhostBox().offset(ifirst);
  int tag_data_width = tags->getGhostBox().numberCells(0);

  for (int i1 = 0; i1 < patch.getBox().numberCells(1); i1++) {
    for (int i0 = 0; i0 < patch.getBox().numberCells(0); i0++) {
      int data_index = data_offset + i0;
      int tags_index = tag_data_offset + i0;

      tags->getPointer()[tags_index] = temporary_tags->getPointer()[data_index];
    }

    data_offset += density0->getGhostBox().numberCells(0);
    tag_data_offset += tag_data_width;
  }
}

void Cleverleaf::debug_knl(hier::Patch& patch)
{
  std::shared_ptr<clever::pdat::CleverCellData<double> > density0(
      SHARED_PTR_CAST(clever::pdat::CleverCellData<double> ,
        patch.getPatchData(d_density, getCurrentDataContext())));

  std::shared_ptr<clever::pdat::CleverCellData<double> > density1(
      SHARED_PTR_CAST(clever::pdat::CleverCellData<double> ,
        patch.getPatchData(d_density, getNewDataContext())));

  std::shared_ptr<clever::pdat::CleverCellData<double> > energy0(
      SHARED_PTR_CAST(clever::pdat::CleverCellData<double> ,
        patch.getPatchData(d_energy, getCurrentDataContext())));

  std::shared_ptr<clever::pdat::CleverCellData<double> > energy1(
      SHARED_PTR_CAST(clever::pdat::CleverCellData<double> ,
        patch.getPatchData(d_energy, getNewDataContext())));

  std::shared_ptr<clever::pdat::CleverCellData<double> > pressure(
      SHARED_PTR_CAST(clever::pdat::CleverCellData<double> ,
        patch.getPatchData(d_pressure, getCurrentDataContext())));

  std::shared_ptr<clever::pdat::CleverCellData<double> > soundspeed(
      SHARED_PTR_CAST(clever::pdat::CleverCellData<double> ,
        patch.getPatchData(d_soundspeed, getCurrentDataContext())));

  std::shared_ptr<clever::pdat::CleverCellData<double> > viscosity(
      SHARED_PTR_CAST(clever::pdat::CleverCellData<double> ,
        patch.getPatchData(d_viscosity, getCurrentDataContext())));

  std::shared_ptr<clever::pdat::CleverNodeData<double> > velocity0(
      SHARED_PTR_CAST(clever::pdat::CleverNodeData<double> ,
        patch.getPatchData(d_velocity, getCurrentDataContext())));

  std::shared_ptr<clever::pdat::CleverNodeData<double> > velocity1(
      SHARED_PTR_CAST(clever::pdat::CleverNodeData<double> ,
        patch.getPatchData(d_velocity, getNewDataContext())));

  std::shared_ptr<clever::pdat::CleverSideData<double> > volume_flux(
      SHARED_PTR_CAST(clever::pdat::CleverSideData<double> ,
        patch.getPatchData(d_volflux, getCurrentDataContext())));

  std::shared_ptr<clever::pdat::CleverSideData<double> > mass_flux(
      SHARED_PTR_CAST(clever::pdat::CleverSideData<double> ,
        patch.getPatchData(d_massflux, getCurrentDataContext())));

  std::shared_ptr<clever::pdat::CleverNodeData<double> > pre_volume(
      SHARED_PTR_CAST(clever::pdat::CleverNodeData<double> ,
        patch.getPatchData(d_workarray1, getCurrentDataContext())));

  std::shared_ptr<clever::pdat::CleverNodeData<double> > post_volume(
      SHARED_PTR_CAST(clever::pdat::CleverNodeData<double> ,
        patch.getPatchData(d_workarray2, getCurrentDataContext())));

  std::shared_ptr<clever::pdat::CleverNodeData<double> > pre_mass(
      SHARED_PTR_CAST(clever::pdat::CleverNodeData<double> ,
        patch.getPatchData(d_workarray3, getCurrentDataContext())));

  std::shared_ptr<clever::pdat::CleverNodeData<double> > post_mass(
      SHARED_PTR_CAST(clever::pdat::CleverNodeData<double> ,
        patch.getPatchData(d_workarray4, getCurrentDataContext())));

  std::shared_ptr<clever::pdat::CleverNodeData<double> > advected_volume(
      SHARED_PTR_CAST(clever::pdat::CleverNodeData<double> ,
        patch.getPatchData(d_workarray5, getCurrentDataContext())));

  std::shared_ptr<clever::pdat::CleverNodeData<double> > post_energy(
      SHARED_PTR_CAST(clever::pdat::CleverNodeData<double> ,
        patch.getPatchData(d_workarray6, getCurrentDataContext())));

  std::shared_ptr<clever::pdat::CleverNodeData<double> > energy_flux(
      SHARED_PTR_CAST(clever::pdat::CleverNodeData<double> ,
        patch.getPatchData(d_workarray7, getCurrentDataContext())));

  const hier::Index ifirst = patch.getBox().lower();
  const hier::Index ilast = patch.getBox().upper();

  int xmin = ifirst(0); 
  int xmax = ilast(0); 
  int ymin = ifirst(1); 
  int ymax = ilast(1); 

  F90_FUNC(debug_kernel, DEBUG_KERNEL)
    (&xmin,
     &xmax,
     &ymin,
     &ymax,
     density0->getPointer(),
     density1->getPointer(),
     energy0->getPointer(),
     energy1->getPointer(),
     pressure->getPointer(),
     soundspeed->getPointer(),
     viscosity->getPointer(),
     velocity0->getPointer(0),
     velocity0->getPointer(1),
     velocity1->getPointer(0),
     velocity1->getPointer(1),
     volume_flux->getPointer(0),
     volume_flux->getPointer(1),
     mass_flux->getPointer(0),
     mass_flux->getPointer(1),
     pre_volume->getPointer(),
     post_volume->getPointer(),
     pre_mass->getPointer(),
     post_mass->getPointer(),
     advected_volume->getPointer(),
     post_energy->getPointer(),
     energy_flux->getPointer());
}

void Cleverleaf::fillLevelIndicator(
    hier::Patch& patch,
    const int level_number)
{
    std::shared_ptr<clever::pdat::CleverCellData<int> > level_indicator(
        SHARED_PTR_CAST(clever::pdat::CleverCellData<int>,
          patch.getPatchData(d_level_indicator, getCurrentDataContext())));

    level_indicator->fillAll(level_number);
}

void Cleverleaf::initializeCallback()
{
  t_fill_boundary = tbox::TimerManager::getManager()->getTimer(
      "Cleverleaf::setPhysicalBoundaryConditions()");
}

void Cleverleaf::finalizeCallback()
{
  t_fill_boundary.reset();
}
