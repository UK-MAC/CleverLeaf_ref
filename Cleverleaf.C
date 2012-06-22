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

    d_celldeltas = new pdat::CellVariable<double>(d_dim, "celldelta", d_dim.getValue());
    d_cellcoords = new pdat::CellVariable<double>(d_dim, "cellcoords", d_dim.getValue());

    d_vertexdeltas = new pdat::NodeVariable<double>(d_dim, "vertexdeltas", d_dim.getValue());
    d_vertexcoords = new pdat::NodeVariable<double>(d_dim, "vertexcoords", d_dim.getValue());
}

void Cleverleaf::registerModelVariables(
        LagrangianEulerianIntegrator* integrator) 
{
    // TODO: Change velocity back to field variable!
    integrator->registerVariable(d_velocity, LagrangianEulerianIntegrator::NORMAL,  d_nghosts, d_grid_geometry);
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

        d_visit_writer->registerPlotQuantity("Massflux",
                "VECTOR",
                vardb->mapVariableAndContextToIndex(
                    d_massflux, d_plot_context));

        d_visit_writer->registerPlotQuantity("Volflux",
                "VECTOR",
                vardb->mapVariableAndContextToIndex(
                    d_volflux, d_plot_context));

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
        tbox::Pointer<pdat::CellData<double> > velocity = patch.getPatchData(d_velocity, getCurrentDataContext());
        tbox::Pointer<pdat::CellData<double> > massflux = patch.getPatchData(d_massflux, getCurrentDataContext());
        tbox::Pointer<pdat::CellData<double> > volflux = patch.getPatchData(d_volflux, getCurrentDataContext());
        tbox::Pointer<pdat::CellData<double> > pressure = patch.getPatchData(d_pressure, getCurrentDataContext());
        tbox::Pointer<pdat::CellData<double> > viscosity = patch.getPatchData(d_viscosity, getCurrentDataContext());
        tbox::Pointer<pdat::CellData<double> > soundspeed = patch.getPatchData(d_soundspeed, getCurrentDataContext());
        tbox::Pointer<pdat::CellData<double> > density = patch.getPatchData(d_density, getCurrentDataContext());
        tbox::Pointer<pdat::CellData<double> > energy = patch.getPatchData(d_energy, getCurrentDataContext());
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

                if (((i >= imin + 32) && (i <= imax - 32)) &&
                        ((j >= jmin + 32) && ( j <= jmax - 32))) {
                    density_data[n1] = 1.0;
                    energy_data[n1] = 2.5;
                } else {
                    density_data[n1] = 0.1;
                    energy_data[n1] = 1.0;
                }
            }
        }

        /*
         * Fill in the volume array...
         */
        const tbox::Pointer<geom::CartesianPatchGeometry> pgeom = patch.getPatchGeometry();
        const double* dxs = pgeom->getDx();
        const double* coords = pgeom->getXLower();
        double xmin = coords[0];
        double ymin = coords[1];
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

                vertexx[ind] = xmin + dx*(xcount-xmin);
                vertexy[ind] = ymin + dy*(ycount-ymin);

                xcount++;
            }

            ycount++;
        }

        double* cellx = cellcoords->getPointer(0);
        double* celly = cellcoords->getPointer(1);

        for(int j = jmin; j <= jmax; j++) {
            for(int i = imin; i <= imax; i++) {
                int vind = POLY2(i,j,vimin,vjmin,vnx);
                int vind2 = POLY2(i+1,j,vimin,vjmin,vnx);
                int vind3 = POLY2(i,j+1,vimin,vjmin,vnx);

                int ind = POLY2(i,j,imin,jmin,nx);

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
}

void Cleverleaf::ideal_gas_knl(
        hier::Patch& patch)
{

    const hier::Index ifirst = patch.getBox().lower();
    const hier::Index ilast = patch.getBox().upper();

    int xmin = ifirst(0); 
    int xmax = ilast(0); 
    int ymin = ifirst(1); 
    int ymax = ilast(1); 

    tbox::Pointer<pdat::CellData<double> > v_density0 = patch.getPatchData(d_density, getCurrentDataContext());

    tbox::Pointer<pdat::CellData<double> > v_energy0 = patch.getPatchData(d_energy, getCurrentDataContext());

    tbox::Pointer<pdat::CellData<double> > v_pressure = patch.getPatchData(d_pressure, getCurrentDataContext());

    tbox::Pointer<pdat::CellData<double> > v_soundspeed = patch.getPatchData(d_soundspeed, getCurrentDataContext());

    double* density = v_density0->getPointer();
    double* energy = v_energy0->getPointer();
    double* pressure = v_pressure->getPointer();
    double* soundspeed = v_soundspeed->getPointer();

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

    tbox::Pointer<pdat::CellData<double> > v_density0 = patch.getPatchData(d_density, getCurrentDataContext());

    tbox::Pointer<pdat::CellData<double> > v_pressure = patch.getPatchData(d_pressure, getCurrentDataContext());

    tbox::Pointer<pdat::CellData<double> > v_viscosity = patch.getPatchData(d_viscosity, getCurrentDataContext());
    tbox::Pointer<pdat::CellData<double> > v_vel0 = patch.getPatchData(d_velocity, getCurrentDataContext());
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

          int n1 = POLY2(j,k,xmin,ymin, (xmax-xmin+1));
          int n2 = POLY2(j+1,k,xmin,ymin, (xmax-xmin+1));
          int n3 = POLY2(j,k-1,xmin,ymin, (xmax-xmin+1));
          int n4 = POLY2(j-1,k,xmin,ymin, (xmax-xmin+1));
          int n5 = POLY2(j,k+1,xmin,ymin, (xmax-xmin+1));
          int n6 = POLY2(j+1,k+1,xmin,ymin, (xmax-xmin+1));
          int n7 = POLY2(j+1,k-1,xmin,ymin, (xmax-xmin+1));

          ugrad=(xvel0[n2]+xvel0[n6])-(xvel0[n1]+xvel0[n5]);

          vgrad=(yvel0[n5]+yvel0[n6])-(yvel0[n1]+yvel0[n2]);

          div=(celldx[n1]*(ugrad)+celldy[n1]*(vgrad));
          strain2=0.5*(xvel0[n5]+xvel0[n6]-xvel0[n1]-xvel0[n2])/celldy[n1]+0.5*(yvel0[n2]+yvel0[n6]-yvel0[n1]-yvel0[n5])/celldx[n1];

          pgradx=(pressure[n2]-pressure[n4])/(celldx[n1]+celldx[n2]);
          pgrady=(pressure[n5]-pressure[n3])/(celldy[n1]+celldy[n5]);

          pgradx2 = pgradx*pgradx;
          pgrady2 = pgrady*pgrady;

          limiter = ((0.5*(ugrad)/celldx[n1])*pgradx2+(0.5*(vgrad)/celldy[n1])*pgrady2+strain2*pgradx*pgrady)/max(pgradx2+pgrady2,1.0e-16);

          pgradx = copysign(max(1.0e-16,abs(pgradx)),pgradx);
          pgrady = copysign(max(1.0e-16,abs(pgrady)),pgrady);
          pgrad = sqrt((pgradx*pgradx)+(pgrady*pgrady));
          xgrad = abs(celldx[n1]*pgrad/pgradx);
          ygrad = abs(celldy[n1]*pgrad/pgrady);
          grad  = min(xgrad,ygrad);
          grad2 = grad*grad;

#ifdef DEBUG
          tbox::pout << "limiter = " << limiter << ",div = " << div << std::endl;
#endif
          if(!((limiter > 0.0) || (div >= 0.0))) {
#ifdef DEBUG
            tbox::pout << "Updating viscosity[" << j-xmin << "][" << k-ymin << "]=" << 2.0*density[n1]*grad2*(limiter*limiter) << std::endl;
#endif
              viscosity[n1]=2.0*density[n1]*grad2*(limiter*limiter);
          } else {
#ifdef DEBUG
            tbox::pout << "[ZERO] Updating viscosity[" << j-xmin << "][" << k-ymin << "]=" << 0.0 << std::endl;
#endif
              viscosity[n1] = 0.0;
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
    tbox::Pointer<pdat::CellData<double> > v_vel0 = patch.getPatchData(d_velocity, getCurrentDataContext());
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
    tbox::Pointer<pdat::CellData<double> > v0 = patch.getPatchData(d_velocity, getCurrentDataContext());
    tbox::Pointer<pdat::CellData<double> > v1 = patch.getPatchData(d_velocity, getNewDataContext());

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

      if ( predict ) {

          for (int k = ymin; k <= ymax; k++) {
              for (int j = xmin; j <= xmax; j++) { 

                  int n1 = POLY2(j,k,xmin,ymin, (xmax-xmin+1));
                  int n2 = POLY2(j+1,k,xmin,ymin, (xmax-xmin+1));
                  int n3 = POLY2(j,k-1,xmin,ymin, (xmax-xmin+1));
                  int n4 = POLY2(j-1,k,xmin,ymin, (xmax-xmin+1));
                  int n5 = POLY2(j,k+1,xmin,ymin, (xmax-xmin+1));
                  int n6 = POLY2(j+1,k+1,xmin,ymin, (xmax-xmin+1));
                  int n7 = POLY2(j+1,k-1,xmin,ymin, (xmax-xmin+1));

                  left_flux=(xarea[n1]*(xvel0[n1]+xvel0[n5]
                              +xvel0[n1]+xvel0[n5]))*0.25*dt*0.5;

                  right_flux=(xarea[n2]*(xvel0[n2]+xvel0[n6]
                              +xvel0[n2]+xvel0[n6]))*0.25*dt*0.5;

                  bottom_flux=(yarea[n1]*(yvel0[n1]+yvel0[n2]
                              +yvel0[n1]+yvel0[n2]))*0.25*dt*0.5;

                  top_flux=(yarea[n5]*(yvel0[n5]+yvel0[n6]
                              +yvel0[n5]+yvel0[n6]))*0.25*dt*0.5;

                  total_flux=right_flux-left_flux+top_flux-bottom_flux;

                  volume_change[POLY2(j-xmin,k-ymin,0,0,(xmax-xmin+1))]=volume[n1]/(volume[n1]+total_flux);

                  min_cell_volume=min(
                          volume[n1]+right_flux-left_flux+top_flux-bottom_flux,
                          min(volume[n1]+right_flux-left_flux,
                          volume[n1]+top_flux-bottom_flux));

                  //Perhaps take these tests out so it will vectorise;
                  if (volume_change[POLY2(j-xmin,k-ymin,0,0,(xmax-xmin+1))] <= 0.0 ) { 
                      cout << "HIT PDV ERROR COND 1 " << endl;
                  }
                  if (min_cell_volume <= 0.0) {
                      cout << "HIT PDV ERROR COND 2" << endl;
                  }

                  recip_volume=1.0/volume[n1];

#ifdef DEBUG
                  cout << pressure[n1] << ", " << density0[n1] << ", " << viscosity[n1] << endl;
#endif
                  energy_change=(pressure[n1]/density0[n1]+viscosity[n1]/density0[n1])*total_flux*recip_volume;

#ifdef DEBUG
                  cout << "total flux=" << total_flux << ", n1=" << n1 << ", recip_volume=" << recip_volume << endl;
                  cout << "energy change=" << energy_change << ", volume_change=" << volume_change[POLY2(j-xmin,k-ymin,0,0,(xmax-xmin+1))] << endl;
#endif
                  energy1[n1]=energy0[n1]-energy_change;

                  density1[n1]=density0[n1]*volume_change[POLY2(j-xmin,k-ymin,0,0,(xmax-xmin+1))];
              }
          }

  } else {

//    DO k=y_min,y_max
//      DO j=x_min,x_max
//
//        left_flux=  (xarea(j  ,k  )*(xvel0(j  ,k  )+xvel0(j  ,k+1)                     &
//                                    +xvel1(j  ,k  )+xvel1(j  ,k+1)))*0.25_8*dt
//        right_flux= (xarea(j+1,k  )*(xvel0(j+1,k  )+xvel0(j+1,k+1)                     &
//                                    +xvel1(j+1,k  )+xvel1(j+1,k+1)))*0.25_8*dt
//        bottom_flux=(yarea(j  ,k  )*(yvel0(j  ,k  )+yvel0(j+1,k  )                     &
//                                    +yvel1(j  ,k  )+yvel1(j+1,k  )))*0.25_8*dt
//        top_flux=   (yarea(j  ,k+1)*(yvel0(j  ,k+1)+yvel0(j+1,k+1)                     &
//                                    +yvel1(j  ,k+1)+yvel1(j+1,k+1)))*0.25_8*dt
//        total_flux=right_flux-left_flux+top_flux-bottom_flux
//
//        volume_change(j,k)=volume(j,k)/(volume(j,k)+total_flux)
//
//        min_cell_volume=MIN(volume(j,k)+right_flux-left_flux+top_flux-bottom_flux &
//                           ,volume(j,k)+right_flux-left_flux                      &
//                           ,volume(j,k)+top_flux-bottom_flux)
// 
//        IF(volume_change(j,k).LE.0.0) THEN ! Perhaps take these tests out so it will vectorise
//          error_condition=1 ! Do I need atomic for OpenMP?
//        ENDIF
//        IF(min_cell_volume.LE.0.0) THEN
//          error_condition=2
//        ENDIF
//
//        recip_volume=1.0/volume(j,k) 
//
//        energy_change=(pressure(j,k)/density0(j,k)+viscosity(j,k)/density0(j,k))*total_flux*recip_volume
//
//        energy1(j,k)=energy0(j,k)-energy_change
//
//        density1(j,k)=density0(j,k)*volume_change(j,k)
//
//        density1(j,k)=density0(j,k)*volume_change(j,k)
//
//      ENDDO
//    ENDDO

  }

      delete volume_change;
}
