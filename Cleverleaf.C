#include "Cleverleaf.h"

void Cleverleaf::Cleverleaf(tbox::Pointer<hier::PatchHierarchy> hierarchy) {

    d_hierarchy = hierarchy;

    /*
     * Register variables
     */
    d_velocity0 = new pdat::CellVariable<double>(d_dim, "velocity", d_dim.getValue());
    d_massflux  = new pdat::CellVariable<double>(d_dim, "massflux", d_dim.getValue());
    d_volflux   = new pdat::CellVariable<double>(d_dim, "volflux", d_dim.getValue());
    d_pressure  = new pdat::CellVariable<double>(d_dim, "pressure", 1);
    d_viscosity  = new pdat::CellVariable<double>(d_dim, "viscosity", 1);
    d_soundspeed  = new pdat::CellVariable<double>(d_dim, "soundspeed", 1);
    d_density  = new pdat::CellVariable<double>(d_dim, "density", 1);
    d_energy  = new pdat::CellVariable<double>(d_dim, "energy", 1);

}

void Cleverleaf::registerModelVariables() {
    hier::VariableDatabase* var_db = hier::VariableDatabase.getDatabase();

    tbox::Pointer<hier::VariableContext> current = var_db->getContext("CURRENT");

     
}

void Cleverleaf::registerVisItDataWriter(tbox::Pointer<appu::VisItDataWriter> writer) {
    d_visit_writer = writer;
}
