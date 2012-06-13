#include "SAMRAI/SAMRAI_config.h"

// Header for application-specific algorithm/data structure object

#include "Clamour.h"

// Headers for major algorithm/data structure objects from SAMRAI

#include "SAMRAI/mesh/BergerRigoutsos.h"
#include "SAMRAI/geom/CartesianGridGeometry.h"
#include "SAMRAI/mesh/GriddingAlgorithm.h"
#include "SAMRAI/mesh/TreeLoadBalancer.h"
#include "SAMRAI/hier/PatchHierarchy.h"
#include "SAMRAI/mesh/StandardTagAndInitialize.h"
#include "SAMRAI/appu/VisItDataWriter.h"

// Headers for basic SAMRAI objects

#include "SAMRAI/tbox/Pointer.h"
#include "SAMRAI/tbox/InputManager.h"


#include <stdio.h>
#include <stdlib.h>
#include <string>
#include <fstream>

using namespace std;
using namespace SAMRAI;

int main(int argc, char* argv[]) {

    // Initialise MPI
    tbox::SAMRAI_MPI::init(&argc, &argv);
    tbox::SAMRAIManager::initialize();
    tbox::SAMRAIManager::startup();

    const tbox::SAMRAI_MPI& mpi(tbox::SAMRAI_MPI::getSAMRAIWorld());

    // Everything wrapped in a block so that the SAMRAI pointers get deallocated properly.
    {
        string input_filename;

        // Error check arguments
        if ((argc != 2)) {
            tbox::pout << "USAGE: " << argv[0] << " <input filename>"
                << endl;
        } else {
            input_filename = argv[1];
        }

        tbox::plog << "input filename = " << input_filename << endl;

        /*
         * Create input database and parse input file.
         */

        tbox::Pointer<tbox::Database> input_db(new tbox::InputDatabase("input_db"));
        tbox::InputManager::getManager()->parseInputFile(input_filename, input_db);

        /*
         * Read in main data.
         */
        tbox::Pointer<tbox::Database> main_db = input_db->getDatabase("Main");

        const tbox::Dimension dim(static_cast<unsigned short>(main_db->getInteger("dim")));

        /*
         * Create data and algorithm objects.
         */
        tbox::Pointer<geom::CartesianGridGeometry> grid_geometry(
                new geom::CartesianGridGeometry(dim,
                    "CartesianGeometry",
                    input_db->getDatabase("CartesianGeometry")));

        tbox::Pointer<hier::PatchHierarchy> patch_hierarchy(
                new hier::PatchHierarchy("PatchHierarchy",
                    /* 
                     * Pass the grid geometry into the patch hierarchy, 
                     * following the Strategy pattern.
                     */
                    grid_geometry,
                    input_db->getDatabase("PatchHierarchy")));

        /*
         * Create the kamra model here to control the maths specific to this
         * program.
         */
        tbox::Pointer<tbox::Database> clamour_db = input_db->getDatabase("Clamour");

        bool vis_me = clamour_db->getBool("vis");
        int visit_number_procs_per_file = 1;
        const std::string visit_dump_dirname = "clamour.visit";

        Clamour* clamour = new Clamour(patch_hierarchy);
        clamour->registerModelVariables();

        tbox::Pointer<mesh::StandardTagAndInitialize> error_detector(
                new mesh::StandardTagAndInitialize(dim,
                    "StandardTagAndInitialize",
                    clamour,
                    input_db->getDatabase("StandardTagAndInitialize")));

        tbox::Pointer<mesh::BergerRigoutsos> box_generator(
                new mesh::BergerRigoutsos(
                    dim,
                    input_db->getDatabaseWithDefault(
                        "BergerRigoutsos",
                        SAMRAI::tbox::Pointer<SAMRAI::tbox::Database>(NULL))));

        tbox::Pointer<mesh::TreeLoadBalancer> load_balancer(
                new mesh::TreeLoadBalancer(dim,
                    "LoadBalancer",
                    input_db->getDatabase("LoadBalancer")));

        load_balancer->setSAMRAI_MPI(
                SAMRAI::tbox::SAMRAI_MPI::getSAMRAIWorld());

        tbox::Pointer<mesh::GriddingAlgorithm> gridding_algorithm(
                new mesh::GriddingAlgorithm(
                    patch_hierarchy,
                    "GriddingAlgorithm",
                    input_db->getDatabase("GriddingAlgorithm"),
                    error_detector,
                    box_generator,
                    load_balancer));

        /*
         * All the SAMRAI components are now created!
         */

        /*
         * TODO: Set up ViSiT writer here.
         */

        tbox::Pointer<appu::VisItDataWriter> visit_data_writer(
                new appu::VisItDataWriter(dim,
                    "Kamra VisIt Writer",
                    visit_dump_dirname,
                    visit_number_procs_per_file));

        clamour->registerVisItDataWriter(visit_data_writer);

        /*
         * Write data file
         */
        //if(vis_me)
        //    visit_data_writer->writePlotData(
        //            patch_hierarchy,
        //            time_integrator->getIntegratorStep(),
        //            time_integrator->getIntegratorTime());

        /*
         * Initialise the hierarchy config and data on all patches
         */

        /*
         * MAIN LOOP
         *
         * Step count and integration time are maintained by algs::TimeRefinementIntegrator
         */

        /*
         * Deallocate objects
         */

        patch_hierarchy.setNull();
        grid_geometry.setNull();
        box_generator.setNull();
        load_balancer.setNull();
        error_detector.setNull();
        gridding_algorithm.setNull();

        if (clamour) delete clamour;

    }

    tbox::SAMRAIManager::shutdown();
    tbox::SAMRAIManager::finalize();
    tbox::SAMRAI_MPI::finalize();
}
