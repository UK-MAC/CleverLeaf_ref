#include "SAMRAI/SAMRAI_config.h"

// Header for application-specific algorithm/data structure objects
#include "Cleverleaf.h"
#include "LagrangianEulerianIntegrator.h"

// Headers for SAMRAI
#include "SAMRAI/algs/TimeRefinementIntegrator.h"
#include "SAMRAI/mesh/BergerRigoutsos.h"
#include "SAMRAI/geom/CartesianGridGeometry.h"
#include "SAMRAI/mesh/GriddingAlgorithm.h"
#include "SAMRAI/mesh/TreeLoadBalancer.h"
#include "SAMRAI/hier/PatchHierarchy.h"
#include "SAMRAI/hier/VariableDatabase.h"
#include "SAMRAI/mesh/StandardTagAndInitialize.h"
#include "SAMRAI/appu/VisItDataWriter.h"
#include "SAMRAI/tbox/InputManager.h"

// Normal headers
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

        boost::shared_ptr<tbox::InputDatabase> input_db(new tbox::InputDatabase("input_db"));
        tbox::InputManager::getManager()->parseInputFile(input_filename, input_db);

        /*
         * Read in main data.
         */
        boost::shared_ptr<tbox::Database> main_db = input_db->getDatabase("Main");

        const tbox::Dimension dim(static_cast<unsigned short>(main_db->getInteger("dim")));

        int vis_dump_interval = main_db->getIntegerWithDefault("vis_dump_interval", 1);
        int field_summary_interval = main_db->getIntegerWithDefault("field_summary_interval", 10);

        /*
         * Create data and algorithm objects.
         */
        boost::shared_ptr<geom::CartesianGridGeometry> grid_geometry(
                new geom::CartesianGridGeometry(dim,
                    "CartesianGeometry",
                    input_db->getDatabase("CartesianGeometry")));

        boost::shared_ptr<hier::PatchHierarchy> patch_hierarchy(
                new hier::PatchHierarchy("PatchHierarchy",
                    /* 
                     * Pass the grid geometry into the patch hierarchy, 
                     * following the Strategy pattern.
                     */
                    grid_geometry,
                    input_db->getDatabase("PatchHierarchy")));

        /*
         * Create the Cleverleaf model here to control the maths specific to this
         * program.
         */
        boost::shared_ptr<tbox::Database> cleverleaf_db = input_db->getDatabase("Cleverleaf");

        bool vis_me = cleverleaf_db->getBool("vis");
        int visit_number_procs_per_file = 1;
        const std::string visit_dump_dirname = "cleverleaf.visit";

        Cleverleaf* cleverleaf = new Cleverleaf(
                patch_hierarchy,
                dim,
                grid_geometry);

        boost::shared_ptr<LagrangianEulerianIntegrator> lagrangian_eulerian_integrator(
                new LagrangianEulerianIntegrator("LagrangianEulerianIntegrator",
                    input_db->getDatabase("LagrangianEulerianIntegrator"),
                    /*
                     * Pass the Cleverleaf model to the integrator,
                     * again following the Strategy pattern.
                     */
                    cleverleaf));

        boost::shared_ptr<mesh::StandardTagAndInitialize> error_detector(
                new mesh::StandardTagAndInitialize(dim,
                    "StandardTagAndInitialize",
                    lagrangian_eulerian_integrator.get(),
                    input_db->getDatabase("StandardTagAndInitialize")));

        boost::shared_ptr<mesh::BergerRigoutsos> box_generator(
                new mesh::BergerRigoutsos(
                    dim,
                    input_db->getDatabaseWithDefault(
                        "BergerRigoutsos",
                        boost::shared_ptr<SAMRAI::tbox::Database>())));

        boost::shared_ptr<mesh::TreeLoadBalancer> load_balancer(
                new mesh::TreeLoadBalancer(dim,
                    "LoadBalancer",
                    input_db->getDatabase("LoadBalancer")));

        load_balancer->setSAMRAI_MPI(
                SAMRAI::tbox::SAMRAI_MPI::getSAMRAIWorld());

        boost::shared_ptr<mesh::GriddingAlgorithm> gridding_algorithm(
                new mesh::GriddingAlgorithm(
                    patch_hierarchy,
                    "GriddingAlgorithm",
                    input_db->getDatabase("GriddingAlgorithm"),
                    error_detector,
                    box_generator,
                    load_balancer));

        boost::shared_ptr<algs::TimeRefinementIntegrator> time_integrator(
                new algs::TimeRefinementIntegrator("TimeRefinementIntegrator",
                    input_db->getDatabase("TimeRefinementIntegrator"),
                    patch_hierarchy,
                    lagrangian_eulerian_integrator,
                    gridding_algorithm));

        /*
         * All the SAMRAI components are now created!
         */

        /*
         * Set up ViSiT writer here.
         */
        boost::shared_ptr<appu::VisItDataWriter> visit_data_writer(
                new appu::VisItDataWriter(dim,
                    "Kamra VisIt Writer",
                    visit_dump_dirname,
                    visit_number_procs_per_file));

        cleverleaf->registerVisItDataWriter(visit_data_writer);

        /*
         * Initialise the hierarchy config and data on all patches
         */
        double dt_now = time_integrator->initializeHierarchy();

        visit_data_writer->writePlotData(patch_hierarchy,
                0,
                0.0);

        /*
         * After creating all objects and initializing their state, we
         * print the input database and variable database contents
         * to the log file.
         */
        tbox::plog << "\nCheck input data and variables before simulation:"
            << endl;
        tbox::plog << "Input database..." << endl;
        input_db->printClassData(tbox::plog);
        tbox::plog << "\nVariable database..." << endl;
        hier::VariableDatabase::getDatabase()->printClassData(tbox::plog);

        double loop_time = time_integrator->getIntegratorTime();
        double loop_time_end = time_integrator->getEndTime();

        /*
         * MAIN LOOP
         *
         * Step count and integration time are maintained by algs::TimeRefinementIntegrator
         */
        while ((loop_time < loop_time_end) &&
                time_integrator->stepsRemaining()) {

            int iteration_num = time_integrator->getIntegratorStep() + 1;

            tbox::pout << "++++++++++++++++++++++++++++++++++++++++++++" << endl;
            tbox::pout << "At begining of timestep # " << iteration_num - 1
                << endl;
            tbox::pout << "Simulation time is " << loop_time << endl;
            tbox::pout << "Current dt is " << dt_now << endl;

            double dt_new = time_integrator->advanceHierarchy(dt_now);

            loop_time += dt_now;
            dt_now = dt_new;

            tbox::pout << "At end of timestep # " << iteration_num - 1 << endl;
            tbox::pout << "Simulation time is " << loop_time << endl;
            tbox::pout << "++++++++++++++++++++++++++++++++++++++++++++" << endl;

            /*
             * Write out visualisation data
             */
            if ((vis_dump_interval > 0)
                    && (iteration_num % vis_dump_interval) == 0) {
                visit_data_writer->writePlotData(patch_hierarchy,
                        iteration_num,
                        loop_time);
            }

            if ((field_summary_interval > 0)
                    && (iteration_num % field_summary_interval) == 0) {
                tbox::pout << "Writing field summary..." << endl;
                lagrangian_eulerian_integrator->printFieldSummary(loop_time, iteration_num);
            }
        }

        /*
         * Deallocate objects
         */
        patch_hierarchy.reset();
        grid_geometry.reset();
        box_generator.reset();
        load_balancer.reset();
        error_detector.reset();
        gridding_algorithm.reset();

        if (cleverleaf) delete cleverleaf;

    }

    tbox::SAMRAIManager::shutdown();
    tbox::SAMRAIManager::finalize();
    tbox::SAMRAI_MPI::finalize();
}
