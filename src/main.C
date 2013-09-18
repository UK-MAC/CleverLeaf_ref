#include "SAMRAI/SAMRAI_config.h"

// Header for application-specific algorithm/data structure objects
#include "Cleverleaf.h"
#include "LagrangianEulerianLevelIntegrator.h"
#include "LagrangianEulerianIntegrator.h"

// Headers for SAMRAI
#include "SAMRAI/mesh/BergerRigoutsos.h"
#include "SAMRAI/geom/CartesianGridGeometry.h"
#include "SAMRAI/mesh/GriddingAlgorithm.h"
#include "SAMRAI/mesh/TreeLoadBalancer.h"
#include "SAMRAI/hier/PatchHierarchy.h"
#include "SAMRAI/hier/VariableDatabase.h"
#include "SAMRAI/mesh/StandardTagAndInitialize.h"
#include "SAMRAI/appu/VisItDataWriter.h"
#include "SAMRAI/tbox/InputManager.h"
#include "SAMRAI/tbox/PIO.h"

// Normal headers
#include <iostream>
#include <cstdlib>
#include <string>
#include <fstream>

using namespace std;
using namespace SAMRAI;

string getFilenameFromPath(const string path);
string createUniqueFilename(const string name, const string extension);
bool fileExists(const string path);

int main(int argc, char* argv[]) {

    tbox::SAMRAI_MPI::init(&argc, &argv);
    tbox::SAMRAIManager::initialize();
    tbox::SAMRAIManager::startup();

    const tbox::SAMRAI_MPI& mpi(tbox::SAMRAI_MPI::getSAMRAIWorld());

    {
        string input_path;

        if ((argc != 2)) {
            tbox::pout << "USAGE: " << argv[0] << " <input filename>"
                << endl;
        } else {
            input_path = argv[1];
        }

        tbox::plog << "Reading input from: " << input_path << endl;

        boost::shared_ptr<tbox::InputDatabase> input_db(new tbox::InputDatabase("input_db"));
        tbox::InputManager::getManager()->parseInputFile(input_path, input_db);

        boost::shared_ptr<tbox::Database> main_db = input_db->getDatabase("Cleverleaf");

        string input_file = getFilenameFromPath(input_path);
        string basename = main_db->getStringWithDefault("basename", input_file);

        const tbox::Dimension dim(static_cast<unsigned short>(main_db->getInteger("dim")));

        int vis_dump_interval = main_db->getIntegerWithDefault("vis_dump_interval", 1);
        int field_summary_interval = main_db->getIntegerWithDefault("field_summary_interval", 10);

        string log_basename = main_db->getStringWithDefault("log_filename", basename);
        string log_filename = createUniqueFilename(log_basename, "log");

        bool log_all_nodes = false;
        log_all_nodes = main_db->getBoolWithDefault("log_all_nodes", log_all_nodes);

        if (log_all_nodes) {
            tbox::PIO::logAllNodes(log_filename);
        } else {
            tbox::PIO::logOnlyNodeZero(log_filename);
        }

        boost::shared_ptr<geom::CartesianGridGeometry> grid_geometry(
                new geom::CartesianGridGeometry(dim,
                    "CartesianGeometry",
                    input_db->getDatabase("CartesianGeometry")));

        boost::shared_ptr<hier::PatchHierarchy> patch_hierarchy(
                new hier::PatchHierarchy("PatchHierarchy",
                    grid_geometry,
                    input_db->getDatabase("PatchHierarchy")));


        Cleverleaf* cleverleaf = new Cleverleaf(
                main_db,
                patch_hierarchy,
                dim,
                grid_geometry);

        boost::shared_ptr<LagrangianEulerianLevelIntegrator> lagrangian_eulerian_level_integrator(
                new LagrangianEulerianLevelIntegrator("LagrangianEulerianLevelIntegrator",
                    input_db->getDatabase("LagrangianEulerianLevelIntegrator"),
                    cleverleaf));

        boost::shared_ptr<mesh::StandardTagAndInitialize> error_detector(
                new mesh::StandardTagAndInitialize(
                    "StandardTagAndInitialize",
                    lagrangian_eulerian_level_integrator.get(),
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

        boost::shared_ptr<LagrangianEulerianIntegrator> lagrangian_eulerian_integrator(
                new LagrangianEulerianIntegrator(
                    input_db->getDatabase("LagrangianEulerianIntegrator"),
                    patch_hierarchy,
                    lagrangian_eulerian_level_integrator,
                    gridding_algorithm));

        int visit_number_procs_per_file = 1;
        std::string visit_dirname = createUniqueFilename(basename, "visit");

        boost::shared_ptr<appu::VisItDataWriter> visit_data_writer(
                new appu::VisItDataWriter(dim,
                    "Cleverleaf VisIt Writer",
                    visit_dirname,
                    visit_number_procs_per_file));

        cleverleaf->registerVisItDataWriter(visit_data_writer);

        double dt_now = lagrangian_eulerian_integrator->initializeHierarchy();

        if(vis_dump_interval > 0) {
            visit_data_writer->writePlotData(patch_hierarchy,
                    0,
                    0.0);
        }

        tbox::plog << "\nCheck input data and variables before simulation:"
            << endl;
        tbox::plog << "Input database..." << endl;
        input_db->printClassData(tbox::plog);
        tbox::plog << "\nVariable database..." << endl;
        hier::VariableDatabase::getDatabase()->printClassData(tbox::plog);

        double loop_time = lagrangian_eulerian_integrator->getIntegratorTime();
        double loop_time_end = lagrangian_eulerian_integrator->getEndTime();

        while ((loop_time < loop_time_end) &&
                lagrangian_eulerian_integrator->stepsRemaining()) {

            int iteration_num = lagrangian_eulerian_integrator->getIntegratorStep() + 1;

            tbox::pout << "++++++++++++++++++++++++++++++++++++++++++++" << endl;
            tbox::pout << "At begining of timestep # " << iteration_num - 1
                << endl;
            tbox::pout << "Simulation time is " << loop_time << endl;
            tbox::pout << "Current dt is " << dt_now << endl;

            double dt_new = lagrangian_eulerian_integrator->advanceHierarchy(dt_now);

            loop_time += dt_now;
            dt_now = dt_new;

            tbox::pout << "At end of timestep # " << iteration_num - 1 << endl;
            tbox::pout << "Simulation time is " << loop_time << endl;
            tbox::pout << "++++++++++++++++++++++++++++++++++++++++++++" << endl;

            if ((vis_dump_interval > 0)
                    && (iteration_num % vis_dump_interval) == 0) {
                visit_data_writer->writePlotData(patch_hierarchy,
                        iteration_num,
                        loop_time);
            }

            if ((field_summary_interval > 0)
                    && (iteration_num % field_summary_interval) == 0) {
            }
        }

        if (vis_dump_interval > 0) {
            visit_data_writer->writePlotData(patch_hierarchy,
                    lagrangian_eulerian_integrator->getIntegratorStep() + 1,
                    loop_time);
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

string getFilenameFromPath(const string path)
{
    std::cout << "path: " << path << std::endl;

    size_t start = path.find_last_of("/");
    start += 1;

    size_t end = path.find_last_of(".");

    if(start == string::npos)
        start = 0;

    if(end == string::npos)
        end = path.length();

    size_t count = end - start;

    std::cout << "start, end, count: " << start << ", " << end << ", " << count << std::endl;

    return path.substr(start, count);
}

string createUniqueFilename(const string name, const string extension)
{
    int unique_id = -1;
    stringstream ss;
    string filename;

    do {
        ss.str("");
        ss.clear();
        unique_id++;
        ss << name << "." << unique_id << "." << extension;
        filename = ss.str();
    } while (fileExists(filename));

    return filename;
}

bool fileExists(const string path)
{
    ifstream ifile(path.c_str());

    if(ifile) {
        return true;
    } else {
        return false;
    }
}
