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
#include <iostream>
#include <cstdlib>
#include <string>
#include <fstream>
#include <unistd.h>

#if defined(_OPENMP)
#include <omp.h>
#endif

#include "hydro/Cleverleaf.h"
#include "hydro/LagrangianEulerianLevelIntegrator.h"
#include "hydro/LagrangianEulerianIntegrator.h"

#include "SAMRAI/mesh/BergerRigoutsos.h"
#include "SAMRAI/geom/CartesianGridGeometry.h"
#include "SAMRAI/mesh/GriddingAlgorithm.h"
#include "SAMRAI/mesh/ChopAndPackLoadBalancer.h"
#include "SAMRAI/mesh/TreeLoadBalancer.h"
#include "SAMRAI/hier/PatchHierarchy.h"
#include "SAMRAI/hier/VariableDatabase.h"
#include "SAMRAI/mesh/StandardTagAndInitialize.h"
#include "SAMRAI/appu/VisItDataWriter.h"
#include "SAMRAI/tbox/InputManager.h"
#include "SAMRAI/tbox/PIO.h"
#include "SAMRAI/tbox/TimerManager.h"

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
    if (argc == 1) {
      tbox::pout << "USAGE: " << argv[0] 
        << " [-l max_level] [-r refinement_ratio] [-f regrid_frequency]" 
        << " <input filename>" << std::endl;
      exit(1);
    }

    string input_path;

    int max_level_number = -1;
    int refinement_ratio = -1;
    int regrid_interval = -1;

    int opt;

    while ((opt = getopt(argc,argv,"l:r:f:i:")) != EOF ) {
      switch (opt) {
        case 'l':
          max_level_number = atoi(optarg);
          tbox::plog << "max_level_number overide using command-line: " 
            << max_level_number << endl;
          break;
        case 'r':
          refinement_ratio = atoi(optarg);
          tbox::plog << "refinement_ratio overide using command-line: " 
            << refinement_ratio << endl;
          break;
        case 'f':
          regrid_interval = atoi(optarg);
          tbox::plog << "regrid_interval overide using command-line: " 
            << regrid_interval << endl;
          break;
        case 'i':
          input_path = optarg;
          break;
        case '?':
          if (optopt == 'i')
            tbox::perr << "Option -i requires input file argument" << std::endl;
          else
            tbox::perr << "Unknown option " << optopt << std::endl;

          tbox::pout << "USAGE: " << argv[0] 
            << " [-l max_level] [-r refinement_ratio] [-f regrid_frequency]" 
            << " <input filename>" << std::endl;

          exit(1);
      }
    }

    if (input_path.empty() && optind < argc) {
      input_path = argv[optind];
    } else {
      tbox::perr << "Input file required" << std::endl;
      tbox::pout << "USAGE: " << argv[0] 
        << " [-l max_level] [-r refinement_ratio] [-f regrid_frequency]" 
        << " <input filename>" << std::endl;
      exit(1);
    }

    tbox::pout << "CleverLeaf version #" << VERSION 
      << " compiled on " << HOST_NAME << std::endl;
    tbox::pout << "Running with " << mpi.getSize() << " tasks";
#if defined(_OPENMP)
#pragma omp parallel
    {
#pragma omp master
      { tbox::pout << " and " << omp_get_num_threads() << " threads"; }
    }
#endif
    tbox::pout << std::endl;

    tbox::plog << "Reading input from: " << input_path << endl;

    boost::shared_ptr<tbox::InputDatabase> input_db(
        new tbox::InputDatabase("input_db"));
    tbox::InputManager::getManager()->parseInputFile(input_path, input_db);

    boost::shared_ptr<tbox::Database> main_db = input_db->getDatabase(
        "Cleverleaf");

    string input_file = getFilenameFromPath(input_path);
    string basename = main_db->getStringWithDefault("basename", input_file);

    const tbox::Dimension dim(
        static_cast<unsigned short>(main_db->getInteger("dim")));

    int vis_dump_interval = main_db->getIntegerWithDefault(
        "vis_dump_interval", 1);
    int field_summary_interval = main_db->getIntegerWithDefault(
        "field_summary_interval", 10);

    string log_basename = main_db->getStringWithDefault(
        "log_filename", basename);
    string log_filename = createUniqueFilename(log_basename, "log");

    bool log_all_nodes = false;
    log_all_nodes = main_db->getBoolWithDefault("log_all_nodes", log_all_nodes);

    if (log_all_nodes) {
      tbox::PIO::logAllNodes(log_filename);
    } else {
      tbox::PIO::logOnlyNodeZero(log_filename);
    }

    if(input_db->isDatabase("TimerManager")) {
      tbox::TimerManager::createManager(input_db->getDatabase("TimerManager"));
    }

    boost::shared_ptr<geom::CartesianGridGeometry> grid_geometry(
        new geom::CartesianGridGeometry(
          dim,
          "CartesianGeometry",
          input_db->getDatabase("CartesianGeometry")));

    /*
     * Overwrite max_levels and ratio_to_coarser if passed as command line
     * arguments.
     */
    if(max_level_number != -1) {
      input_db->getDatabase("PatchHierarchy")
        ->putInteger("max_levels", max_level_number);
    }

    if(refinement_ratio != -1) {
      std::vector<int> ratio_vector(dim.getValue());

      for(int i = 0; i < dim.getValue(); i++) {
        ratio_vector[i] = refinement_ratio;
      }

      input_db->getDatabase("PatchHierarchy")
        ->getDatabase("ratio_to_coarser")
        ->putIntegerVector("level_1", ratio_vector);
    }

    boost::shared_ptr<hier::PatchHierarchy> patch_hierarchy(
        new hier::PatchHierarchy(
          "PatchHierarchy",
          grid_geometry,
          input_db->getDatabase("PatchHierarchy")));

    Cleverleaf* cleverleaf = new Cleverleaf(
        main_db,
        patch_hierarchy,
        dim,
        grid_geometry);

    boost::shared_ptr<LagrangianEulerianLevelIntegrator>
      lagrangian_eulerian_level_integrator(
          new LagrangianEulerianLevelIntegrator(
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

    bool DEV_use_chop_and_pack = false;
    DEV_use_chop_and_pack = main_db->getBoolWithDefault(
        "DEV_use_chop_and_pack",
        DEV_use_chop_and_pack);

    boost::shared_ptr<mesh::LoadBalanceStrategy> load_balancer;

    if (DEV_use_chop_and_pack) {
      load_balancer.reset(new mesh::ChopAndPackLoadBalancer(
            dim,
            "LoadBalancer",
            input_db->getDatabase("LoadBalancer")));
    } else {
      mesh::TreeLoadBalancer* tree_load_balancer = new mesh::TreeLoadBalancer(
          dim,
          "LoadBalancer",
          input_db->getDatabase("LoadBalancer"));
      tree_load_balancer->setSAMRAI_MPI(SAMRAI::tbox::SAMRAI_MPI::getSAMRAIWorld());

      load_balancer.reset(tree_load_balancer);
    }

    boost::shared_ptr<mesh::GriddingAlgorithm> gridding_algorithm(
        new mesh::GriddingAlgorithm(
          patch_hierarchy,
          "GriddingAlgorithm",
          input_db->getDatabase("GriddingAlgorithm"),
          error_detector,
          box_generator,
          load_balancer));

    /*
     * Overwrite regrid_interval option if necessary.
     */
    if(regrid_interval != -1) {
      input_db->getDatabase("LagrangianEulerianIntegrator")
        ->putInteger("regrid_interval", regrid_interval);
    }

    boost::shared_ptr<LagrangianEulerianIntegrator>
      lagrangian_eulerian_integrator(
          new LagrangianEulerianIntegrator(
            input_db->getDatabase("LagrangianEulerianIntegrator"),
            patch_hierarchy,
            lagrangian_eulerian_level_integrator,
            gridding_algorithm));

    int visit_number_procs_per_file = 1;
    std::string visit_dirname = createUniqueFilename(basename, "visit");

    boost::shared_ptr<appu::VisItDataWriter> visit_data_writer(
        new appu::VisItDataWriter(
          dim,
          "Cleverleaf VisIt Writer",
          visit_dirname,
          visit_number_procs_per_file));

    cleverleaf->registerVisItDataWriter(visit_data_writer);

    double dt_now = lagrangian_eulerian_integrator->initializeHierarchy();

    if(vis_dump_interval > 0) {
      visit_data_writer->writePlotData(patch_hierarchy, 0, 0.0);
    }

    tbox::plog << "\nCheck input data and variables before simulation:" << endl;
    tbox::plog << "Input database..." << endl;
    input_db->printClassData(tbox::plog);
    tbox::plog << "\nVariable database..." << endl;
    hier::VariableDatabase::getDatabase()->printClassData(tbox::plog);

    double loop_time = lagrangian_eulerian_integrator->getIntegratorTime();
    double loop_time_end = lagrangian_eulerian_integrator->getEndTime();

    while ((loop_time < loop_time_end) &&
        lagrangian_eulerian_integrator->stepsRemaining()) {

      int iteration_num = lagrangian_eulerian_integrator->getIntegratorStep()+1;

      tbox::pout << "++++++++++++++++++++++++++++++++++++++++++++" << endl;
      tbox::pout << "At begining of timestep # " << iteration_num-1 << endl;
      tbox::pout << "Simulation time is " << loop_time << endl;
      tbox::pout << "Current dt is " << dt_now << endl;

      double dt_new = lagrangian_eulerian_integrator->advanceHierarchy(dt_now);

      loop_time += dt_now;
      dt_now = dt_new;

      if ((field_summary_interval > 0)
          && (iteration_num % field_summary_interval) == 0) {
        lagrangian_eulerian_integrator->printFieldSummary();
      }

      tbox::pout << "At end of timestep # " << iteration_num - 1 << endl;
      tbox::pout << "Simulation time is " << loop_time << endl;
      tbox::pout << "++++++++++++++++++++++++++++++++++++++++++++" << endl;

      if ((vis_dump_interval > 0) && (iteration_num % vis_dump_interval) == 0) {
        visit_data_writer->writePlotData(
            patch_hierarchy,
            iteration_num,
            loop_time);
      }
    }

    if ((vis_dump_interval > 0) && 
        (lagrangian_eulerian_integrator->getIntegratorStep()
         % vis_dump_interval) != 0)
    {
      visit_data_writer->writePlotData(
          patch_hierarchy,
          lagrangian_eulerian_integrator->getIntegratorStep(),
          loop_time);
    }

    if ((field_summary_interval > 0)
        && (lagrangian_eulerian_integrator->getIntegratorStep()
          % field_summary_interval) != 0) 
    {
      lagrangian_eulerian_integrator->printFieldSummary();
    }

    if (main_db->getBoolWithDefault("DEV_check_result", false)) {
      double final_ke = lagrangian_eulerian_integrator->printFieldSummary();
      const double expected_ke = 2.07982664378;
      double ke_difference = fabs((100.0*(final_ke/expected_ke))-100.0);

      tbox::pout << "++++++++++++++++++++++++++++++++++++++++++++" << endl;
      tbox::pout << "Testing kinetic energy..." << std::endl;
      tbox::pout << std::setprecision(12)
        << "    expected: " << expected_ke << std::endl
        << "    actual:   " << final_ke << std::endl;

      if(ke_difference < 0.001) {
        tbox::pout << "Test PASSED." << std::endl;
      } else {
        tbox::pout << "Test FAILED." << std::endl;
      }
      tbox::pout << "++++++++++++++++++++++++++++++++++++++++++++" << endl;
    }

    tbox::TimerManager::getManager()->print(tbox::plog);

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
  size_t start = path.find_last_of("/");
  start += 1;

  size_t end = path.find_last_of(".");

  if(start == string::npos)
    start = 0;

  if(end == string::npos)
    end = path.length();

  size_t count = end - start;

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
