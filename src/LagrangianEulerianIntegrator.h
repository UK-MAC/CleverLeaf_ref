/*
 * Copyright 2013 David Beckingsale.
 * 
 * This file is part of CleverLeaf.
 * 
 * CleverLeaf is free software: you can redistribute it and/or modify it under 
 * the terms of the GNU General Public License as published by the 
 * Free Software Foundation, either version 3 of the License, or (at your option) 
 * any later version.
 * 
 * CleverLeaf is distributed in the hope that it will be useful, but 
 * WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or 
 * FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more 
 * details.
 * 
 * You should have received a copy of the GNU General Public License along with 
 * CleverLeaf. If not, see http://www.gnu.org/licenses/.
 */ 

#ifndef included_LagrangianEulerianIntegrator
#define included_LagrangianEulerianIntegrator

#include "LagrangianEulerianLevelIntegrator.h"

#include "SAMRAI/mesh/GriddingAlgorithmStrategy.h"
#include "SAMRAI/hier/PatchHierarchy.h"
#include "SAMRAI/tbox/Database.h"
#include "SAMRAI/tbox/Timer.h"
#include "SAMRAI/tbox/StartupShutdownManager.h"

#include "boost/shared_ptr.hpp"

#include <vector>

using namespace SAMRAI;

class LagrangianEulerianIntegrator {
    public:
        LagrangianEulerianIntegrator(
                const boost::shared_ptr<tbox::Database>& input_db,
                const boost::shared_ptr<hier::PatchHierarchy>& hierarchy,
                const boost::shared_ptr<LagrangianEulerianLevelIntegrator>& level_integrator,
                const boost::shared_ptr<mesh::GriddingAlgorithmStrategy>& gridding_algorithm);

        double initializeHierarchy();

        double advanceHierarchy(const double dt);

        double getIntegratorTime() const
        {
            return d_integrator_time;
        }

        double getStartTime() const
        {
            return d_start_time;
        }

        double getEndTime() const
        {
            return d_end_time;
        }

        int getIntegratorStep() const
        {
            return d_integrator_step;
        }

        bool stepsRemaining() const;

        const boost::shared_ptr<hier::PatchHierarchy> getPatchHierarchy() const
        {
            return d_patch_hierarchy;
        }

        boost::shared_ptr<LagrangianEulerianLevelIntegrator> getLevelIntegrator() const
        {
            return d_level_integrator;
        }

        boost::shared_ptr<mesh::GriddingAlgorithmStrategy> getGriddingAlgorithm() const
        {
            return d_gridding_algorithm;
        }
    private:
        void initializeLevelData(const int level_number);

        void getMinHeirarchyDt(const bool initial_time);

        double d_dt;
        double d_grow_dt;
        double d_max_dt;
        bool d_fix_dt;

        double d_start_time;
        double d_end_time;

        double d_integrator_time;

        int d_integrator_step;
        int d_end_step;

        int d_regrid_interval;
        std::vector<int> d_tag_buffer;

        boost::shared_ptr<hier::PatchHierarchy> d_patch_hierarchy;
        boost::shared_ptr<LagrangianEulerianLevelIntegrator> d_level_integrator;
        boost::shared_ptr<mesh::GriddingAlgorithmStrategy> d_gridding_algorithm;

        static boost::shared_ptr<tbox::Timer> t_initialize_hierarchy;
        static boost::shared_ptr<tbox::Timer> t_advance_hierarchy;

        static void initializeCallback();
        static void finalizeCallback();

        static tbox::StartupShutdownManager::Handler s_initialize_handler;
};

#endif
