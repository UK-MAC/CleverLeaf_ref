/*
 * Copyright 2013 David Beckingsale.
 * 
 * This file is part of CleverLeaf.
 * 
 * CleverLeaf is free software: you can redistribute it and/or modify it under
 * the terms of the GNU General Public License as published by the Free
 * Software Foundation, either version 3 of the License, or (at your option)
 * any later version.
 * 
 * CleverLeaf is distributed in the hope that it will be useful, but WITHOUT
 * ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
 * FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for
 * more details.
 * 
 * You should have received a copy of the GNU General Public License along with
 * CleverLeaf. If not, see http://www.gnu.org/licenses/.
 */ 
#ifndef CLEVER_GEOM_CARTESIANCELLDOUBLEMASSWEIGHTEDAVERAGE_H_
#define CLEVER_GEOM_CARTESIANCELLDOUBLEMASSWEIGHTEDAVERAGE_H_

#include <string>

#include "SAMRAI/hier/CoarsenOperator.h"
#include "SAMRAI/hier/Box.h"
#include "SAMRAI/hier/IntVector.h"
#include "SAMRAI/hier/Patch.h"
#include "SAMRAI/tbox/StartupShutdownManager.h"
#include "SAMRAI/tbox/Timer.h"
#include "SAMRAI/tbox/TimerManager.h"

using namespace SAMRAI;

namespace clever {
namespace geom {

class CartesianCleverCellDoubleMassWeightedAverage:
  public hier::CoarsenOperator
{
  public:
    explicit CartesianCleverCellDoubleMassWeightedAverage(
        const tbox::Dimension& dim);

    virtual ~CartesianCleverCellDoubleMassWeightedAverage();

    bool findCoarsenOperator(const boost::shared_ptr<hier::Variable>& var,
        const std::string& op_name) const;

    int getOperatorPriority() const;

    hier::IntVector getStencilWidth(const tbox::Dimension &dim) const;

    void coarsen(
        hier::Patch& coarse,
        const hier::Patch& fine,
        const int dst_component,
        const int src_component,
        const hier::Box& coarse_box,
        const hier::IntVector& ratio) const;
  private:
    static boost::shared_ptr<tbox::Timer> t_coarsen;

    static void initializeCallback();
    static void finalizeCallback();

    static tbox::StartupShutdownManager::Handler s_initialize_handler;
};

}
}
#endif // CLEVER_GEOM_CARTESIANCELLDOUBLEMASSWEIGHTEDAVERAGE_H_
