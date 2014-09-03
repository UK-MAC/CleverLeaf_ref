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
#include "LagrangianEulerianPatchStrategy.h"

LagrangianEulerianPatchStrategy::LagrangianEulerianPatchStrategy(
    const tbox::Dimension& dim):
  RefinePatchStrategy(),
  d_dim(dim)
{
}

void LagrangianEulerianPatchStrategy::setCurrentDataContext(
    boost::shared_ptr<hier::VariableContext> context)
{
  d_current_data_context = context;
}

void LagrangianEulerianPatchStrategy::setNewDataContext(
    boost::shared_ptr<hier::VariableContext> context)
{
  d_new_data_context = context;
}

boost::shared_ptr<hier::VariableContext>
LagrangianEulerianPatchStrategy::getCurrentDataContext()
{
  return d_current_data_context;
}

boost::shared_ptr<hier::VariableContext>
LagrangianEulerianPatchStrategy::getNewDataContext()
{
  return d_new_data_context;
}

const tbox::Dimension& LagrangianEulerianPatchStrategy::getDim() const
{
  return d_dim;
}

hier::IntVector LagrangianEulerianPatchStrategy::getRefineOpStencilWidth(
    const tbox::Dimension &dim ) const 
{
  return hier::IntVector(dim,0);
}

void LagrangianEulerianPatchStrategy::preprocessRefine(
    hier::Patch& fine,
    const hier::Patch& coarse,
    const hier::Box& fine_box,
    const hier::IntVector& ratio)
{
  /*
   * Dummy implementation!
   */
}

void LagrangianEulerianPatchStrategy::postprocessRefine(
    hier::Patch& fine,
    const hier::Patch& coarse,
    const hier::Box& fine_box,
    const hier::IntVector& ratio)
{
  /*
   * Dummy implementation!
   */
}

void LagrangianEulerianPatchStrategy::setExchangeFlag(const int exchange)
{
  d_which_exchange = exchange;
}
