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
#ifndef CLEVERLEAF_MACROS_H_
#define CLEVERLEAF_MACROS_H_

#include <memory>

#ifdef DEBUG
#define SHARED_PTR_CAST(TYPE,VAR) std::dynamic_pointer_cast<TYPE >(VAR)
#else
#define SHARED_PTR_CAST(TYPE,VAR) std::static_pointer_cast<TYPE >(VAR)
#endif // DEBUG

#ifdef DEBUG
#define PTR_CAST(TYPE,VAR) dynamic_cast<TYPE>(VAR)
#else
#define PTR_CAST(TYPE,VAR) static_cast<TYPE>(VAR)
#endif // DEBUG

#endif // CLEVERLEAF_MACROS_H_
