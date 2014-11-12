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

#include "boost/version.hpp"

#if BOOST_VERSION > 105200

#ifdef DEBUG
#define SHARED_PTR_CAST(TYPE,VAR) boost::dynamic_pointer_cast<TYPE >(VAR)
#else
#define SHARED_PTR_CAST(TYPE,VAR) boost::static_pointer_cast<TYPE >(VAR)
#endif // DEBUG

#else

#ifdef DEBUG
#define SHARED_PTR_CAST(TYPE,VAR) VAR, boost::detail::dynamic_cast_tag()
#else
#define SHARED_PTR_CAST(TYPE,VAR) VAR, boost::detail::static_cast_tag()
#endif // DEBUG

#endif // BOOST_VERSION > 105200


#ifdef DEBUG
#define PTR_CAST(TYPE,VAR) dynamic_cast<TYPE>(VAR)
#else
#define PTR_CAST(TYPE,VAR) static_cast<TYPE>(VAR)
#endif // DEBUG

#endif // CLEVERLEAF_MACROS_H_
