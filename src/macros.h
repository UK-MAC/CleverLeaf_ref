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
