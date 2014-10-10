/* isinf.h
 *
 * Mattias Schevenels
 * July 2008
 */

#ifndef _ISINF_H_
#define _ISINF_H_

#ifndef __GNUC__
#include <limits>
inline bool isinf(const double& x) 
{ 
	return (x>=0 ? x : -x) == std::numeric_limits<double>::infinity(); 
}
#endif

#endif
