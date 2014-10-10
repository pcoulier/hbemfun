/* ripemd128.h
 *
 * Mattias Schevenels
 * May 2009
 */

#ifndef _RIPEMD128_H_
#define _RIPEMD128_H_

#include <string>

void rmdstring(const std::string& s, std::string& hash);
/* Compute RMD-128 hash for string s.
 */

void rmdfile(const std::string& file, std::string& hash);
/* Compute RMD-128 hash for specified file.
 */

#endif
