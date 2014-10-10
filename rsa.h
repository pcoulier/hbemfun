/* rsa.h
 *
 * Mattias Schevenels
 * May 2009
 */

#ifndef _RSA_H_
#define _RSA_H_

#include <string>

void rsadec(const std::string& in, const std::string& exp, const std::string& mod, std::string& out);
/* RSA decoder.
 */

void rsaenc(const std::string& in, const std::string& exp, const std::string& mod, std::string& out);
/* RSA encoder.
 */

#endif
