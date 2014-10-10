/* edtlicense.h
 *
 * Mattias Schevenels
 * July 2009
 */

#ifndef _BEMFUNLICENSE_H_
#define _BEMFUNLICENSE_H_

int bemfunlicense(const char* const action);
/* Check BEMFUN license and show license information.
 * If no license file is found, BEMFUN will run in trial mode.
 *
 * mode   "VerifyAlways"   Force license verification.
 *        "VerifyOnce"     Verify license if not verified yet.
 *        "Reset"          Reset license verification status.
 *
 * The return value is the BEMFUN license status:
 *
 *   -1   Unverified
 *    0   OK
 *    1   License file not found
 *    2   Expired
 *    3   Incorrect MAC address
 *    4   Invalid license file
 */

#endif
