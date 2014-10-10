#ifndef _GREENROTATE3D_
#define _GREENROTATE3D_
double fminstep(double(*funToMinimize)(std::valarray<double>&, const void** const),
                std::valarray<double>& x, const std::valarray<double>& xRes,
                const std::valarray<double>& xTol, const int& nMaxIter,
                const void** const varargin);
#endif
