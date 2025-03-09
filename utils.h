//
// Created by Yuzhen Feng on 19/2/2025.
//

#ifndef UTILS_H
#define UTILS_H

#include <cmath>

namespace math2 {
    inline double opt_pow(const double &_base, const double &_exp) {
        if (static_cast<int>(_exp) == _exp) {
            if (_exp == 4.) {
                return _base * _base * _base * _base;
            }

            if (_exp == 3.) {
                return _base * _base * _base;
            }

            if (_exp == 0.) {
                return 1.;
            }

            if (_exp < 0.) {
                if (_base != 0.) {
                    return std::pow(_base, _exp);
                }
                return 0.;
            }

            double tmp = _base;
            unsigned short int i = 0;
            while (i++ < _exp - 1)
                tmp *= _base;
            return tmp;
        }

        return std::pow(_base, _exp);
    }
}

#endif //UTILS_H
