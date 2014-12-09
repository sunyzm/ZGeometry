#include "ColorMap.h"

namespace ZGeom {

double ColorMap::parula[] = {

double ColorMap::jet[] = { 0, 0, 0.51563,
    0, 0, 0.53125,
    0, 0, 0.54688,
    0, 0, 0.5625,
    0, 0, 0.57813,
    0, 0, 0.59375,
    0, 0, 0.60938,
    0, 0, 0.625,
    0, 0, 0.64063,
    0, 0, 0.65625,
    0, 0, 0.67188,
    0, 0, 0.6875,
    0, 0, 0.70313,
    0, 0, 0.71875,
    0, 0, 0.73438,
    0, 0, 0.75,
    0, 0, 0.76563,
    0, 0, 0.78125,
    0, 0, 0.79688,
    0, 0, 0.8125,
    0, 0, 0.82813,
    0, 0, 0.84375,
    0, 0, 0.85938,
    0, 0, 0.875,
    0, 0, 0.89063,
    0, 0, 0.90625,
    0, 0, 0.92188,
    0, 0, 0.9375,
    0, 0, 0.95313,
    0, 0, 0.96875,
    0, 0, 0.98438,
    0, 0, 1,
    0, 0.015625, 1,
    0, 0.03125, 1,
    0, 0.046875, 1,
    0, 0.0625, 1,
    0, 0.078125, 1,
    0, 0.09375, 1,
    0, 0.10938, 1,
    0, 0.125, 1,
    0, 0.14063, 1,
    0, 0.15625, 1,
    0, 0.17188, 1,
    0, 0.1875, 1,
    0, 0.20313, 1,
    0, 0.21875, 1,
    0, 0.23438, 1,
    0, 0.25, 1,
    0, 0.26563, 1,
    0, 0.28125, 1,
    0, 0.29688, 1,
    0, 0.3125, 1,
    0, 0.32813, 1,
    0, 0.34375, 1,
    0, 0.35938, 1,
    0, 0.375, 1,
    0, 0.39063, 1,
    0, 0.40625, 1,
    0, 0.42188, 1,
    0, 0.4375, 1,
    0, 0.45313, 1,
    0, 0.46875, 1,
    0, 0.48438, 1,
    0, 0.5, 1,
    0, 0.51563, 1,
    0, 0.53125, 1,
    0, 0.54688, 1,
    0, 0.5625, 1,
    0, 0.57813, 1,
    0, 0.59375, 1,
    0, 0.60938, 1,
    0, 0.625, 1,
    0, 0.64063, 1,
    0, 0.65625, 1,
    0, 0.67188, 1,
    0, 0.6875, 1,
    0, 0.70313, 1,
    0, 0.71875, 1,
    0, 0.73438, 1,
    0, 0.75, 1,
    0, 0.76563, 1,
    0, 0.78125, 1,
    0, 0.79688, 1,
    0, 0.8125, 1,
    0, 0.82813, 1,
    0, 0.84375, 1,
    0, 0.85938, 1,
    0, 0.875, 1,
    0, 0.89063, 1,
    0, 0.90625, 1,
    0, 0.92188, 1,
    0, 0.9375, 1,
    0, 0.95313, 1,
    0, 0.96875, 1,
    0, 0.98438, 1,
    0, 1, 1,
    0.015625, 1, 0.98438,
    0.03125, 1, 0.96875,
    0.046875, 1, 0.95313,
    0.0625, 1, 0.9375,
    0.078125, 1, 0.92188,
    0.09375, 1, 0.90625,
    0.10938, 1, 0.89063,
    0.125, 1, 0.875,
    0.14063, 1, 0.85938,
    0.15625, 1, 0.84375,
    0.17188, 1, 0.82813,
    0.1875, 1, 0.8125,
    0.20313, 1, 0.79688,
    0.21875, 1, 0.78125,
    0.23438, 1, 0.76563,
    0.25, 1, 0.75,
    0.26563, 1, 0.73438,
    0.28125, 1, 0.71875,
    0.29688, 1, 0.70313,
    0.3125, 1, 0.6875,
    0.32813, 1, 0.67188,
    0.34375, 1, 0.65625,
    0.35938, 1, 0.64063,
    0.375, 1, 0.625,
    0.39063, 1, 0.60938,
    0.40625, 1, 0.59375,
    0.42188, 1, 0.57813,
    0.4375, 1, 0.5625,
    0.45313, 1, 0.54688,
    0.46875, 1, 0.53125,
    0.48438, 1, 0.51563,
    0.5, 1, 0.5,
    0.51563, 1, 0.48438,
    0.53125, 1, 0.46875,
    0.54688, 1, 0.45313,
    0.5625, 1, 0.4375,
    0.57813, 1, 0.42188,
    0.59375, 1, 0.40625,
    0.60938, 1, 0.39063,
    0.625, 1, 0.375,
    0.64063, 1, 0.35938,
    0.65625, 1, 0.34375,
    0.67188, 1, 0.32813,
    0.6875, 1, 0.3125,
    0.70313, 1, 0.29688,
    0.71875, 1, 0.28125,
    0.73438, 1, 0.26563,
    0.75, 1, 0.25,
    0.76563, 1, 0.23438,
    0.78125, 1, 0.21875,
    0.79688, 1, 0.20313,
    0.8125, 1, 0.1875,
    0.82813, 1, 0.17188,
    0.84375, 1, 0.15625,
    0.85938, 1, 0.14063,
    0.875, 1, 0.125,
    0.89063, 1, 0.10938,
    0.90625, 1, 0.09375,
    0.92188, 1, 0.078125,
    0.9375, 1, 0.0625,
    0.95313, 1, 0.046875,
    0.96875, 1, 0.03125,
    0.98438, 1, 0.015625,
    1, 1, 0,
    1, 0.98438, 0,
    1, 0.96875, 0,
    1, 0.95313, 0,
    1, 0.9375, 0,
    1, 0.92188, 0,
    1, 0.90625, 0,
    1, 0.89063, 0,
    1, 0.875, 0,
    1, 0.85938, 0,
    1, 0.84375, 0,
    1, 0.82813, 0,
    1, 0.8125, 0,
    1, 0.79688, 0,
    1, 0.78125, 0,
    1, 0.76563, 0,
    1, 0.75, 0,
    1, 0.73438, 0,
    1, 0.71875, 0,
    1, 0.70313, 0,
    1, 0.6875, 0,
    1, 0.67188, 0,
    1, 0.65625, 0,
    1, 0.64063, 0,
    1, 0.625, 0,
    1, 0.60938, 0,
    1, 0.59375, 0,
    1, 0.57813, 0,
    1, 0.5625, 0,
    1, 0.54688, 0,
    1, 0.53125, 0,
    1, 0.51563, 0,
    1, 0.5, 0,
    1, 0.48438, 0,
    1, 0.46875, 0,
    1, 0.45313, 0,
    1, 0.4375, 0,
    1, 0.42188, 0,
    1, 0.40625, 0,
    1, 0.39063, 0,
    1, 0.375, 0,
    1, 0.35938, 0,
    1, 0.34375, 0,
    1, 0.32813, 0,
    1, 0.3125, 0,
    1, 0.29688, 0,
    1, 0.28125, 0,
    1, 0.26563, 0,
    1, 0.25, 0,
    1, 0.23438, 0,
    1, 0.21875, 0,
    1, 0.20313, 0,
    1, 0.1875, 0,
    1, 0.17188, 0,
    1, 0.15625, 0,
    1, 0.14063, 0,
    1, 0.125, 0,
    1, 0.10938, 0,
    1, 0.09375, 0,
    1, 0.078125, 0,
    1, 0.0625, 0,
    1, 0.046875, 0,
    1, 0.03125, 0,
    1, 0.015625, 0,
    1, 0, 0,
    0.98438, 0, 0,
    0.96875, 0, 0,
    0.95313, 0, 0,
    0.9375, 0, 0,
    0.92188, 0, 0,
    0.90625, 0, 0,
    0.89063, 0, 0,
    0.875, 0, 0,
    0.85938, 0, 0,
    0.84375, 0, 0,
    0.82813, 0, 0,
    0.8125, 0, 0,
    0.79688, 0, 0,
    0.78125, 0, 0,
    0.76563, 0, 0,
    0.75, 0, 0,
    0.73438, 0, 0,
    0.71875, 0, 0,
    0.70313, 0, 0,
    0.6875, 0, 0,
    0.67188, 0, 0,
    0.65625, 0, 0,
    0.64063, 0, 0,
    0.625, 0, 0,
    0.60938, 0, 0,
    0.59375, 0, 0,
    0.57813, 0, 0,
    0.5625, 0, 0,
    0.54688, 0, 0,
    0.53125, 0, 0,
    0.51563, 0, 0,
    0.5, 0, 0 };

} // end of namespace