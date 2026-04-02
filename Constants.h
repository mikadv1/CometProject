#pragma once

// Параметры Солнца
const double J2_SUN = 2.21e-7;
const double R_SUN = 695700.0 / 149597870.7;    // 0.00465047 а.е.
const double GM_SUN = 295912208288.32e-15;      // 2.9591220828832e-4
const double SUN_POLE_X = 0.12235349;
const double SUN_POLE_Y = -0.42307208;
const double SUN_POLE_Z = 0.8977971;

// Физические константы (EPM2021, ИПА РАН)
// Единицы: а.е.^3/день^2

const double GM_MERCURY = 49124.8045e-15;       // 4.91248045e-11
const double GM_VENUS = 724345.2333e-15;        // 7.243452333e-10
const double GM_EARTH = 888769.2464e-15;        // 8.887692464e-10
const double GM_MOON = 10931.894556e-15;        // 1.0931894556e-11
const double GM_MARS = 95495.4870e-15;          // 9.54954870e-11
const double GM_JUPITER = 282534582.5972e-15;   // 2.825345825972e-7
const double GM_SATURN = 84597060.7325e-15;     // 8.45970607325e-8
const double GM_URANUS = 12920265.7963e-15;     // 1.29202657963e-8
const double GM_NEPTUNE = 15243573.4789e-15;    // 1.52435734789e-8
const double GM_PLUTO = 2175.0991e-15;          // 2.1750991e-12
const double GM_CERES = 139.64409e-15;
const double GM_PALLAS = 30.8296e-15;
const double GM_VESTA = 38.5475e-15;
const double GM_ERIS = 2450.4387e-15;
const double GM_HAUMEA = 596.1653e-15;
const double GM_MAKEMAKE = 461.3361e-15;

struct ObjectInfo {
    int code;       
    double gm;      
    const char* name;
};

// Массив всех планет для итерации
const ObjectInfo OBJECTS[] = {
    { 1,       GM_MERCURY,  "Mercury"  },
    { 2,       GM_VENUS,    "Venus"    },
    { 399,     GM_EARTH,    "Earth"    },
    { 301,     GM_MOON,     "Moon"     },
    { 4,       GM_MARS,     "Mars"     },
    { 5,       GM_JUPITER,  "Jupiter"  },
    { 6,       GM_SATURN,   "Saturn"   },
    { 7,       GM_URANUS,   "Uranus"   },
    { 8,       GM_NEPTUNE,  "Neptune"  },
    { 9,       GM_PLUTO,    "Pluto"    },
    { 2000001, GM_CERES,    "Ceres"    },
    { 2000002, GM_PALLAS,   "Pallas"   },
    { 2000004, GM_VESTA,    "Vesta"    },
    { 2136199, GM_ERIS,     "Eris"     },
    { 2136108, GM_HAUMEA,   "Haumea"   },
    { 2136472, GM_MAKEMAKE, "Makemake" }
};

const int NUM_OBJECTS = sizeof(OBJECTS) / sizeof(OBJECTS[0]);