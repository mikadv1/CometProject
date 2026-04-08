#include "EphemerisWrapper.h"
#include "ephaccess.h"

static EphAccess* g_eph = nullptr;

bool initEphemeris(const char* filePath) {
    g_eph = ephCreate();
    if (!g_eph) return false;

    int res = ephLoadFile(g_eph, filePath);
    if (res != EPH_OK) {
        return false;
    }

    ephSetDistanceUnits(g_eph, EPH_AU);
    ephSetTimeUnits(g_eph, EPH_DAY);

    return true;
}

StateVector getBodyState(int bodyCode, double jd) {
    double pos[3], vel[3];
    double date0 = (double)(long long)jd;
    double date1 = jd - date0;

    ephCalculateRectangular(g_eph, bodyCode, EPH_SSB, date0, date1, pos, vel);

    return { { pos[0], pos[1], pos[2] }, { vel[0], vel[1], vel[2] } };
}

double TTmTDB(double jd_tdb) {
    double diff = 0.0;
    double date0 = (double)(long long)jd_tdb;
    double date1 = jd_tdb - date0;

    int res = ephCalculateTimeDiff(g_eph, EPH_TT_MINUS_TDB, date0, date1, &diff);

    if (res != EPH_OK) {
        return 0.0;
    }

    return diff;
}

void cleanupEphemeris() {
    if (g_eph) {
        ephDestroy(g_eph);
        g_eph = nullptr;
    }
}