#include "FileIO.h"
#include <cstdio>
#include <cmath>

bool loadTrajectoryFromCSV(const char* filename, Trajectory& traj, double t0, double tend, double dt) {
    FILE* file = fopen(filename, "r");
    if (!file) {
        printf("[ERROR] Cannot open %s\n", filename);
        return false;
    }

    int nPoints = static_cast<int>((tend - t0) / dt) + 1;
    traj.allocate(nPoints, t0, dt);

    char line[200];
    fgets(line, sizeof(line), file); 

    int count = 0;
    while (fgets(line, sizeof(line), file)) {
        double jd, x, y, z, vx, vy, vz;
        if (sscanf(line, "%lf,%*[^,],%lf,%lf,%lf,%lf,%lf,%lf",
            &jd, &x, &y, &z, &vx, &vy, &vz) != 7) {
            printf("[WARN] Parse error at line %d\n", count + 2);
            printf("%s", line);
            continue;
        }
        traj.t[count] = jd;
        traj.state[count].r = Vector3(x, y, z);
        traj.state[count].v = Vector3(vx, vy, vz);
        count++;
    }
    fclose(file);

    printf("Expected %d points, loaded %d\n", nPoints, count);

    return count > 0;
}

bool writeResiduals(const std::vector<Observation>& obs, const char* outFile, const Trajectory& traj) {
    FILE* fout = fopen(outFile, "w");
    if (!fout) {
        printf("[ERROR] Cannot open output files\n");
        return false;
    } 

    fprintf(fout, "JD,dRA,dDec,sigmaRA,sigmaDec,Code\n");
    int count = 0;
    for (const Observation& ob : obs) {
        ReductionResult res;
        reduceObservation(ob, traj, res);
        //printf("%s %.6f\n", ob.code, res.jd_tdb);
        fprintf(fout, "%.8f,%.4f,%.4f,%.4f,%.4f,%s\n", ob.jd_utc, res.dRA, res.dDec, ob.sigma_ra, ob.sigma_dec, ob.code);
        count++;
    }

    fclose(fout);

    printf("Processed %d observations\n", count);
    return true;
}

bool loadObservations(const char* obsFile, std::vector<Observation>& obs_vector) {
    FILE* fin = fopen(obsFile, "r");
    if (!fin) {
        printf("[ERROR] Cannot open input file\n");
        if (fin) fclose(fin);
        return false;
    }

    char line[256];
    fgets(line, sizeof(line), fin);

    int count = 0;

    while (fgets(line, sizeof(line), fin)) {
        double jd = 0, ra = 0, dec = 0, sigma_ra, sigma_dec, x = 0, y = 0, z = 0;
        char code[4] = {}; char type = '\0';

        if (sscanf(line, "%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%3[^,],%c",
            &jd, &ra, &dec, &sigma_ra, &sigma_dec, &x, &y, &z, code, &type) != 10) {
            printf("Error while reading line %d %s\n",count + 2, line);
            continue;
        }

        Observation observation;
        observation.jd_utc = jd;
        observation.ra = ra;
        observation.dec = dec;
        observation.sigma_ra = sigma_ra;
        observation.sigma_dec = sigma_dec;
        observation.r_station_itrf = Vector3(x, y, z);
        observation.type = type;
        strncpy(observation.code, code, 3);
        observation.code[3] = '\0';  

        obs_vector.push_back(observation);

        count++;
    }

    fclose(fin);

    printf("Load %d observations\n", count);
    return true;
}
