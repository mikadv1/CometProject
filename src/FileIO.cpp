#include "FileIO.h"
#include "Constants.h"
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

bool processObservations(const char* obsFile, const char* outFile, const Trajectory& traj) {
    FILE* fin = fopen(obsFile, "r");
    FILE* fout = fopen(outFile, "w");
    if (!fin || !fout) {
        printf("[ERROR] Cannot open input/output files\n");
        if (fin) fclose(fin);
        if (fout) fclose(fout);
        return false;
    }

    fprintf(fout, "JD,dRA,dDec,Code\n");

    char line[256];
    fgets(line, sizeof(line), fin);  

    int count = 0;

    while (fgets(line, sizeof(line), fin)) {
        double jd, ra, dec, x, y, z;
        int code; char mode;

        if (sscanf(line, "%lf,%lf,%lf,%lf,%lf,%lf,%d,%c",
            &jd, &ra, &dec, &x, &y, &z, &code, &mode) != 8) {
            continue;
        }
        if (mode == 'S') continue;
        if (jd < traj.t0 || jd > traj.tend) continue;

        Vector3 r_station(x, y, z);
        ReductionResult res;
        Observation obs = { jd, ra, dec, r_station };

        reduceObservation(obs, traj, res);

        fprintf(fout, "%.8f,%.4f,%.4f,%d\n", jd, res.dRA, res.dDec, code);

        count++;
    }

    fclose(fin);
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
        double jd = 0, ra = 0, dec = 0, x = 0, y = 0, z = 0;
        int code = 0; char type = '\0';

        if (sscanf(line, "%lf,%lf,%lf,%lf,%lf,%lf,%d,%c",
            &jd, &ra, &dec, &x, &y, &z, &code, &type) != 8) {
            printf("Error while reading line %d %s\n",count + 2, line);
            printf("%.6f %.3f %.3f %.3f %.3f %.3f %d %c\n", jd, ra, dec, x, y, z, code, type);
            continue;
        }
        if (type == 'S') continue;

        Vector3 r_station(x, y, z);
        obs_vector.push_back({ jd, ra, dec, r_station });

        count++;
    }

    fclose(fin);

    printf("Load %d observations\n", count);
    return true;
}
