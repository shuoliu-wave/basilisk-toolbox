/* 
 * This script converts Basilisk dump files into .vtu format
 * and writes a single .pvd manifest for visualization in ParaView.
 * 
 * Usage:
 *   qcc -O2 -Wall dump2vtu.c -o dump2vtu -lm
 *   ./dump2vtu
 */
#include "grid/quadtree.h"
#include "navier-stokes/centered.h"
#include "two-phase.h"
#include "output_vtu.h"
#include "output_pvd.h"

scalar step[];  //step
scalar dissrate[];  //dissipation rate in water

char filename_dump[80], filename_vtu[80];

int main() {
  run();
}

event logfile (i = 0) {
  system("mkdir -p vtu");
  int n = 0;
  double timebgn = 0.0;
  double timestp = 0.1;
  double timeend = 2.0;

  for(double timeload = timebgn; timeload <= timeend + 1e-9; timeload += timestp) {
    sprintf(filename_dump,"../dumps/dump%g", timeload);
    if (!restore (file = filename_dump)) {
      fprintf (stderr,
               "[SKIP] Failed to restore '%s' (file missing/corrupt/incompatible)."
               "Skip processing for t=%g.\n", filename_dump, timeload);
      continue;
    }

    scalar l[];
    foreach()
      l[] = level;

    scalar omega[];
    vorticity (u, omega);

    sprintf (filename_vtu, "vtu/%.1f", timeload);
    output_vtu ((scalar *) {f,l,omega,dissrate,step}, (vector *) {u,g}, filename_vtu);

    bool firstTimeWritten = false;
    char pvd_name[80];
    char vtu_file[80];
    sprintf (vtu_file, "vtu/%.1f.vtu", timeload);
    sprintf(pvd_name,"output.pvd");
    FILE *fp = fopen(pvd_name, "r+");
    if( (n == 0) ||  (fp == NULL) ) {
      fp = fopen(pvd_name,"w");
      firstTimeWritten = true;
    }
    output_pvd(vtu_file, timeload, fp, firstTimeWritten);
    fclose(fp);

    fprintf (stderr, "restore '%s', output '%s.vtu', success!\n", filename_dump, filename_vtu);
    n++;
  }
}
