/*  This script converts Basilisk dump files into .htg format for visualization in ParaView.
    Note that in the exported files, the x–y axes are swapped in 2D and the x–z axes are swapped in 3D, including for vector components such as u and g.  */
#include "grid/quadtree.h"
#include "navier-stokes/centered.h"
#include "two-phase.h"
#include "output_htg.h"

vector U[], G[];
scalar step[];  //step
scalar dissrate[];  //dissipation rate in water

char filename_dump[80], filename_htg[80], path[]="htg";

int main() {
  run();
}

event logfile (i = 0) {
  system("mkdir -p htg");
  int i = 0;
  double timebgn = 0.0;
  double timestp = 0.05;
  double timeend = 2.0;

  for(double timeload = timebgn; timeload <= timeend + 1e-9; timeload += timestp) {
    sprintf(filename_dump,"../dumps/dump%g", timeload);
    if (!restore (file = filename_dump)) {
      fprintf (stderr,
               "[SKIP] Failed to restore '%s' (file missing/corrupt/incompatible)."
               "Skip processing for t=%g.\n", filename_dump, timeload);
      continue;
    }

    foreach() {
      foreach_dimension() {
        U.x[] = u.x[];
        G.x[] = g.x[];
      }
    }

    scalar l[];
    foreach()
      l[] = level;

    scalar omega[];
    vorticity (U, omega);

    sprintf (filename_htg, "%.1f", timeload);
    output_htg ((scalar *) {f,l,omega,dissrate,step}, (vector *) {U,G}, path, filename_htg, i, timeload);
    fprintf (stderr, "restore '%s', output '%s/%s.htg', success!\n", filename_dump, path, filename_htg);
    i++;
  }
}
