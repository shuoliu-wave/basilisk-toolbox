/* 
 * Export water-phase velocity fields from Basilisk dumps into per-time .dat files.
 * 
 * Usage:
 *   qcc -O2 -Wall output_velocity.c -o output_velocity -lm
 *   ./output_velocity
 *
 * Notes:
 *   - Reads files named like dump0 (no trailing ".0"). If your files are dump0.0, adjust the pattern accordingly.
 *   - Writes per-rank files named velocity-<pid>, then concatenates them into V-<t>.dat.
 *   - Only cells with f >= 0.5 (water) and step == 0 (fluid domain) are written.
 */
#include "grid/quadtree.h"
#include "navier-stokes/centered.h"
#include "two-phase.h"

scalar step[];  // solid mask if present (0 in fluid)
char filename_dump[80], names[80], command[80];

int main() {
  run();
}

event export_velocity(i = 0) {
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

    sprintf(names, "velocity-%d", pid()); 
    FILE * fp = fopen (names, "w");

    scalar omega[];
    vorticity (u, omega);

    foreach_leaf() 
      if (f[] >= 0.5 && step[] == 0.0) {
        fprintf (fp, "%g %g %g %g %g\n", x, y, u.x[], u.y[], omega[]); 
      }
    fclose(fp);

    // Concatenate per-rank files into a single output for this time.
    sprintf(command, "LC_V=C cat velocity* > V-%g.dat", timeload);
    system(command);
    fprintf (stderr, "restore '%s', output 'V-%g.dat', success!\n", filename_dump, timeload);
  }
}
