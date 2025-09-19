/* 
 * Export air-water interface coordinates from Basilisk dumps into per-time .dat files.
 * 
 * Usage:
 *   qcc -O2 -Wall output_interface.c -o output_interface -lm
 *   ./output_interface
 *
 * Notes:
 *   - Read files without extra 0, i.e., dump0 not dump0.0, if dump0.0 is the case, revise accordingly.
 *   - Writes per-rank files named interface-<pid>, then concatenates them into ALL-<t>.dat.
 */
#include "grid/quadtree.h"
#include "navier-stokes/centered.h"
#include "two-phase.h"

scalar step[];  // solid mask if present (0 in fluid)
char filename_dump[80], names[80], command[80];

int main() {
  run();
}

event export_interface(i = 0) {
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
 
    sprintf(names, "interface-%d", pid()); 
    FILE * fp = fopen (names, "w"); 

    output_facets (f,fp);
    fclose(fp);

    // Concatenate per-rank files into a single output for this time.
    sprintf(command, "LC_ALL=C cat interface* > ALL-%g.dat", timeload);
    system(command);
    fprintf (stderr, "restore '%s', output 'ALL-%g.dat', success!\n", filename_dump, timeload);
  }
}
