# Basilisk dump → ParaView VTK Unstructured Grid Files (.vtu) exporter

## Purpose
- Batch-read dump%g files from ../dumps/ (time-stamped), and write .vtu files into vtu/ for ParaView.
 
## Exported fields
- Scalars:
  - f         : VOF volume fraction
  - l         : AMR refinement level (alias of `level`)
  - omega     : vorticity magnitude computed from u
  - dissrate  : dissipation rate in water (restored if present, else 0)
  - step      : solid/step mask (restored if present, else 0)
- Vectors:
  - u         : velocity
  - g         : body force/gravity
 
## How to use
1) Set the export time window:
   timebgn, timestp, timeend
   e.g. 0 → 2 by 0.1.
2) Adjust the dump path pattern if needed (`filename_dump`).
3) Compile & run (example):
   ```bash
   qcc -O2 -Wall dump2vtu.c -o dump2vtu -lm
   ./dump2vtu
   ```
   The program creates vtu/ and writes <time>.vtu files there.
4) Open vtu/*.vtu directly in ParaView.
 
 ## Notes & assumptions
 - Compile in the same dimension (2D/3D) as the original simulation.
 - `restore()` recovers fields that exist in the dump; missing fields remain at their default (e.g., zero).
 
 ## Logging
 - Each successful export logs to stderr:
   restored dump<time>, wrote <time>.htg
 
 ## Troubleshooting
 - If dumps are not found or permissions fail, check relative paths and file access rights.
 - To export additional fields, append them to the argument lists of output_vtu().