# Basilisk dump → ParaView HyperTreeGrid (.htg) exporter

## Purpose
- Batch-read dump%g files from ../dumps/ (time-stamped), and write .htg files into htg/ for ParaView.
- The output preserves the AMR structure via HyperTreeGrid and includes key flow fields.
 
## Exported fields
- Scalars:
  - f         : VOF volume fraction
  - l         : AMR refinement level (alias of `level`)
  - omega     : vorticity magnitude computed from U
  - dissrate  : dissipation rate in water (restored if present, else 0)
  - step      : solid/step mask (restored if present, else 0)
- Vectors:
  - U         : copy of velocity `u`
  - G         : copy of body force/gravity `g`

## Coordinate / component remapping (applied consistently to vectors)
- 2D: swap x ↔ y
- 3D: swap x ↔ z
- Vector components are permuted in the same way as coordinates.
 
## How to use
1) Set the export time window:
   timebgn, timestp, timeend
   e.g. 0 → 2 by 0.1.
2) Adjust the dump path pattern if needed (`filename_dump`).
3) Compile & run (example):
   ```bash
   qcc -O2 -Wall dump2htg.c -o dump2htg -lm
   ./dump2htg
   # The program creates htg/ and writes <time>.htg files there.
4) Open htg/*.htg directly in ParaView.
   `Set view direction to +Z`, then `Rotate 90° clockwise` to align the display with the simulation's physical axes.
 
 ## Notes & assumptions
 - Compile in the same dimension (2D/3D) as the original simulation.
 - `restore()` recovers fields that exist in the dump; missing fields remain at their default (e.g., zero).
 - The actual HTG writing and axis/component permutation are handled by output_htg.h.
 
 ## Logging
 - Each successful export logs to stderr:
   restored dump<time>, wrote <time>.htg
 
 ## Troubleshooting
 - If dumps are not found or permissions fail, check relative paths and file access rights.
 - To export additional fields, append them to the argument lists of output_htg().




