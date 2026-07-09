Versioned Laser Simulation Runs
================================

1. Prepare a run folder:
   ./prepare_laser_run.sh RUN_LABEL [--events N] [--queue M] [--include path/to/extra]

2. Inspect or edit the copied files inside runs/RUN_LABEL/ as needed.

3. Submit jobs from the new directory:
   cd runs/RUN_LABEL
   ./submit.sh

4. Outputs and logs stay under runs/RUN_LABEL/ (logs live in log/ and ROOT files in the run root).

Each run folder keeps an immutable snapshot, so you can start another batch elsewhere without touching the files used by jobs that are already in flight.
