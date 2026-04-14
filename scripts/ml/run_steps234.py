"""Run Steps 2, 3, and 4 of the ML pipeline (skip Step 1 — varforge done).

Use this to resume after varforge has completed.
"""

import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent))

from run_ml_pipeline import step2_kam, step3_samples, step4_training_data

if __name__ == "__main__":
    kam_failures = step2_kam(workers=4)
    sample_failures = step3_samples()
    training_failure = step4_training_data()

    print(f"\n=== Steps 2-4 complete ===", flush=True)
    print(f"Kam failures: {kam_failures}", flush=True)
    print(f"Sample failures: {sample_failures}", flush=True)
    print(f"Training data failure: {training_failure}", flush=True)
