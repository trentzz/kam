"""Run Steps 3 and 4 only (rebuild sample dirs + training data)."""

import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent))

from run_ml_pipeline import step3_samples, step4_training_data

if __name__ == "__main__":
    sample_failures = step3_samples()
    training_failure = step4_training_data()

    print(f"\n=== Steps 3-4 complete ===", flush=True)
    print(f"Sample failures: {sample_failures}", flush=True)
    print(f"Training data failure: {training_failure}", flush=True)
