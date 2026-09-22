"""Pytest configuration: make the ``src`` layout importable without installation."""

# =============================================================================
# IMPORTS
# =============================================================================
import sys
from pathlib import Path

# =============================================================================
# CORE LOGIC
# =============================================================================
SRC_DIR = (Path(__file__).parent / "../src").resolve()
if str(SRC_DIR) not in sys.path:
    sys.path.insert(0, str(SRC_DIR))
