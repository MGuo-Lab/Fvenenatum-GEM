"""Make the repo-root scripts (flat modules, not a package) importable in tests."""

import os
import sys

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
