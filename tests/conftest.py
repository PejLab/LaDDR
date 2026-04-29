import sys
from pathlib import Path
from types import SimpleNamespace

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

sys.modules.setdefault("statsmodels", SimpleNamespace(api=SimpleNamespace()))
sys.modules.setdefault("statsmodels.api", SimpleNamespace())
