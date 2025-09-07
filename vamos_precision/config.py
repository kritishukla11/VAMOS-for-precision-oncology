import yaml
from pathlib import Path

class Config:
    def __init__(self, path: str | Path | None = None):
        self.data = {}
        if path:
            self.load(path)

    def load(self, path: str | Path):
        path = Path(path)
        with open(path, "r") as f:
            self.data = yaml.safe_load(f) or {}

    def get(self, key, default=None):
        return self.data.get(key, default)

    def __getitem__(self, key):
        return self.data[key]

    def __repr__(self):
        return f"<Config {self.data}>"
