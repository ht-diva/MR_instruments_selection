import pandas as pd
from pathlib import Path


def ws_path(file_path):
    Path(config.get("workspace_path")).mkdir(parents=True, exist_ok=True)
    return str(Path(config.get("workspace_path"), file_path))


def dst_path(file_path):
    return str(Path(config.get("dest_path"), file_path))
