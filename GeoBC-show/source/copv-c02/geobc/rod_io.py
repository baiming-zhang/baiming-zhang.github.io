"""Atomic JSON products with retry for transient Windows reader sharing locks."""
import json
from pathlib import Path
import time


def write_json(path, data):
    path=Path(path)
    path.parent.mkdir(parents=True,exist_ok=True)
    temp=path.with_suffix(path.suffix+'.tmp')
    temp.write_text(json.dumps(data,indent=2,allow_nan=False)+'\n',encoding='utf8')
    for attempt in range(101):
        try:
            temp.replace(path)
            return
        except PermissionError:
            if attempt==100:raise
            time.sleep(.05)
