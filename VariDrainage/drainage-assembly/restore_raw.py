"""Restore verified raw CSV bytes from downloaded ZIP parts; refuses overwrite.
Usage: python restore_raw.py REPORT_DIRECTORY NEW_OUTPUT_DIRECTORY
REPORT_DIRECTORY contains download_manifest.json and downloads/*.zip.
"""
from pathlib import Path,PurePosixPath
import sys,json,zipfile,hashlib
def main():
    source=Path(sys.argv[1]);out=Path(sys.argv[2]);manifest=json.loads((source/'download_manifest.json').read_text(encoding='utf-8'))
    for part in manifest['parts']:
        p=source/part['path'];assert hashlib.sha256(p.read_bytes()).hexdigest()==part['sha256'],p
    for f in manifest['files']:
        relative=PurePosixPath(f['path']);assert not relative.is_absolute() and '..' not in relative.parts
        target=out/Path(*relative.parts);assert not target.exists(),target;data=b''
        for c in f['chunks']:
            with zipfile.ZipFile(source/'downloads'/c['archive']) as z:chunk=z.read(c['entry'])
            assert hashlib.sha256(chunk).hexdigest()==c['sha256'];assert len(data)==c['offset'];data+=chunk
        assert len(data)==f['bytes'] and hashlib.sha256(data).hexdigest()==f['sha256'];target.parent.mkdir(parents=True,exist_ok=True);target.write_bytes(data)
    print('Verified restored files:',len(manifest['files']))
if __name__=='__main__':main()
