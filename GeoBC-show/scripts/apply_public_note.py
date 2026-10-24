"""Apply only two bounded plain-text task notes from a public GitHub issue.

Issue content never becomes executable code, HTML, a path, or a shell argument.
"""
import json,os,re,sys,datetime
from pathlib import Path
IDS={f'R{i:02}' for i in range(1,9)}|{f'C{i:02}' for i in range(1,7)}|{'M01','M02','M03','E01','F01'}
def parse(issue):
    body=issue.get('body') or ''
    if len(body)>15000:raise ValueError('Issue body too long')
    m=re.fullmatch(r'<!-- geobc-task-note-v1 -->\s*```json\s*\n(.*?)\n```\s*',body,re.S)
    if not m:raise ValueError('Only the structured public note form is accepted')
    data=json.loads(m[1]);assert isinstance(data,dict) and set(data)=={'task','nickname','summary','review'}
    assert isinstance(data['task'],str) and data['task'] in IDS
    assert issue.get('title')=='[GeoBC note] '+data['task']
    for key,limit in [('nickname',40),('summary',2000),('review',4000)]:
        value=data[key];assert isinstance(value,str) and 0<len(value.strip())<=limit
        assert not any(ord(c)<32 and c not in '\n\t\r' for c in value)
        data[key]=value.strip()
    return data
def tests():
    def issue(d):return dict(title='[GeoBC note] R01',body='<!-- geobc-task-note-v1 -->\n```json\n'+json.dumps(d)+'\n```')
    valid=dict(task='R01',nickname='Test reviewer',summary='Updated geometry',review='<script>alert(1)</script>')
    assert parse(issue(valid))==valid # retained as text, escaped by textContent at render time
    for bad in [dict(valid,task='../../README'),dict(valid,summary=''),dict(valid,review='a'*4001),dict(valid,arbitrary='field')]:
        try:parse(issue(bad))
        except (AssertionError,ValueError):pass
        else:raise AssertionError('Invalid public input accepted')
    try:parse(dict(title='[GeoBC note] R01',body='$(echo unsafe)'))
    except ValueError:pass
    else:raise AssertionError('Unstructured input accepted')
    print('Public note parser checks passed')
def main():
    if '--self-test' in sys.argv:tests();return
    event=json.loads(Path(os.environ['GITHUB_EVENT_PATH']).read_text(encoding='utf-8'))
    if 'issue' not in event:
        inputs=event['inputs'];data={k:inputs[k] for k in ['task','nickname','summary','review']}
        event['issue']=dict(title='[GeoBC note] '+data['task'],body='<!-- geobc-task-note-v1 -->\n```json\n'+json.dumps(data)+'\n```',user=dict(login=os.environ['GITHUB_ACTOR']),html_url='https://github.com/'+os.environ['GITHUB_REPOSITORY']+'/actions/runs/'+os.environ['GITHUB_RUN_ID'],updated_at=datetime.datetime.now(datetime.timezone.utc).isoformat())
    data=parse(event['issue'])
    path=Path('GeoBC-show/data/task-notes.json');notes=json.loads(path.read_text(encoding='utf-8'));assert data['task'] in notes['notes']
    notes['notes'][data['task']]=dict(summary=data['summary'],review=data['review'],nickname=data['nickname'],author=event['issue']['user']['login'],issue=event['issue']['html_url'],updated=event['issue']['updated_at'])
    path.write_text(json.dumps(notes,ensure_ascii=False,indent=2)+'\n',encoding='utf-8')
if __name__=='__main__':main()
