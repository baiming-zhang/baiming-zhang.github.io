"""Read-only scheduled monitor; never launches/restarts or alters native data."""
from pathlib import Path
import json,csv,subprocess,datetime
R=Path(__file__).resolve().parents[3];O=R/'results/uniform5_pilot_20261113'
def main():
    if not (O/'status.json').exists():print(json.dumps(dict(state='not_prepared',action='No native job; do not start or restart.')));return
    s=json.loads((O/'status.json').read_text());out=dict(checked_utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),batch_state=s['state'],terminal_cases=[dict(name=c['name'],returncode=c['returncode'],native_gate=c['native_gate']) for c in s['cases']],full_C1_complete=False)
    if (O/'live.json').exists():
        live=json.loads((O/'live.json').read_text());out['live']=live
        if live['state']=='running':
            pid=int(live['native_pid']);cmd='[Console]::OutputEncoding = [System.Text.Encoding]::UTF8; Get-CimInstance Win32_Process -Filter "ProcessId='+str(pid)+'" | Select-Object ProcessId,CreationDate,CommandLine | ConvertTo-Json -Compress';raw=subprocess.check_output(['powershell','-NoProfile','-Command',cmd],text=True,encoding='utf-8');proc=json.loads(raw) if raw.strip() else None;out['observed_process']=proc;out['identity_matches']=bool(proc and live['output_directory'].casefold() in proc['CommandLine'].casefold() and 'tf_local_assembly_bench.exe' in proc['CommandLine'].casefold());root=Path(live['output_directory']);history=root/'history.csv';steps=[]
            if history.exists():
                with history.open(encoding='utf-8') as stream:
                    for row in csv.DictReader(stream):
                        try:
                            if None in row or any(v is None or v=='' for v in row.values()):continue
                            [float(v) for v in row.values()];steps.append(int(row['step']))
                        except ValueError:continue
            out['last_complete_history_step']=max(steps) if steps else 0
            if proc is None:out['action']='Re-read live/status and inspect SAME parent/handle; a transient missing PID is not permission to restart.'
    print(json.dumps(out,ensure_ascii=True,indent=2))
if __name__=='__main__':main()
