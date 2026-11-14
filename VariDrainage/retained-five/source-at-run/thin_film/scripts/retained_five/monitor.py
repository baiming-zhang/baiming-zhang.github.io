"""Read-only native PID/command/output identity check; no starts or retries."""
from pathlib import Path
import json,subprocess
from run import O,read,progress

def main():
    state=read(O/'status.json');live=read(O/'live.json') if (O/'live.json').exists() else {}
    result=dict(batch_state=state['state'],terminal_cases=[dict(name=c['name'],returncode=c['returncode'],native_gate=c['native_gate']) for c in state['cases']],live=live,full_C1_complete=False)
    if live.get('state')=='running':
        pid=int(live['native_pid'])
        cmd=['powershell','-NoProfile','-Command',f'[Console]::OutputEncoding = [System.Text.UTF8Encoding]::new(); Get-CimInstance Win32_Process -Filter "ProcessId={pid}" | Select-Object ProcessId,CommandLine | ConvertTo-Json -Compress']
        query=subprocess.run(cmd,capture_output=True,text=True,encoding='utf-8',errors='replace',timeout=20,creationflags=subprocess.CREATE_NO_WINDOW)
        result['query_returncode']=query.returncode
        try:proc=json.loads(query.stdout)
        except json.JSONDecodeError:proc=None
        result['observed_process']=proc
        result['identity_matches']=bool(proc and proc['ProcessId']==pid and 'tf_retained_face_bench.exe' in proc['CommandLine'] and live['output_directory'] in proc['CommandLine'])
        result.update(progress(Path(live['output_directory'])))
    print(json.dumps(result,indent=2,ensure_ascii=True))

if __name__=='__main__':main()
