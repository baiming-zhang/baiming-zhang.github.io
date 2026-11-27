"""Start one bounded, resumable timing waiter; preserve an existing one."""
from pathlib import Path
import argparse, datetime, fcntl, json, os, subprocess, sys

R=Path(__file__).resolve().parent
S=R/'stage3_runtime_symmetry'

def main():
    assert os.name=='posix','Run this small launcher with the Linux experiment Python'
    parser=argparse.ArgumentParser();parser.add_argument('--stable',action='store_true');args=parser.parse_args()
    phase='runtime_stable' if args.stable else 'runtime_formal'
    OUT=S/(phase+'_results')
    OUT.mkdir(parents=True,exist_ok=True)
    guard=(OUT/'launcher.lock').open('a+')
    try:fcntl.flock(guard,fcntl.LOCK_EX|fcntl.LOCK_NB)
    except BlockingIOError:print('Existing launcher retained');return
    done=OUT/'completion.json'
    if done.exists() and json.loads(done.read_text()).get('state')=='complete':
        print(json.dumps(dict(state='complete',new_process_started=False)));return
    pidfile=S/(phase+'.pid')
    if pidfile.exists():
        pid=int(pidfile.read_text())
        command=Path(f'/proc/{pid}/cmdline')
        if command.exists() and phase+'_waiter.py' in command.read_bytes().decode(errors='replace'):
            print(json.dumps(dict(state='existing_waiter',pid=pid,new_process_started=False)));return
    script=S/(phase+'_waiter.py')
    assert script.exists() and (R/'deployment_selection_locked.json').exists()
    with (OUT/'controller.out').open('ab') as log,(OUT/'controller.err').open('ab') as err:
        proc=subprocess.Popen([sys.executable,'-B','-u',str(script),'--max-wait-seconds','300'],stdout=log,stderr=err,stdin=subprocess.DEVNULL,start_new_session=True,cwd=R.parent)
    pidfile.write_text(str(proc.pid))
    entry=dict(state='bounded_waiter_started',phase=phase,pid=proc.pid,max_wait_seconds=300,started=datetime.datetime.now(datetime.timezone(datetime.timedelta(hours=-7))).isoformat())
    with (OUT/'launch_history.jsonl').open('a') as log:log.write(json.dumps(entry)+'\n')
    print(json.dumps(entry))

if __name__=='__main__':main()
