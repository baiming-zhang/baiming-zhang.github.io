# Monitored actual uniform5 pilot launch

Plan/driver:1b1c2d8, frozen protocol:6f68d6a. Qualified operator/control results:
423f019; original assembly source5d00f7e; control drivere098cb3.
Native executable tf_local_assembly_bench, SHA256 in protocol.json.
Launched ONLY after native automation create/view and ACTIVE30min config
verification. monitor_at_launch.json and launch_snapshot.json are immutable
snapshots. Local status.json/live.json/monitor JSONL remain mutable and are
NOT copied into the public page as finished data. No live frame archiving.

Local monitor: python thin_film/scripts/uniform5_pilot/monitor.py(read-only).
Live native PTY handle is retained by the main agent; never restart on a stale
file, delayed output or observation timeout. Parent Popen handle is authority.
Captured headers/CPP/driver/binary must remain unchanged during native work;
new independent Python analyses and result documentation can be added safely.
All scientific acceptance remains open until terminal data and audits prove it.
