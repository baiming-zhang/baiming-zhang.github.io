# Independent audits of the new 5um terminal cases

Prepared before any new5um case is terminal. Reuse retained_reference's audit
wrapper with exactly these changes: campaign input namespace and auditor module
name, import the frozen runner ENV through run.base, per-case cap600s, planned
case count2 and matching coverage key. Native formulas and gates unchanged.
The underlying global_space/audit_ready.py remains unmodified, with the same
explicit immutable-snapshot path and reused24candidate-fixture metadata adapter.
No new native executions or fixture reruns in this processing stage.

Require exact terminal passing record,N2112,160/320steps,T10us, complete file set
and all raw/closed-log hashes. Refuse live cases before reading frame contents.
Audit only disjoint immutable terminal copies. One600s attempt per case; no
automatic retry of an unsealed attempt. Successful seals are hash-verified and
reused.18scalar+18mechanical saved audits expected for the complete pair, not
assumed. Partial coverage cannot qualify the pair or the failed case.

Full pair comparisons wait for BOTH new native cases and independent audits.
Use original moving-overlap metric/time2%/space10% and same-candidate10um
references at matching dt. Retain all failures and original50/200us/physics/
exchange/TFE/global/experimental/fullC1 requirements.
