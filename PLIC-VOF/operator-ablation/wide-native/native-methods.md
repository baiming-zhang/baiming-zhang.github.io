# Native local-operator comparison

The primary comparison uses the original compiled Aphros and SimLiquid local operators versus the fixed neural weights compiled with TorchScript, all in float64 on the same pinned CPU core and one thread. The raw input is the identical float32-quantized VOF promoted to double. Normal estimation, PLIC reconstruction, fragment construction and particle fitting are inside the native operator boundary. This measures regular-interior local dependency slices, not complete simulators.

All 2,560 shuffled exclusive cases passed: 396,800 measured invocations and 51,200 warmups across N128/N512, two physical suites and five batch sizes. The 199 preliminary native-first CPU cases are excluded. No significant-preemption case needed retry. The three fixed anchor medians were 1.496236, 1.641562 and 1.822417 ms: 21.8% drift. Small timing differences are not a robust ordering; no normalization was invented.

The measured CPU64 neural implementation does not demonstrate speedup over these native local operators. For independent N512, batch4096, Aphros/SimLiquid flux cost 0.04276/0.1088 ms with flux-density MAE 0.00120311. The fastest three-seed group reaching that native error is control 11x11, 27.05 +/- 2.26 ms and MAE 0.00106158 +/- 0.0000418 (between-seed sample SD).

For independent N512 curvature, Aphros/SimLiquid cost 7.693/6.428 ms with physical curvature MAE 0.218996. None of the ten new three-seed groups meets that error. The mature original single fixed reference does meet it: MAE 0.187899263, batch4096 median 53.612645 ms, p05/p95 49.8781873/56.3643823 ms. Its batch1 median is 0.040958 ms. This single fixed reference has no between-seed SD. At N128 its MAE 0.07219657 does not reach the native 0.0616697 target.

These statements concern these measured implementations and hardware. They do not establish a universal ordering of neural and numerical algorithms. CPU32 and CUDA32 are separate deployment layers, and the mature ONNX retains native double intermediate arithmetic despite float32 inputs. No CPU/GPU ratio is presented as an algorithmic speedup.

The graph suite preserves 20 native particle failures per 3,100 original queries at each N. Full-task frontiers exclude incomplete outputs; conditional common-finite frontiers are separately labelled. No missing output is replaced by zero and no unreachable target is extrapolated. The independent suite has full native coverage.

Review `comparison/cpu64_key_scenes.csv` for the compact N128/N512, batch1/4096 table; `comparison/readable_native_iso_examples.json` for exact model/fixture/subset/N/batch keys and both historical-inclusive and three-seed-only comparisons; and `comparison/qualification_summary.json` for measurement qualification. The two compact figures and their captions are in `comparison/figures/`. All three seed values, call quantiles and adjacent stencil-width changes are retained in the CSV/JSON indices.
