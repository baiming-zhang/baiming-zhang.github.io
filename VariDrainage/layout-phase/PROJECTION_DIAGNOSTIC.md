# Initial field representation check, observer only

Registered after initial_continuum.json, before this new computation.
The prior nodal constrained-preparation diagnostic passed its analytic ODE
identity; nodal edge RMS difference is about0.71nm in all recorded layouts.
Do not overwrite that diagnostic, its plan hash or the native state.

For each frozen background layout and guarded split20, compare INITIAL h
against uniform20 using the exact original film-cell overlap geometry.
Define observed delta = h_native_a - h_native_b on each overlap; define a
sampling counterpart delta_c = h_circle(x_a)-h_circle(x_b), using the model's
analytic preparation circle only in this observer. Report RMS(delta),
RMS(delta_c), RMS(delta-delta_c) and their norm triangle bound.

The final remainder isolates deviations from the sampled preparation circle
within this algebraic decomposition. It is not a model change, an acceptance
gate, a dynamical causation test, or a continuum proof. Never replace native
initial h by circle values in the running native experiment. All old gates and
failed diagnostics remain unchanged. This helps distinguish a nodal shape
diagnostic from comparison of discontinuous cell representations.

Also verify the logical initial-site mapping in x<180um, checking position,
h,V,S,theta before reporting matched-node dynamic diagnostics. Values that
only agree within floating precision must not be called byte-identical.
