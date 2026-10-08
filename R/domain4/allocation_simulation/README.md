# D4 allocation-only simulation

The scripts `diagnostics.R`, `run_simulation.R`, `test_diagnostics.R`, and
`verify_outputs.R` implement the later post hoc component simulation.
They generate allocation data, not mortality or RRT effect estimates. The
31 settings use 500 formal repeats each; failures and zero-treated repeats
remain in the denominator. Public aggregate outputs are under
`data/simulation/domain4_allocation/`.

From this directory, run `Rscript test_diagnostics.R`,
`Rscript run_simulation.R --smoke`, `Rscript run_simulation.R --full`, and
`Rscript verify_outputs.R`. The scripts write a local `output/` directory;
do not commit that directory without a new aggregate-only review.

This does not replace the historical D4 G5/G6 longitudinal controls. It
does not evaluate the full support criterion involving 30-day deaths.
