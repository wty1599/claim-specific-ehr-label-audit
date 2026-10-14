# Calibration-definition addendum

Date: 2026-10-08. Recorded after inspection of the archived calibration helper
and before corrected-temperature Raw-model evaluation. This is a clarification
of the authorized post hoc completion, not a new analysis standard.

The initial specification remains unmodified:
`RAW_MODEL_COMPLETION_SPEC.md`, SHA-256
`2CA190DCAF14FEE54A76606C73CFA80D9A872E13259B408A30BEE5799789E905`.

Its phrase "joint calibration slope/intercept" is imprecise for the archived
function. With LP = logit(p), using the original probability clamp of 1e-6:

- `calibration_slope` is the slope from `glm(y ~ LP, family = binomial())`,
  with a jointly estimated intercept.
- `calibration_intercept` is the intercept from
  `glm(y ~ 1, offset = LP, family = binomial())`: calibration-in-the-large
  (CITL), with the LP coefficient fixed at one.

The completion preserves these original definitions and the Figure 4 display
convention. The separately named `calibration_joint_intercept` field may remain
in the aggregate output but must not replace the original displayed intercept.
No additional intervals, model selection, or new decision standard are added.
The requested patient bootstrap AUC intervals and paired AUC contrasts are
unchanged. This addendum does not modify the initial specification or its hash.
