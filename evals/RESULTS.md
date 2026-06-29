# Eval results — candidate-gene interpretation

Faithfulness + calibration of the candidate-gene interpretation step (**not** causal
correctness — see [`README.md`](README.md)). Higher *gene existence* = fewer
hallucinated genes; lower *citation fabrication* = fewer non-resolving PMIDs; higher
*calibration* = better hedging of an under-determined call. Means are over items the
model actually answered; refusals and unparseable replies are reported separately, not
scored as faithfulness failures.

## Summary

| model | answered | refused | gene existence | citation fabrication | calibration |
|---|---|---|---|---|---|
| `claude-opus-4-8` | 4/4 | 0 | 1.00 | 0.15 | 0.81 |
| `claude-sonnet-4-6` | 2/4 | 2 | 1.00 | 0.38 | 1.00 |
| `claude-haiku-4-5-20251001` | 4/4 | 0 | 1.00 | 0.41 | 1.00 |

## Per item

### `claude-opus-4-8`

| item | genes | status | existence | fabrication | calibration |
|---|---|---|---|---|---|
| eb9 | 50 | ok | 1.00 | 0.40 | 1.00 |
| eb5 | 70 | ok | 1.00 | 0.20 | 1.00 |
| bs4 | 67 | ok | 1.00 | 0.00 | 0.50 |
| fw22 | 62 | ok | 1.00 | 0.00 | 0.75 |

### `claude-sonnet-4-6`

| item | genes | status | existence | fabrication | calibration |
|---|---|---|---|---|---|
| eb9 | 50 | ok | 1.00 | 0.42 | 1.00 |
| eb5 | 70 | no_response | — | — | — |
| bs4 | 67 | no_response | — | — | — |
| fw22 | 62 | ok | 1.00 | 0.33 | 1.00 |

### `claude-haiku-4-5-20251001`

| item | genes | status | existence | fabrication | calibration |
|---|---|---|---|---|---|
| eb9 | 50 | ok | 1.00 | 0.33 | 1.00 |
| eb5 | 70 | ok | 1.00 | 0.38 | 1.00 |
| bs4 | 67 | ok | 1.00 | 0.43 | 1.00 |
| fw22 | 62 | ok | 1.00 | 0.50 | 1.00 |

