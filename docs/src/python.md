# Python package `pyentropyscaling`

`pyentropyscaling` exposes the functionality of `EntropyScaling.jl` to Python via [`juliacall`](https://juliapy.github.io/PythonCall.jl/stable/juliacall/), which provides the bridge to Julia under the hood. Unless otherwise noted, the public Python API mirrors the Julia API described in this documentation.

## Installation

Install from PyPI:

```
pip install pyentropyscaling
```

## Example

```python
import pyentropyscaling as es
import pyclapeyron
import numpy as np

pure = es.RefpropRESModel("water")

eta = es.viscosity(pure, 1e5, 300)

mix = es.ChapmanEnskogModel(["butane", "methanol"])
lmbd = es.thermal_conductivity(mix, 300)
```

## Notes & limitations

- Inputs: arrays passed to `EntropyScaling.jl` functions must be NumPy arrays.
- Outputs: returned arrays are `juliacall` array objects; they can be converted to NumPy arrays with `np.array(...)` (see example above).
- Julia macros and mutating (in-place, `!`) functions are not directly callable from Python.
- Please report any issues on the `EntropyScaling.jl` GitHub repository.
