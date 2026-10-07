# SAD sizing workflow

Open `run_sad_sizing.jl` in VS Code and select **Run SAD sizing workflow** from the Run and Debug menu. The runner activates the TASOPT project automatically, so it does not depend on the terminal's working directory.

Open `run_sad_sizing.jl` in VS Code and edit the model lists and `CONFIG` near the top.
The current optimizer limits and sweep grids are reduced for quick workflow tests.

- Every model in `BUILD_MODELS` is rebuilt and overwrites its existing TOML.
- An omitted prerequisite is loaded from `models`; the run stops if it is missing.
- An empty payload-analysis list disables that analysis.
- `SUMMARY_MODELS` controls the comparison rows independently of the sweep lists.
- `RUN_ENGINE_DECK` controls the optional CFM56 engine-deck export.
- `RUN_SUMMARY` controls the final comparison CSV.

## SO2 payload tank in the radius trade

The radius-trade payload grid is SO2 mass, not total payload mass. For every
grid point, the workflow sizes one longitudinal aluminum tank with hemispherical
ends, adds the calculated tank mass to the mission payload, and checks the
combined mass against the radius-dependent structural payload capacity. It also
requires the tank and its radial installation clearance to fit inside the
fuselage skin and requires the tank length to fit within the cylindrical cabin.

The tank inputs and the required design SO2 payload are configured under
`CONFIG.radius_trade` in `run_sad_sizing.jl`. Radius points that cannot carry the
configured design SO2 payload plus its tank are excluded when selecting the
reduced-radius model. Among the feasible radius points, the selected aircraft is
the one with the highest service ceiling at the design SO2 payload; lower MTOW
breaks ties within the altitude-search tolerance. The downstream payload-Mach
sweep also treats its payload grid as SO2 mass and adds the corresponding tank
mass before each mission evaluation. Liquid properties are specified at 303.15 K
with 4% ullage. The design pressure is the NIST Antoine-correlation saturation
pressure plus 5 bar. The aluminum barrel and hemispherical heads are sized from
their pressure membrane stresses using the configured allowable stress and weld
efficiency. The aspect ratio increases from its configured minimum when required
to satisfy fuselage diameter and radial clearance. Tank mass is bare shell mass;
no installation allowance is applied.

Run `run_sad_sizing.jl` directly from VS Code. The workflow always applies composites and interior removal before installing the CFM56.

Plots remain separate:

```powershell
python "sad sizing/plot_sweeps.py"
python "sad sizing/plot_payload_range.py"
```

Generated TOMLs are stored in `models`. CSV outputs and figures are stored under `results`.
