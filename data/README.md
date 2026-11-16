# Batch fermentation data

Batch fermentations of *Fusarium venenatum* A3/5 on glucose and fructose, published in Vinestock & Guo (2026),
Experiments & modelling of batch fermentation of *Fusarium venenatum* on glucose-fructose mixtures,
*Syst. Control Trans.* 5, 776–782, https://doi.org/10.69997/sct.178293.

## Culture conditions

- 10,000 spores per mL in 50 mL of Rank Hovis McDougall (RHM) medium in 250 mL baffled shake flasks, 28 °C, 220 rpm.
- Medium per flask: 46 mL RHM made with 30 g L⁻¹ KH2PO4 (other salts, trace elements and biotin as in RHM: NH4Cl
  4.4 g L⁻¹, K2SO4 0.3 g L⁻¹, MgSO4·7H2O 0.25 g L⁻¹), 0.5 mL 5 % polypropylene glycol, 0.5 mL 1 % Tween 80 and
  3 mL 25 % sugar stock. Measured initial pH 6.5–6.7.
- Humidified headspace with 200 % of the O2 of air; kLa 96 h⁻¹, which gives a maximum oxygen transfer rate of
  46.9 mmol L⁻¹ h⁻¹.
- Four conditions with 15 g L⁻¹ total sugar, three flasks (A1–A3) each:

| Folder | Condition | Glucose (g L⁻¹) | Fructose (g L⁻¹) |
|---|---|---|---|
| `flask/Glc15` | Glc 15 | 15 | 0 |
| `flask/Glc10Fru5` | Glc 10 + Fru 5 | 10 | 5 |
| `flask/Glc5Fru10` | Glc 5 + Fru 10 | 5 | 10 |
| `flask/Fru15` | Fru 15 | 0 | 15 |

## Files

Each condition folder holds three workbooks.

- `workbook.xlsx`
  - `HPLC`: one row per sampling time point (TP0, TP1, ...), with date, elapsed time (h) and sampled volume, and for
    each flask the peak area and concentration (g L⁻¹) of glucose, fructose, glycerol, ethanol, acetate and
    "lactate".
  - `Calibration`: HPLC calibration factors (peak area per g L⁻¹) of glycerol, ethanol, acetate and lactic acid.
  - `Parameters`: culture settings (temperature, headspace, shaking, volume, initial sugars, kLa, injection volume).
  - `Initial Data`, `Final Data`: pH, broth mass and volume, volume removed by sampling, dry biomass (mg),
    biomass concentration (g L⁻¹) and yield (g g⁻¹) of each flask.
- `Processed_Presens_Data.xlsx`: online readings of the PreSens shake flask readers, about every 15 min, positions
  1–3 being flasks A1–A3: dissolved O2 (% of air saturation), pH, biomass (g L⁻¹) and O2 uptake rate
  (g L⁻¹ h⁻¹), with absolute timestamps. Biomass is inferred from reflectance and calibrated to the start and end
  points of each flask; pH is calibrated to the start and end values.
- `Processed_HPLC_Data.xlsx`: the HPLC sampling times, used to put online readings and HPLC samples on one time axis.

The tidy CSV files in `flask/` hold the same measurements for all four conditions, with time in hours from the first
HPLC sample of each fermentation. They are written by `analysis/scripts/export_flask_data.py`.

| File | Columns |
|---|---|
| `hplc.csv` | `condition`, `flask`, `time_h`, `compound` (glucose, fructose, glycerol, ethanol, acetate), `conc_gl` |
| `online.csv` | `condition`, `flask`, `time_h`, `biomass_gl`, `ph`, `do_pct_as` (dissolved O2, % air saturation), `our_g_l_h` |
| `final.csv` | `condition`, `flask`, `biomass_gl` (end-point dry weight), `yield_gg`, `end_volume_ml`, `end_ph` |

## Notes

- HPLC: Bio-Rad Aminex HPX-87H, 0.001 M H2SO4 at 0.4 mL min⁻¹ and 30 °C, refractive-index detection. Sugar
  concentrations are peak areas scaled to the nominal initial concentration of each run. In `hplc.csv`, glycerol,
  ethanol and acetate are recomputed from peak areas with the calibration factor of each compound.
- The workbook column "lactate" holds a peak at 16.1 min, which is not lactic acid and is tentatively assigned to
  arabitol (Supplementary Note 3 of the paper). Its values in the workbooks use the lactic acid calibration and are
  therefore not concentrations. The re-integrated peak areas are in `results/flask_hplc/unidentified.csv`, and the
  peak is left out of `hplc.csv`.
- Yield is the dry biomass divided by the effective final volume (initial volume minus all sampled volumes) and by
  the nominal initial sugar concentration.
- Online readings start up to 1 h after the first HPLC sample. Time in the CSV files and in the analyses is computed
  from the absolute timestamps, so both data types share the same zero.
