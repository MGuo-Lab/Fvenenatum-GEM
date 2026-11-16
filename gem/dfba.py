"""Dynamic flux balance analysis (dFBA) over several substrates.

``run_dfba_multi`` advances biomass and extracellular concentrations by an explicit Euler step, solving an FBA
problem at every step. Substrate uptake bounds follow Michaelis-Menten kinetics, an optional rule blocks every
substrate other than a primary one (glucose) while the primary is above a threshold, and an optional function bounds
specific O2 uptake, for example by the maximum oxygen transfer rate divided by the biomass concentration.

The analyses in ``analysis/`` call it through ``fbaval.simulate.run_dfba``.
"""

from __future__ import annotations

import cobra

from dfba_kinetics import (
    CardinalPH,
    OverflowParams,
    apply_overflow_lb,
    calibrate_net_cation,
    make_acid,
    reset_overflow_lb,
)


def mm_kinetics(vmax: float, km: float, conc_mm: float) -> float:
    if conc_mm <= 0:
        return 0.0
    return vmax * conc_mm / (km + conc_mm)


# ── dFBA loop ────────────────────────────────────────────────────────────────

def _organic_export_weights(model: cobra.Model, biomass_id: str) -> dict:
    """Export variable -> C atoms per mmol, for every open boundary except
    CO2 and biomass. A metabolite without a formula counts as 1 C."""
    weights = {}
    for r in model.boundary:
        if r.id == biomass_id or r.upper_bound <= 0:
            continue
        met = next(iter(r.metabolites))
        el = met.elements if met.formula else {"C": 1}
        if el.get("C", 0) > 0 and not set(el) <= {"C", "O"}:
            weights[r.forward_variable] = float(el["C"])
    return weights


def _nutrient_need_per_mu(model, biomass_rxn, mu, rids) -> dict[str, float]:
    """Least uptake of each mineral-nutrient exchange that growth at `mu`
    requires, per unit growth (the biomass composition, found by LP)."""
    need = {}
    with model:
        biomass_rxn.bounds = (mu, mu)
        for rid in rids:
            r = model.reactions.get_by_id(rid)
            model.objective = model.problem.Objective(r.flux_expression, direction="max")
            v = model.slim_optimize(error_value=float("nan"))
            if v == v:
                need[rid] = max(0.0, -v) / mu
    return need


def _growth_capped_fluxes(model, biomass_rxn, substrates, mu, export_c, caps=None):
    """Fluxes with uptake held at its kinetic rate (the current lower bound)
    and growth fixed at `mu`, choosing the distribution that exports the
    least organic carbon; mineral uptake held to ``caps`` (exchange id ->
    largest uptake). None if infeasible."""
    tol = 1e-9
    with model:
        for rid, cap in (caps or {}).items():
            r = model.reactions.get_by_id(rid)
            r.lower_bound = max(r.lower_bound, -cap * (1 + 1e-6) - tol)
        for s in substrates:
            if not s.get("on_demand"):
                v = s["rxn"].lower_bound
                s["rxn"].bounds = (v - tol, v + tol)
        biomass_rxn.bounds = (max(mu - tol, 0.0), mu + tol)
        model.objective = model.problem.Objective(
            sum(c * var for var, c in export_c.items()), direction="min")
        sol2 = model.optimize()
        if sol2.status != "optimal":
            return None
        return sol2.fluxes


def run_dfba_multi(
    model: cobra.Model,
    biomass_rxn: cobra.Reaction,
    substrates: list[dict],
    biomass_0: float,
    volume: float,
    t_end: float,
    dt: float,
    byproduct_ids: list[str],
    diauxic_primary_id: str | None = None,
    diauxic_threshold_mm: float = 0.5,
    pH_params: CardinalPH | None = None,
    overflow_params: OverflowParams | None = None,
    overflow_byproduct_keys: dict[str, str] | None = None,
    buffer_citrate_mM: float = 0.0,
    initial_pH: float = 6.5,
    medium_phosphate_mM: float = 50.0,
    primary_substrate_key: str | None = None,
    mu_max_lit: float = 0.30,
    o2_cap_fn=None,
    o2_rxn_id: str = "BOUNDARY_O2e",
    continue_on_infeasible: bool = False,
    acid_keys: dict[str, str] | None = None,
    nh4_rxn_id: str | None = "BOUNDARY_NH4e",
    pi_rxn_id: str | None = "BOUNDARY_PIe",
    so4_rxn_id: str | None = "BOUNDARY_SLFe",
    ionic_strength_M: float = 0.0,
    temp_c: float = 25.0,
    pH_coupling: str = "scale_all",
) -> list[dict]:
    """Generalised dFBA across an arbitrary list of substrates.

    Each entry of ``substrates`` is a dict with keys: ``rxn``, ``mw``,
    ``vmax``, ``km``, ``conc_mm`` (initial). The ``rxn.lower_bound`` is
    set every step from MM kinetics. Diauxic shift, when enabled, blocks
    every substrate other than ``diauxic_primary_id`` (on-demand nutrients
    excepted) whenever that primary is above ``diauxic_threshold_mm``.

    Optional kinetic add-ons (`dfba_kinetics`):

    * ``pH_params`` — when given, every step the medium pH is solved from
      acid speciation + buffer + phosphate, and Rosso 1995 cardinal-pH
      γ(pH) ∈ [0,1] slows the culture. ``pH_coupling`` picks how:

      - ``"scale_all"`` (default, the original behaviour): the whole flux
        vector is scaled by γ, uptake included, so acid stress slows
        consumption as much as growth.
      - ``"growth_cap"``: only growth is slowed. Substrate uptake is fixed at
        its kinetic rate, growth is fixed at γ·μ_fba, and a second LP
        minimises the organic carbon exported (C atoms, CO₂ and biomass
        excluded). Carbon that does not go to biomass is respired as far as
        O₂ allows and spilled as byproducts beyond that (uncoupled
        catabolism; Hanly & Henson 2011 BBE 108:376). Only carbon is
        uncoupled: uptake of the mineral nutrients (``nh4_rxn_id``,
        ``pi_rxn_id``, ``so4_rxn_id``) is capped at what the growth of the
        step needs, the minimal uptake per unit growth found once by LP, so a
        nutrient cannot be taken up only to dispose of electrons (e.g.
        sulfate reduced to sulfite and exported). This also holds when
        growth is limited by something other than pH (γ = 1 without
        ``pH_params``), e.g. a nutrient running out. A substrate dict with
        ``"on_demand": True`` (a nutrient such as NH4) is only bounded by
        its kinetic rate, not forced. If the second LP fails the step falls
        back to ``"scale_all"`` and the row carries ``coupling_fallback``.

      ``mu_max_lit`` is not used by either mode.

    * ``overflow_params`` + ``overflow_byproduct_keys`` — when given,
      every step a *minimum* secretion flux is imposed on each fermentative
      byproduct exchange, proportional to the previous step's primary
      substrate uptake (Vemuri 2007 / Basan 2015 phenomenology).

    * ``buffer_citrate_mM`` / ``medium_phosphate_mM`` / ``initial_pH`` —
      define the medium's buffer chemistry.  Net cation charge is
      calibrated *once* at t=0 to match ``initial_pH`` against the
      buffer-only acids; that calibrated cation is held constant for the
      run (no base feed model).

    * ``o2_cap_fn(t_h, biomass_gL) -> qO2_max`` — when given, the O2
      exchange lower bound is set every step to ``-qO2_max``
      (mmol/gDW/h), e.g. ``OTR_max / X`` for a kLa-limited vessel.

    * ``continue_on_infeasible`` — when True, a step whose LP is not
      optimal (e.g. maintenance cannot be met after substrate depletion)
      is recorded with zero fluxes and integration continues to ``t_end``
      instead of stopping.

    * ``acid_keys`` — exchange id -> acid name (``PKA_TABLE`` key) for the
      byproducts that count in the pH balance. Without it the overflow
      mapping is used, as before, so acids only counted with overflow on.

    * ``nh4_rxn_id`` / ``pi_rxn_id`` / ``so4_rxn_id`` — uptake exchanges
      whose cumulative assimilation shifts the pH balance (NH4+ acidifies,
      SO4^2- uptake alkalinises, phosphate leaves the buffer). Ignored when
      None or absent from the model. ``ionic_strength_M`` / ``temp_c``
      switch on Davies-corrected pKas.
    """
    if pH_coupling not in ("scale_all", "growth_cap"):
        raise ValueError(f"unknown pH_coupling {pH_coupling!r}")
    biomass = biomass_0
    t = 0.0
    results: list[dict] = []
    export_c = (_organic_export_weights(model, biomass_rxn.id)
                if pH_coupling == "growth_cap" else {})
    need_per_mu: dict[str, float] | None = None   # mineral uptake per unit growth

    byprod_conc = {rid: 0.0 for rid in byproduct_ids}

    # ── pH state setup ───────────────────────────────────────────────────
    pH_active = pH_params is not None
    overflow_active = (overflow_params is not None
                       and overflow_byproduct_keys is not None)

    # rid → human-readable acid key for pH calculation (acetate/pyruvate/…)
    rid_to_acid_key: dict[str, str] = {}
    if acid_keys is not None:
        rid_to_acid_key = dict(acid_keys)
    elif overflow_active:
        for k, rid in overflow_byproduct_keys.items():
            if rid is not None:
                rid_to_acid_key[rid] = k

    # cumulative nutrient assimilation (mM of medium) for the pH balance
    uptake_ids = {k: rid for k, rid in (("nh4", nh4_rxn_id), ("pi", pi_rxn_id),
                                         ("so4", so4_rxn_id))
                  if rid is not None and rid in model.reactions}
    assimilated = {k: 0.0 for k in uptake_ids}

    net_cation_mM = 0.0
    if pH_active:
        # baseline acids in fresh medium = buffer + phosphate only
        base_acids = []
        if buffer_citrate_mM > 0:
            base_acids.append(make_acid("citrate", buffer_citrate_mM, ionic_strength_M, temp_c))
        if medium_phosphate_mM > 0:
            base_acids.append(make_acid("phosphate", medium_phosphate_mM, ionic_strength_M, temp_c))
        base_acids = [a for a in base_acids if a is not None]
        net_cation_mM = calibrate_net_cation(initial_pH, base_acids)

    # ── primary substrate for overflow trigger ──────────────────────────
    primary_rxn_id: str | None = None
    if overflow_active:
        if primary_substrate_key is not None:
            for s in substrates:
                if s["name"] == primary_substrate_key:
                    primary_rxn_id = s["rxn"].id
                    break
        if primary_rxn_id is None and substrates:
            # fall back: the first substrate (e.g. glucose-only run)
            primary_rxn_id = substrates[0]["rxn"].id

    def medium_pH() -> float:
        """pH of the medium as it stands (acids secreted, nutrients taken up)."""
        current_byprod_mM: dict[str, float] = {}
        for rid, conc in byprod_conc.items():
            key = rid_to_acid_key.get(rid)
            if key:
                current_byprod_mM[key] = current_byprod_mM.get(key, 0.0) + max(0.0, conc)
        from dfba_kinetics import step_pH_state  # local to avoid cycle
        return step_pH_state(current_byprod_mM, buffer_citrate_mM,
                             medium_phosphate_mM, net_cation_mM,
                             ammonium_assimilated_mM=assimilated.get("nh4", 0.0),
                             phosphate_assimilated_mM=assimilated.get("pi", 0.0),
                             sulfate_assimilated_mM=assimilated.get("so4", 0.0),
                             ionic_strength_M=ionic_strength_M, temp_c=temp_c)

    last_q_S = 0.0  # |uptake| of the primary substrate in the previous LP
    while t <= t_end:
        # ── 1. uptake bounds (MM kinetics + diauxic) ───────────────────
        primary_above = False
        if diauxic_primary_id is not None:
            for s in substrates:
                if s["rxn"].id == diauxic_primary_id:
                    primary_above = s["conc_mm"] > diauxic_threshold_mm
                    break
        for s in substrates:
            # MM rate, capped so one Euler step cannot take more substrate
            # than is left (otherwise the clamp below creates carbon)
            v = min(mm_kinetics(s["vmax"], s["km"], s["conc_mm"]),
                    s["conc_mm"] * volume / (biomass * dt))
            if (diauxic_primary_id is not None and primary_above
                    and s["rxn"].id != diauxic_primary_id and not s.get("on_demand")):
                s["rxn"].lower_bound = 0.0
            else:
                s["rxn"].lower_bound = -v

        # ── 2. overflow byproduct lb (uses last step's q_S) ─────────────
        if overflow_active:
            apply_overflow_lb(model, overflow_byproduct_keys,
                              last_q_S, overflow_params)

        if o2_cap_fn is not None:
            model.reactions.get_by_id(o2_rxn_id).lower_bound = -abs(
                o2_cap_fn(t, biomass / volume))

        # ── 3. solve LP ────────────────────────────────────────────────
        sol = model.optimize()
        if sol.status != "optimal":
            # If overflow lb made the LP infeasible (e.g. uptake=0 but
            # forced byproduct flux > 0), retry without overflow once.
            if overflow_active:
                reset_overflow_lb(model, overflow_byproduct_keys)
                sol = model.optimize()
            if sol.status != "optimal":
                if not continue_on_infeasible:
                    print(f"t={t:.2f}h: solver status '{sol.status}', stopping.")
                    break
                row = {"time": round(t, 4), "biomass": biomass, "mu": 0.0,
                       "mu_fba": 0.0, "pH": medium_pH() if pH_active else float("nan"),
                       "gamma_pH": 1.0, "infeasible": True,
                       "nh4_assimilated_mM": assimilated.get("nh4", 0.0)}
                for s in substrates:
                    row[f"v_{s['rxn'].id}"] = 0.0
                    row[f"c_{s['rxn'].id}"] = s["conc_mm"]
                for rid in byproduct_ids:
                    row[f"v_{rid}"] = 0.0
                    row[f"c_{rid}"] = byprod_conc[rid]
                results.append(row)
                t += dt
                continue

        mu_fba = sol.fluxes[biomass_rxn.id]

        # ── 4. pH from current medium state, then γ(pH) multiplier ──────
        gamma = 1.0
        pH_now = float("nan")
        if pH_active:
            pH_now = medium_pH()
            gamma = pH_params.gamma(pH_now)

        mu_eff = mu_fba * gamma

        row = {"time": round(t, 4), "biomass": biomass, "mu": mu_eff,
               "mu_fba": mu_fba, "pH": pH_now, "gamma_pH": gamma,
               "nh4_assimilated_mM": assimilated.get("nh4", 0.0),
               "pi_assimilated_mM": assimilated.get("pi", 0.0),
               "so4_assimilated_mM": assimilated.get("so4", 0.0)}

        # flux vector for the Euler step and the factor applied to it
        fluxes, scale = sol.fluxes, gamma
        if pH_coupling == "growth_cap":
            if need_per_mu is None and mu_eff > 1e-9:
                need_per_mu = _nutrient_need_per_mu(model, biomass_rxn, mu_eff,
                                                    list(uptake_ids.values()))
            caps = {rid: n * mu_eff for rid, n in (need_per_mu or {}).items()}
            capped = _growth_capped_fluxes(model, biomass_rxn, substrates,
                                           mu_eff, export_c, caps)
            if capped is None:
                row["coupling_fallback"] = True
            else:
                fluxes, scale = capped, 1.0

        # nutrient uptake this step (uptake fluxes are negative)
        for k, rid in uptake_ids.items():
            assimilated[k] += max(0.0, -fluxes[rid]) * biomass * dt * scale / volume

        # Substrate flux + concentration
        for s in substrates:
            v = fluxes[s["rxn"].id]
            row[f"v_{s['rxn'].id}"] = v
            row[f"c_{s['rxn'].id}"] = s["conc_mm"]

        # Byproduct cumulative concentration (scaled the same way as the
        # biomass update, so stoichiometry stays consistent)
        for rid in byproduct_ids:
            v = fluxes[rid]
            byprod_conc[rid] += v * biomass * dt * scale / volume
            row[f"v_{rid}"] = v
            row[f"c_{rid}"] = byprod_conc[rid]

        results.append(row)

        # ── 5. Euler update (every flux scaled by `scale`) ──────────────
        # Substrates and byproducts use the same pre-step biomass as the
        # growth term, so the update is mass-balanced for any dt.
        for s in substrates:
            v = fluxes[s["rxn"].id]
            s["conc_mm"] += v * biomass * dt * scale / volume
            s["conc_mm"] = max(s["conc_mm"], 0.0)
        biomass += mu_eff * biomass * dt

        # Track q_S for next-step overflow lb
        if primary_rxn_id is not None:
            last_q_S = abs(fluxes[primary_rxn_id])

        t += dt

    return results
