"""Generate the data files for the cold-plate resistance figures (Figs 4.3, 4.4
and 4.8), including the propagated measurement uncertainty.

Run from this folder:  python make_rcp_data.py
Writes rcp_*.dat files read by figures/galinstan_rcp.tex, water_rcp.tex and
galinstan_rcp_cfd.tex, and prints a summary of the uncertainties.

Model: epsilon-NTU for an isothermal plate (Section 3.3.2), Shah & London mean
Nu for laminar flow and Gnielinski for turbulent water (Appendix B).
Uncertainty: first-order (root-sum-square) propagation of independent errors
(Appendix C).
"""
import math

# --- Properties and geometry (Table 3.5) -----------------------------------
GAL = dict(rho=6370.0, cp=366.0, k=16.5, mu=0.0024)
WAT = dict(rho=1000.0, cp=4182.0, k=0.60, mu=0.0010)
D, L_PASS, N_PASS = 0.008, 0.152, 4
A = math.pi * D * N_PASS * L_PASS            # heated area, 0.0153 m^2
R_PACK = 0.24199                             # calibrated package resistance (C/W)

# --- Measurement uncertainty (Appendix C) -----------------------------------
U_T = 2.2          # thermocouple accuracy per reading (C)
U_QIN = 0.02       # relative uncertainty of the supply power (assumed)
Q_NOFLOW = 12.0    # heat load of the no-flow tests that give R_misc (W)


def nu_laminar(xs):
    """Shah & London mean Nu, constant wall temperature."""
    if xs <= 0.005:
        return 1.615 * xs ** (-1 / 3) - 0.7
    if xs < 0.03:
        return 1.615 * xs ** (-1 / 3) - 0.2
    return 3.657 + 0.0499 / xs


def nu_gnielinski(re, pr):
    f = (0.790 * math.log(re) - 1.64) ** -2
    return (f / 8) * (re - 1000) * pr / (1 + 12.7 * (f / 8) ** 0.5 * (pr ** (2 / 3) - 1))


def r_cp_model(vdot_mls, fluid, regime=None):
    """Cold-plate resistance from the epsilon-NTU model (C/W)."""
    v = vdot_mls * 1e-6
    re = 4 * fluid["rho"] * v / (math.pi * D * fluid["mu"])
    pr = fluid["mu"] * fluid["cp"] / fluid["k"]
    if regime is None:
        regime = "lam" if re < 2300 else "turb" if re >= 4000 else None
    if regime == "lam":
        nu = nu_laminar(L_PASS / (D * re * pr))
    elif regime == "turb":
        nu = nu_gnielinski(re, pr)
    else:
        return None
    ha = nu * fluid["k"] / D * A
    mcp = fluid["rho"] * v * fluid["cp"]
    return 1 / (mcp * (1 - math.exp(-ha / mcp)))


def r_cp_ats(vdot_mls):
    """ATS-TCP-1001 datasheet fit for water (C/W)."""
    return 0.0881 * vdot_mls ** -0.415


# --- Measured steady-state data (Thermal_Network_Model.xlsx, Inputs) --------
# name, coolant, nominal flow (mL/s), T_in, T_out, T_res (mean L/R), T_A, R_misc, Q_in
TESTS = [
    ("gal_full", "gal", 1.4, 29.80, 47.63, 59.25, 20.50, 3.3267, 72.0),
    ("gal_low", "gal", 0.75, 26.50, 43.60, 51.10, 20.50, 3.3267, 40.0),
    ("wat_58", "wat", 58.3, 33.10, 33.50, 49.22, 21.57, 2.95, 72.0),
    ("wat_29", "wat", 29.2, 32.00, 32.70, 48.60, 21.57, 2.95, 72.0),
    ("wat_12", "wat", 12.2, 33.68, 33.31, 50.12, 21.57, 2.95, 72.0),
    ("wat_8", "wat", 8.3, 33.94, 33.91, 50.53, 21.57, 2.95, 72.0),
]


def propagate(t):
    name, cool, vnom, tin, tout, tres, ta, rmisc, qin = t
    u_dt = math.sqrt(2) * U_T                          # any temperature difference
    # heat lost to ambient and its uncertainty
    qmisc = (tres - ta) / rmisc
    u_rmisc_rel = math.hypot(u_dt / (rmisc * Q_NOFLOW), U_QIN)
    u_qmisc = qmisc * math.hypot(u_dt / (tres - ta), u_rmisc_rel)
    qplate = qin - qmisc
    u_qplate = math.hypot(U_QIN * qin, u_qmisc)
    # total resistance resistor-to-inlet and inferred cold-plate resistance
    dt = tres - tin
    r_tot = dt / qplate
    u_rtot = r_tot * math.hypot(u_dt / dt, u_qplate / qplate)
    r_cp = r_tot - R_PACK
    # flow: Galinstan from the energy balance, water from the meter
    if cool == "gal":
        fl = GAL
        dtc = tout - tin
        v = qplate / (fl["rho"] * fl["cp"] * dtc) * 1e6
        u_v = v * math.hypot(u_qplate / qplate, u_dt / dtc)
    else:
        v, u_v = vnom, 0.0
    return dict(name=name, cool=cool, v=v, u_v=u_v, qplate=qplate,
                u_qplate=u_qplate, r_tot=r_tot, u_rtot=u_rtot, r_cp=r_cp,
                u_rcp=u_rtot)


def write(fname, header, rows):
    with open(fname, "w") as f:
        f.write(header + "\n")
        for r in rows:
            f.write(" ".join(f"{x:.5g}" for x in r) + "\n")


def sweep(lo, hi, n, fn):
    out = []
    for i in range(n + 1):
        v = lo * (hi / lo) ** (i / n)
        r = fn(v)
        if r is not None:
            out.append((v, r))
    return out


if __name__ == "__main__":
    res = [propagate(t) for t in TESTS]
    gal = [r for r in res if r["cool"] == "gal"]
    wat = [r for r in res if r["cool"] == "wat"]
    write("rcp_meas_galinstan.dat", "v rcp u_rcp u_v",
          [(r["v"], r["r_cp"], r["u_rcp"], r["u_v"]) for r in gal])
    write("rcp_meas_water.dat", "v rcp u_rcp",
          [(r["v"], r["r_cp"], r["u_rcp"]) for r in wat])
    write("rcp_model_galinstan.dat", "v rcp",
          sweep(0.3, 5.5, 120, lambda v: r_cp_model(v, GAL)))
    v_lam_max = 2300 * math.pi * D * WAT["mu"] / (4 * WAT["rho"]) * 1e6
    v_turb_min = 4000 * math.pi * D * WAT["mu"] / (4 * WAT["rho"]) * 1e6
    write("rcp_model_water_lam.dat", "v rcp",
          sweep(5.0, v_lam_max * 0.9999, 60, lambda v: r_cp_model(v, WAT, "lam")))
    write("rcp_model_water_turb.dat", "v rcp",
          sweep(v_turb_min, 75.0, 80, lambda v: r_cp_model(v, WAT, "turb")))
    write("rcp_ats_fit.dat", "v rcp", sweep(31.5, 75.0, 40, r_cp_ats))
    write("rcp_ats_extrap.dat", "v rcp", sweep(5.0, 31.5, 40, r_cp_ats))
    write("rcp_ats_points.dat", "v rcp",
          [(1.89271 * 1000 / 60, 0.021), (3.78541 * 1000 / 60, 0.016)])
    # CHT results, final iteration on the 427k-cell mesh
    write("rcp_cfd.dat", "v rcp",
          [(0.75, 0.5615), (1.4, 0.2820), (2.8, 0.1181), (5.0, 0.0499)])

    print(f"A = {A:.5f} m^2; laminar water up to {v_lam_max:.2f} mL/s, "
          f"turbulent from {v_turb_min:.2f} mL/s")
    print("test        V(mL/s)        Q_plate (W)    R_tot (C/W)      R_cp (C/W)")
    for r in res:
        print(f"{r['name']:9s} {r['v']:6.2f}+-{r['u_v']:4.2f}  "
              f"{r['qplate']:5.1f}+-{r['u_qplate']:3.1f}  "
              f"{r['r_tot']:.3f}+-{r['u_rtot']:.3f}  "
              f"{r['r_cp']:.3f}+-{r['u_rcp']:.3f}")
    for v in (0.75, 1.4, 2.8, 5.0):
        print(f"model Galinstan {v} mL/s: {r_cp_model(v, GAL):.4f}")
    for v in (8.3, 12.2):
        print(f"model laminar water {v} mL/s: {r_cp_model(v, WAT, 'lam'):.4f}")
