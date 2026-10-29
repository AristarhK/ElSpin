"""Independent regression/limit tests. Python 3.10+, standard library only."""
from __future__ import annotations
import argparse
import csv
import math
from pathlib import Path
import subprocess
import tempfile


def read_csv(path):
    with Path(path).open(encoding="utf-8") as stream:
        return [{key: float(value) for key, value in row.items()} for row in csv.DictReader(stream)]


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--exe", default="./electrospinning_debug")
    args = parser.parse_args()
    exe = Path(args.exe).resolve()
    passed = []
    proc = subprocess.run([str(exe), "--self-test"], text=True, capture_output=True, check=True)
    print(proc.stdout.strip())
    passed.append("analytical planar field, manufactured Poisson refinement, charge deposition")
    with tempfile.TemporaryDirectory(prefix="spinning-tests-") as temp:
        root = Path(temp)

        def run(name, text, good=True):
            work = root / name
            work.mkdir()
            cfg = work / "case.nml"
            cfg.write_text("&spinning\nNr=24, Nz=32, Nj=101, stretch_r=5, stretch_z=3,\n" + text + "\n/\n")
            proc = subprocess.run([str(exe), "--config", str(cfg)], cwd=work,
                                  text=True, capture_output=True, timeout=180)
            if good:
                if proc.returncode:
                    raise AssertionError(f"{name} failed:\n{proc.stdout}\n{proc.stderr}")
                assert (work / "run_status.txt").read_text().startswith("NUMERICALLY_CONVERGED")
                jet = read_csv(work / "jet_profile.csv")
                with (work / "diagnostics.csv").open() as stream:
                    diagnostics = {row["metric"]: float(row["value"]) for row in csv.DictReader(stream)}
                assert all(math.isfinite(v) for row in jet for v in row.values())
                assert all(row["r_m"] > 0 and row["V_m_per_s"] > 0 for row in jet)
                return jet, diagnostics, work
            assert proc.returncode != 0, f"{name}: invalid run was accepted"
            assert not (work / "jet_profile.csv").exists()
            return proc, work

        quiet = "V_voltage=0, Jet_current_A=0, Conductivity=0, Gamma=0, g_z=0, Cd_drag=0,\n" \
                "h_conv=0, evaporation=.false., joule_heat=.false., viscous_heat=.false.,\n"
        jet, diag, _ = run("constant", quiet + "Gamma=0.055, space_charge=.false.,\n")
        u0, r0 = jet[0]["V_m_per_s"], jet[0]["r_m"]
        assert max(abs(row["V_m_per_s"] / u0 - 1) for row in jet) < 2e-10
        assert max(abs(row["r_m"] / r0 - 1) for row in jet) < 2e-10
        assert max(abs(row["T_C"] - 25) for row in jet) < 2e-9
        passed.append("uniform force-free jet; no artificial thinning or cooling")

        gravity_errors = []
        for nodes in (201, 401):
            text = quiet + f"g_z=9.81, rho_air=0, Viscosity=1e-10, Nj={nodes}, stretch_jet=0,\n" \
                "Q_pump=3.141592653589793e-10, jet_radius_in=1e-5, space_charge=.false.,\n"
            jet, _, _ = run(f"gravity{nodes}", text)
            # The final node carries the zero-stress outflow boundary: compare interior.
            row = jet[-2]
            exact = math.sqrt(1 + 2 * 9.81 * row["z_m"])
            gravity_errors.append(abs(row["V_m_per_s"] - exact))
        assert gravity_errors[1] < 0.6 * gravity_errors[0], gravity_errors
        assert gravity_errors[1] < 1e-3
        passed.append("inviscid gravity limit; first-order axial mesh convergence")

        jet, diag, _ = run("evaporation", quiet + "evaporation=.true., k_mass=0.5,\n")
        assert jet[-1]["solvent_mass_flow_kg_s"] < jet[0]["solvent_mass_flow_kg_s"]
        assert jet[-1]["T_C"] < jet[0]["T_C"]
        assert all(jet[i]["solvent_mass_flow_kg_s"] <= jet[i-1]["solvent_mass_flow_kg_s"]
                   for i in range(1, len(jet)))
        loss = sum(row["evaporated_segment_kg_s"] for row in jet)
        assert abs(loss - (jet[0]["solvent_mass_flow_kg_s"] - jet[-1]["solvent_mass_flow_kg_s"])) < 1e-23
        assert diag["thermal_global_relative_residual"] < 1e-8
        assert diag["volume_flow_relative_error"] < 1e-12
        passed.append("evaporation, latent cooling, cumulative solvent and energy balances")

        electric = "Jet_current_A=3e-9, stretch_r=11, Nr=64, Nz=96, Nj=201,\n"
        positive, dp, wp = run("positive", electric + "V_voltage=3500,\n")
        negative, dn, wn = run("negative", electric + "V_voltage=-3500,\n")
        for p, n in zip(positive, negative):
            for key in ("r_m", "V_m_per_s", "T_C"):
                assert abs(p[key] - n[key]) < 2e-10 * max(abs(p[key]), 1e-6)
            assert abs(p["lambda_C_m"] + n["lambda_C_m"]) < 1e-19
            assert abs(p["E_axial_V_m"] + n["E_axial_V_m"]) < 1e-3
        for key in ("current_relative_error", "volume_flow_relative_error", "polymer_relative_error",
                    "charge_deposition_relative_error"):
            assert dp[key] < 1e-12, (key, dp[key])
        assert dp["poisson_relative_residual"] < 5e-10
        assert dp["coupling_relative_error"] < 2e-5
        passed.append("voltage-polarity symmetry; finite coupled solution and conservation checks")

        _, low, _ = run("viscosity_low", electric + "V_voltage=3500, Viscosity=0.02,\n")
        _, high, _ = run("viscosity_high", electric + "V_voltage=3500, Viscosity=2.0,\n")
        assert abs(low["velocity_end_m_s"] - high["velocity_end_m_s"]) > 1e-3
        passed.append("viscosity changes the momentum solution")

        _, no_polarization, _ = run("no_polarization", electric + "V_voltage=3500, polarization_force=.false.,\n")
        assert abs(dp["velocity_end_m_s"]-no_polarization["velocity_end_m_s"]) > 1e-3
        passed.append("dielectric Maxwell force changes the momentum solution")

        proc, work = run("nonconverged", electric + "max_outer=1,\n", good=False)
        assert "FAILED" in (work / "run_status.txt").read_text()
        passed.append("iteration exhaustion reports failure and exports no new final profile")
        for name, text in (("bad_humidity", "Humidity=120,"), ("bad_viscosity", "Viscosity=0,"),
                           ("bad_fraction", "polymer_fraction=1,"), ("bad_nan", "Gamma=NaN,")):
            run(name, text, good=False)
        passed.append("invalid/non-finite input rejection")

        # Original format: check actual parser and optional defaults on an unforced case.
        old = root / "legacy"
        old.mkdir()
        (old / "config_spinning.txt").write_text(
            "0\n0.1\n0.000003\n0.055\n0.2\n1.38888888888889e-10\n40\n0\n1000\n25\n40\n")
        proc = subprocess.run([str(exe), "--preset", "quick"], cwd=old, text=True,
                              capture_output=True, timeout=180)
        if proc.returncode:
            raise AssertionError(f"Legacy config failed: {proc.stdout}\n{proc.stderr}")
        assert "legacy alpha_whip is ignored" in proc.stdout
        assert (old / "jet_profile.csv").exists()
        passed.append("original 11-line config compatibility")

    print("\nVALIDATION PASSED")
    for index, description in enumerate(passed, 1):
        print(f"{index}. {description}")


if __name__ == "__main__":
    main()
