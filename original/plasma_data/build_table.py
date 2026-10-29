"""Convert the pinned BOLSIG+/Phelps air tables; standard-library Python only.

The original dataset is from MD-CWI/afivo-streamer, blob
e3c3027260ed11716fe8e377b8070cce72cc99fa. Artificial dry air: 80% N2 / 20% O2.
Reaction C27 has a mandatory 1e-6 unit conversion specified in the source.
"""
from pathlib import Path
import csv
import hashlib


def table(lines, name):
    start = lines.index(name) + 1
    factor = 1.0
    while not lines[start].startswith("---"):
        if lines[start].startswith("FACTOR:"):
            factor = float(lines[start].split(":", 1)[1])
        start += 1
    result = []
    for line in lines[start + 1:]:
        if line.startswith("---"):
            break
        if line.strip():
            x, y = map(float, line.split())
            result.append((x, factor*y))
    return result


def main():
    root = Path(__file__).resolve().parent
    source = root / "air_light_example_v0.txt"
    data = source.read_bytes()
    git_sha = hashlib.sha1(f"blob {len(data)}\0".encode()+data).hexdigest()
    assert git_sha == "e3c3027260ed11716fe8e377b8070cce72cc99fa", git_sha
    lines = data.decode().splitlines()
    names = ["Mean energy (eV)", "Mobility *N (1/m/V/s)", "Diffusion coefficient *N (1/m/s)",
             "C25 N2 Ionization 15.60 eV", "C26 N2 Ionization 18.80 eV",
             "C43 O2 Ionization 12.06 eV", "C28 O2 Attachment", "C27 O2 Attachment"]
    tables = [table(lines, n) for n in names]
    assert all([r[0] for r in t] == [r[0] for r in tables[0]] for t in tables)
    target = root.parent / "air_transport.csv"
    with target.open("w", newline="") as stream:
        writer = csv.writer(stream)
        writer.writerow(["EN_Td", "mean_energy_eV", "muN_1_m_V_s", "DN_1_m_s",
                         "ki_m3_s", "ka2_m3_s", "ka3_m6_s"])
        for rows in zip(*tables):
            en = rows[0][0]
            e, mu, diff, i1, i2, io, a2, a3 = [r[1] for r in rows]
            values = [en, e, mu, diff, .8*(i1+i2)+.2*io, .2*a2, .2**2*a3]
            writer.writerow([f"{x:.16e}" for x in values])
    print(target)


if __name__ == "__main__":
    main()
