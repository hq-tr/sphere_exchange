import numpy as np
import matplotlib.pyplot as plt
from argparse import ArgumentParser

ap = ArgumentParser()
ap.add_argument("--N_el", "-e", type=int, help="number of electrons")
ap.add_argument("--N_orb", "-o", type=int, help="number of orbitals of the **ground state**")
ap.add_argument("--v2", "-v", action="store_true", default=False,help="use this tag for 'v2' codes")
ap.add_argument("--line","-l", type=float, default=-999999.0, help="reference line")
ap.add_argument("--shift", "-s", type=float, default=0)

aa = ap.parse_args()

Ne = aa.N_el
No = aa.N_orb

if aa.v2:
	appen = "_v2"
else:
	appen = ""

with open(f"output_Lz{appen}/{Ne}e{No}_Lz.dat") as f:
	data = [list(map(float, x.split())) for x in f.readlines()]
	theta2 = np.array([x[0] for x in data])
	Lz2   = np.array([-x[1] for x in data])


with open(f"output_gap{appen}/{Ne}e{No}_gap.dat") as f:
	data = [list(map(float, x.split())) for x in f.readlines()]
	theta2 = np.array([x[0] for x in data])
	gap   = np.array([x[1] for x in data])


fig = plt.figure(figsize=(12,5))
ax1 = fig.add_subplot(121)
ax2 = fig.add_subplot(122)
Omega2  = (1-np.cos(theta2))/2 # This is actually Solid angle / 4π
#Lz2_mod = (Lz2 - Lz2[0])*(No-1)/(No+1)

m, b = np.polyfit(Omega2[-5:], Lz2[-5:],1) # fit y = mx + b for the last 4 points of data
print(f"m = {m}")
Lz2_mod = Lz2 - Lz2[0] - 2*Omega2*(Ne/2) # Here Nϕ = Nₒ
#Lz2_mod = Lz2 - Lz2[0] - Omega2*m

ax1.plot(Omega2, gap)
ax1.set_xlabel(r"$\Omega/4\pi$")
ax1.set_ylabel(r"$\Delta E$")

ax2.plot(Omega2, Lz2_mod, label="Braiding phase")
ax2.set_xlabel(r"$\Omega/4\pi$")
ax2.set_ylabel(r"$\langle L_z\rangle$")

if abs(aa.line)<2:
	ref = aa.line * np.ones(Omega2.size)
	ax2.plot(Omega2, ref, "k--", label=f"{aa.line}")
	plt.legend()

print(Lz2_mod)

print(gap)

plt.suptitle(f"N_e = {Ne}, N_o = {No}")

plt.savefig(f"plots/Lz_{Ne}e{No}.svg")
