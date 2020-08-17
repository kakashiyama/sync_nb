# sync_nb (ver.0)
  
Original ver. created by Kazumi Kashiyama on 17/08/20.

An implicit solver for synchrotron nebula emission. The code calculates the time evolution of the electron energy distribution in a pulsar wind nebula (PWN) and the synchrotron specrum for a simplified supernova dynamics. The electron injection spectrum is so far assumed to be a fixed-shape broken power law. Only synchrotron self absorption (SSA) is taken into account as the photon absorption process.


---
## calculating nebula spectra

`$ ./run`

---

## plotting the electron and photon spectra

`$ gnuplot plt/plt_e_dis.txt`

`$ gnuplot plt/plt_ph_spec.txt`

---

## input.dat
**list of input parameters**
1. Bp [G] : dipole magnetic field 
2. Prot0 [s] : rotation period at birth
3. Mej [Msun] : ejecta mass 
4. tmin [s] : start time
5. tmax [s] : finish time
6. fac_dt : time resolution
7. epsB : megnetic field amplification efficiency
8. epse : electron acceleration efficiency
9. gam_b : peak Lorentz factor of the electron injection
10. gam_max : maximum Lorentz factor 
11. p1 : low energy electron spectral index
12. p2 : high energy electron spectral index
13. Nbin_e : number of mesh for electron energy distribution
14. gam_ph_min [Me*c^2] : minimum photon energy
15. gam_ph_max [Me*c^2] : maximum photon energy
16. Nbin_ph : number of mesh for photon spectrum
17. output_file_num : number of output file

