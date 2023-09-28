# entLSM
numerical results on the entanglement LSM theorem

## aklt
A spin-1/2 chain is coupled to the spin-3/2 chain such that the total system is an AKLT.

`mpo_M.jl` generates the tensor that can be further assembled into an MPO, which is done in `gen_rho.jl` and `gen_rho_16.jl`, the latter for the case $L=16$. `ent_SA.jl` and `ent_SL.jl` calculate various entanglement entropies. This is used to obtain the conditional mutual information, plotted by `plot_mi.py`. `mpo_corr.jl`, `dmrg_corr.jl` and `plot_corr.py` calculate and plot the correlation functions. `ed_rho.jl` and `plot_spectrum.py` deal with the entanglement spectrum. Finally, `ed_twist.jl`, `dmrg_twist.jl` and `plot_twist.py` are used to calculate and plot the entanglement energy difference between twisted and untwisted states.

## ssb
A spin-1/2 Majumbdar–Ghosh chain is decohered by coupling to spin-3/2 modes.

`ed_H.jl` and `ed_K.jl` calculate the spectrum of the original and entanglement Hamiltonians, respectively. The results are plotted using `plot.py`.
