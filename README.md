Berkeley ESPM QuantEco lab group's study on calibrated posterior predictive p-values in ecology.

`simulate_data.R`: functions for simulating data under different model variations

`occupancy_model.R`: nimble models and functions for calculating discrepancy measures and running models

- `_latent` extension includes models that do not condition the generation of $y^{rep}$ on the latent state, $z$

`main.R`: run MCMC and calculate p values for model variations

- `_latent` extension includes running models that do not condition the generation of $y^{rep}$ on the latent state, $z$

[Link to meeting notes](https://docs.google.com/document/d/1ZPJUCIgU28_Fm_4gV6fgKoFihSasZjHoTyGMsj6vVMU/edit?tab=t.0#heading=h.v38mx9tc4ayu)