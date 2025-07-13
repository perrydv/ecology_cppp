
Berkeley ESPM QuantEco lab group’s study on calibrated posterior
predictive p-values in ecology.

[Link to meeting
notes](https://docs.google.com/document/d/1ZPJUCIgU28_Fm_4gV6fgKoFihSasZjHoTyGMsj6vVMU/edit?tab=t.0#heading=h.v38mx9tc4ayu)

## File structure

    ecology_cppp
    │   README.md
    │   utils.R # general auxiliary functions relevant for all model types    
    │
    └───cppp
    │   │
    │   └─────newcomb_example # Sally's toy example
    │   │     │   model.R
    │   │     │   light.txt
    │   │
    │   └─────sally_code # code from Paganin et al.
    │         │   1_runCPPP.R
    │         │   calculateCPPP_original.R
    │         │   registeredDiscrepancies.R
    │   
    └───figures
    │   │
    │   └─────occupancy # occupancy figures
    │
    └───occupancy
    │   │
    │   │  calculate_CPPP.R # runs simulations to calculate cppp, ppp, coverage; generates plots
    │   │  minimal_model.R # occupancy model specifications
    │   │  occupancy_discrepancy_functions.R # occupancy model-specific discrepancy functions
    │   │  occupancy_model.R # occupancy model specifications that calculate posterior predictive samples in model code (let's phase this out)
    │   │  simulate_data.R # simulates detection/non-detection data for different occu model specifications

## Coding style

We adhere to Tidystyle, except for a few exceptions outlined in
`.lintr`, notably that we allow for camel case to ensure consistency
with `nimble` style.

To check that code adheres to our style, use `lintr`:

``` r
lintr::lint_dir()
```
