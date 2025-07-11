Berkeley ESPM QuantEco lab group's study on calibrated posterior predictive p-values in ecology.

`simulate_data.R`: functions for simulating data under different model variations

`occupancy_model.R`: nimble models and functions for calculating discrepancy measures and running models

`main_null.R`: run MCMC and calculate p values for model variations under the *null* hypothesis that the model generated the data.

`main_alt.R`: run MCMC and calculate p values for model variations under the *alternative* hypothesis that the model did not generate the data.

[Link to meeting notes](https://docs.google.com/document/d/1ZPJUCIgU28_Fm_4gV6fgKoFihSasZjHoTyGMsj6vVMU/edit?tab=t.0#heading=h.v38mx9tc4ayu)

I followed an example in another repository and wrapped the directory structure within a pair of triple backticks (```):

```
ecology_cppp
│   README.md # readme
│   file001.txt    
│
└───folder1
│   │   file011.txt
│   │   file012.txt
│   │
│   └───subfolder1
│       │   file111.txt
│       │   file112.txt
│       │   ...
│   
└───folder2
    │   file021.txt
    │   file022.txt
```