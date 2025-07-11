Berkeley ESPM QuantEco lab group's study on calibrated posterior predictive p-values in ecology.

`simulate_data.R`: functions for simulating data under different model variations

`occupancy_model.R`: nimble models and functions for calculating discrepancy measures and running models

`main_null.R`: run MCMC and calculate p values for model variations under the *null* hypothesis that the model generated the data.

`main_alt.R`: run MCMC and calculate p values for model variations under the *alternative* hypothesis that the model did not generate the data.

[Link to meeting notes](https://docs.google.com/document/d/1ZPJUCIgU28_Fm_4gV6fgKoFihSasZjHoTyGMsj6vVMU/edit?tab=t.0#heading=h.v38mx9tc4ayu)

├── data/ # Raw and processed data files
│ ├── raw/ # Unmodified original data
│ └── processed/ # Cleaned and transformed data
├── scripts/ # R or Python scripts for analysis
│ ├── analysis.R # Main analysis script
│ └── utils.R # Helper functions
├── results/ # Output files: figures, tables, etc.
│ ├── figures/ # Plots and visualizations
│ └── tables/ # Summary tables
├── docs/ # Documentation, markdowns, or rendered reports
├── app/ # Code for Shiny or web application
├── .gitignore # Git ignore rules
├── README.md # Project overview and structure
└── requirements.txt # Python dependencies (or use renv.lock for R)