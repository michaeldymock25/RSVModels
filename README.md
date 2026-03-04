# RSVModels

This R package contains implementations for a range of respiratory syncytial virus (RSV) dynamic epidemiological models (hereafter known as transmission models) including for base, maternal vaccination and monoclonal antibody scenarios. Some (or all) of these models will be used in upcoming work evaluating the value of information for hypothetical trial designs.

These models and their implementation are largely based off the models developed in [Hogan et al. (2017)](https://doi.org/10.1016/j.vaccine.2017.09.043) ([implementation](https://github.com/abhogan/RSV_model)), [Giannini et al. (2024)](https://doi.org/10.1186/s12879-024-09400-2) ([implementation](https://github.com/fionagi/rsvmod)) and [Giannini et al. (2025)](https://doi.org/10.1016/j.vaccine.2025.127155) ([implementation](https://github.com/fionagi/rsvmod.imm)). Credit is due to the authors of the respective models.

## Structure

The repository is structured as follows:

### Data

Contains data required for modelling. Examples include *auspop.rds* and *mixing.csv* which contain the aggregate Australian population in 2021 and the computed social contact mixing matrix (see below), respectively.

### R

Contains a series of R scripts. These scripts include *parameters.R*, *models.R* and *functions.R* which contain the distributions of the transmission model parameters, the transmission models themselves, and accessory functions required for data manipulation (e.g., to generate the social contact mixing matrix), respectively. The *main.R* script loads the required packages, sources the required scripts and allows the user the run the models.

## References

1. Hogan AB, Campbell PT, Blyth CC et al. Potential impact of a maternal vaccine for RSV: A mathematical modelling study. Vaccine. 2017; 35 (45): 6172-9.
2. Giannini F, Hogan AB, Sarna M et al. Modelling respiratory syncytial virus age-specific risk of hospitalisation in term and preterm infants. BMC Infectious Diseases. 2024; 24: 510.
3. Giannini F, Hogan AB, Cameron E et al. Estimating the impact of Western Australia's first respiratory syncytial virus immunisation program for all infants: A mathematical modelling study. Vaccine. 2025; 56: 127155.
