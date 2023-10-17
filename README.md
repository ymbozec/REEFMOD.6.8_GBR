# REEFMOD.6.8_GBR

This repository contains the scripts of ReefMod-GBR, a coral individual-based model that simulates coral populations across 2,300 km of Australia's Great Barrier Reef (GBR).

The first version (Feb-Apr 2022) was used to forecast possible futures of the GBR under climate change scenarios (CMIP-5) and assess the possible benefits of reef restoration techniques for the Reef Restoration and Adaptation Program (RRAP: https://gbrrestoration.org/).

The last update (March 2023) integrates CMIP-6 climate change scenarios downscaled at 10-km resolution provided by McWhorter et al. (2022): the projection of heat stress (maximum annual Degree Heating Weeks) from five global circulation models (CNRM-ESM2-1, EC-Earth3-Veg, IPSL-CM6A-LR, MRI-ESM2-0, UKESM1-0-LL) under five carbon emission scenarios (SSP1-119, SSP1-2.6, SSP2-4.5, SSP3-7.0, SSP5-8.5). This update also integrates the revised area estimates of coral reef habitat for 3,806 individual reefs across the GBR based on the 3D surface area of geomorphic classes (Roelfsema et al. 2021) extended to inshore reefs (Castro-Sanguino et al. 2023).

Either CMIP-5 or CMIP-6 models can be selected from the front script MAIN_REEFMOD_GBR.m. The script is designed to run the model under one climate change scenario (ie, one CMIP-5 or CMIP-6 model with the available scenario of carbon emission RCP/SSP) chosen by the user. The number of repeat simulations, ie, the stochastic simulation of the same warming scenario under different projections of cyclones and other random events (including the initialisation of coral cover and CoTS density, the magnitude of coral mortality, the forcing scheme of water quality,...) can be set with 'NB_SIMULATIONS' (set to 20 in RRAP).

The model reconstructs coral trajectories across the GBR between 2008-2022 and forecasts possible coral futures (2022-2100) based on temporally- and spatially-explicit forcing of water quality, cyclones, heat stress (coral bleaching) and the simulated population dynamics of the coral-eating (crown-of-thorns) starfish. Options for simulating management interventions on any given reef include:
- the outplanting of corals of specified species group, size and heat tolerance at a specified density
- the enrichment of coral larvae of specified species group with a specified density
- the reduction of heat stress through Solar Radiation Management (fogging)
- the consolidation of lose coral rubble to increase survival of coral recruits

This version also simulates the CoTS control program in space and time from 2019 onwards, with a specific number of boats (default: 5 boats) and the list of priority reefs as currently (2021) defined by GBRMPA.

Yves-Marie Bozec, The University of Queensland (y.bozec@uq.edu.au)

Citation: Bozec, Y.-M., K. Hock, R. A. Mason, M. E. Baird, C. Castro-Sanguino, S. A. Condie, M. Puotinen, A. Thompson, and P. J. Mumby. 2022. Cumulative impacts across Australia’s Great Barrier Reef: A mechanistic evaluation. Ecological Monographs 92(1), e01494
https://esajournals.onlinelibrary.wiley.com/doi/abs/10.1002/ecm.1494
