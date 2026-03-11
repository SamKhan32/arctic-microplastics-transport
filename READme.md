# arctic-microplastics-transport

Undergraduate research project modeling microplastics and PFAS transport in the Arctic Ocean using [MITgcm](https://mitgcm.org/). This repository collects all experiments, from an initial barotropic baseline through progressively more complete configurations.

## Background

This project is part of Dr. Roger Tipton's NSF-sponsored research on PFAS chemical distribution in the Arctic Ocean, tracing contamination patterns between Alaska and Norway to establish baseline data before increased Arctic shipping opens the region to broader pollution pathways.

I joined this work as an Office of Undergraduate Research scholar at UNC Charlotte, paired with Dr. Tipton's lab group. The modeling work here — grid setup, ERA5 wind forcing, and passive tracer transport — is my contribution to that broader research effort.

## What this is

A wind-forced barotropic Arctic circulation model built as a baseline for passive tracer advection. The goal is to simulate how surface contaminants move under realistic wind stress, starting from a single-layer barotropic configuration and building toward a more complete representation of Arctic surface circulation.

The immediate output is a velocity field over the Arctic cap that can advect a passive tracer — used here as a proxy for surface microplastics and PFAS — to produce transport visualizations.

## Setup

- **Grid:** ECCO LLC270 Arctic cap (tile 7, ~67.5°N–90°N), extracted from the ECCO Version 5 nctiles grid
- **Forcing:** ERA5 monthly wind stress climatology (1991–2020), interpolated to the LLC270 grid
- **Bathymetry:** IBCAO-informed, via ECCO Depth field (negated, land set to 0)
- **Model:** MITgcm, barotropic (Nr=1), curvilinear configuration
- **HPC:** Orion cluster, SLURM, GCC + OpenMPI

## Status

Work in progress. Presented at NCUR 2026.

## Acknowledgments

This work was supported by the UNC Charlotte Office of Undergraduate Research. NSF-sponsored research group led by Dr. Roger Tipton, UNC Charlotte.

ECCO grid data from [ecco.jpl.nasa.gov](https://ecco.jpl.nasa.gov). ERA5 forcing data from ECMWF via the Copernicus Climate Data Store.