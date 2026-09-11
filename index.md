---
layout: page
title: Start
menubar: menu
permalink: /index
---

*This document is being updated continuously*


# GAMERA

is a project launched by the [Max Planck Insitute for Nuclear Physics in Heidelberg (MPIK)](https://www.mpi-hd.mpg.de/mpi/en/),
an open-source C++/python package which handles the spectral modelling of non-thermally emitting astrophysical sources in a simple and modular way. It allows the user to devise time-dependent models of leptonic and hadronic particle populations in a general astrophysical context (including SNRs, PWNs and AGNs) and to compute their subsequent photon emission. 

GAMERA is written in C++ and can be wrapped to python as the **GA**MERA **P**ython **Pa**ckage (GAPPA)

The software is listed in the Astrophysics Source Code Library <a href="https://ascl.net/2203.007"><img src="https://img.shields.io/badge/ascl-2203.007-blue.svg?colorB=262255" alt="ascl:2203.007" /></a>



NEWS
====

**11-09-2026**

The documentation has now a new layout with a sidebar for better navigation through the individual pages and tutorials.


**07-01-2021**

After some time, we have added an example on how to use the python wapper of GAMERA in order to fit data. This example script can guide the user in fitting real astrophysical data using models produced by GAMERA.

Jump [here]({{ '/docs/fitting_data' | relative_url }}) if you want to quickly check it out.

**13-10-2020**
- implemented the gamma-gamma absorption effect (absorption only, no secondaries)
- implemented the nuclear enhancement factor for pi0 emission
- ensured internal consistency between pi0 and bremsstrahlung emission when using custom ambient compositions
- ionization losses made optional

Documentation & Tutorials
=========================
[Download & Installation]({{ '/docs/download_installation' | relative_url }})

[General Info]({{ '/docs/documentation' | relative_url }})

[Tutorials]({{ '/docs/tutorials_main' | relative_url }})










For questions and comments, feel free to write to `mbreuhaus@mpifr-bonn.mpg.de`


