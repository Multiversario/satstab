# satstab
--------
Tools for determining the stability limit of a satellites and sub-satellites

This is a repository for tools to determine the stability limit, or critical semimajor axis (a_c), of exomoons and submoons. Boundaries (a_sat for exomoons and a_sub for submoons) are provided in terms of AUs. 

This summarized data in contour_moon.txt and contour_submoon.txt contains the summarized data resulting from N-bodt simulations performed with REBOUND (Rein 2012, Rein 2015). Explicitly, the files contain e_p, e_sat and a_crit (in units of Hill radius). 

The simulations for exomoons considered a Neptune-like exomoon orbiting a Jupiter-like planet, at a timescale of 10^5 yr, where orbital eccentricity of the planet and exomoon are varied between [0.0 - 0.5] in steps of 0.1. Orbits were established to be co-planar, with the argument of pericenter and ascending node set to zero. The planet is given an initial mean anomaly of 0 degrees, whilst for the exomoon’s 20 values of initial mean anomaly were randomly selected from a uniform distribution from 0 to 180 degrees. A similar procedure was used for submoons with a Neptune-like exomoon host (e.g., Kepler1625b-I). 

After cloning the repository, the tool for determining a_c for a given set of ecentricity parameters (e_p, e_sat) as well as the planet's semi-major axis (a_p) in AU, mass (m_p) and stellar mass (m_star) in solar masses. The stability limit can be determined simply by running **'python get_ac.py a_p m_p m_star e_p e_sat sub'**, where e_p and e_sat are floats between [0, 0.5] and sub is a flag to indicate whether to evaluate for an exomoon(sub=0) or submoon(sub=1). Additionally, a_p, m_p and m_star are positive floats.

#Attribution
--------
If you find this useful for your research, please cite this work using the information below. 

```
@article{RF2020,
author = {{Rosario-Franco}, M. and {Quarles}, B. and {Musielak}, Z. and {Cutz}, M.},
title = "{Stability Limits of Circumbinary Planets: Is There a Pile-up in the Kepler CBPs?}",
journal = {\aj},
archivePrefix = "arXiv",
notes = {Submitted}

}
```
