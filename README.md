# satstab
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.3751872.svg)](https://doi.org/10.5281/zenodo.3751872)

Tools for determining the orbital stability limit of a satellites and sub-satellites

This is a repository with tools to determine the stability limit, or critical semimajor axis (a_c), of exomoons and submoons. Boundaries (a_sat for exomoons and a_sub for submoons) are provided in terms of astronomical units (AUs). 

The summarized data in contour_moon.txt and contour_submoon.txt contains the summarized data resulting from N-body simulations performed with REBOUND (Rein 2012, Rein 2015). Explicitly, the files contain e_p, e_sat and a_crit (in units of Hill radius). 

The simulations for exomoons considered a Neptune-like exomoon orbiting a Jupiter-like planet, at a timescale of 10^5 yr, where orbital eccentricity of the planet and exomoon are varied between [0.0 - 0.5] in steps of 0.1. Orbits were established to be co-planar, with the argument of pericenter and ascending node set to zero. The planet is given an initial mean anomaly of 0 degrees, whilst for the exomoon’s 20 values of initial mean anomaly were randomly selected from a uniform distribution from 0 to 180 degrees. A similar procedure was used for submoons with a Neptune-like exomoon host (e.g., Kepler1625b-I). 

After cloning the repository, the tool for determining a_c for a given set of ecentricity parameters (e_p, e_sat) as well as the planet's semi-major axis (a_p) in AU, mass (m_p) and stellar mass (m_star) in solar masses. The stability limit can be determined simply by running **'python get_ac.py a_p m_p m_star e_p e_sat sub'**, where e_p and e_sat are floats between [0, 0.5] and sub is a flag to indicate whether to evaluate for an exomoon(sub=0) or submoon(sub=1). Additionally, a_p, m_p and m_star are positive floats.

# Attribution
--------
If you find this useful for your research, please cite this work using the information below. 

```
@ARTICLE{2020AJ....159..260R,
       author = {{Rosario-Franco}, Marialis and {Quarles}, Billy and
         {Musielak}, Zdzislaw E. and {Cuntz}, Manfred},
        title = "{Orbital Stability of Exomoons and Submoons with Applications to Kepler 1625b-I}",
      journal = {\aj},
     keywords = {Orbits, Exoplanet dynamics, Exoplanets, Natural satellites (Extrasolar), Astrophysics - Earth and Planetary Astrophysics},
         year = 2020,
        month = jun,
       volume = {159},
       number = {6},
          eid = {260},
        pages = {260},
          doi = {10.3847/1538-3881/ab89a7},
archivePrefix = {arXiv},
       eprint = {2005.06521},
 primaryClass = {astro-ph.EP},
       adsurl = {https://ui.adsabs.harvard.edu/abs/2020AJ....159..260R},
      adsnote = {Provided by the SAO/NASA Astrophysics Data System}
}

```
