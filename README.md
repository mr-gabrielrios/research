## SWISHE: suppressed wind-induced surface heat exchange

#### Data availability

A frozen data state is in effect as of 2024-03-07. Run data is available for HIRAM, AM2.5, and FLOR. These runs were performed in the fall of 2023 but were not initialized properly, either with improper diagnostic tables (AM2.5, HIRAM) or initial conditions. 

The shortcomings are minor for AM2.5 and HIRAM, since I’m just missing some diagnostics that I might be able to derive from the data available. The shortcomings are larger for FLOR, since the initial conditions used do not represent a climate system in quasi-equilibrium. 

To address this, new runs are in process for all models, with FLOR and HIRAM taking priority due to their ability to resolve the coupled system and TCs in high-fidelity, respectively. However, due to computing constraints on Tiger, the new runs risk being unavailable for sufficient analysis in time for Generals.

To address this, the old runs for AM2.5, FLOR, and HIRAM are being used. Data is available from 0101-0150 for AM2.5 and HiRAM and 2001-2140 for FLOR. Note that AM2.5 and HIRAM are bitwise reproducible by Wenchang Yang, FLOR is not.

| Model | Configuration | Years     | Location                                                                                      |
|-------|---------------|-----------|-----------------------------------------------------------------------------------------------|
| AM2.5 | Control       | 0101-0150 | /projects/GEOCLIM/gr7610/MODEL_OUT/AM2.5/v1/CTL1990s_tigercpu_intelmpi_18_540PE               |
| AM2.5 | SWISHE        | 0101-0150 | /projects/GEOCLIM/gr7610/MODEL_OUT/AM2.5/v1/CTL1990s_swishe_tigercpu_intelmpi_18_540PE        |
| HIRAM | Control       | 0101-0150 | /projects/GEOCLIM/gr7610/MODEL_OUT/HIRAM/v1/CTL1990s_tigercpu_intelmpi_18_540PE               |
| HIRAM | SWISHE        | 0101-0150 | /projects/GEOCLIM/gr7610/MODEL_OUT/HIRAM/v1/CTL1990s_swishe_tigercpu_intelmpi_18_540PE        |
| FLOR  | Control       | 2001-2125 | /projects/GEOCLIM/gr7610/MODEL_OUT/FLOR/v1/CTL1990s_v201905_tigercpu_intelmpi_18_576PE        |
| FLOR  | SWISHE        | 2001-2125 | /projects/GEOCLIM/gr7610/MODEL_OUT/FLOR/v1/CTL1990s_v201905_swishe_tigercpu_intelmpi_18_576PE |
