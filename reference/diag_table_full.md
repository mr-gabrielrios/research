|Module|Field|Long Name|Units|Number of Axis|Time Axis|Missing Value|Min Value|Max Value|AXES LIST|
|------|-----|---------|-----|--------------|---------|-------------|---------|---------|--------|
|dynamics|bk|vertical coordinate sigma value|none|1|F||||phalf|
|dynamics|pk|pressure part of the hybrid coordinate|pascal|1|F||||phalf|
|dynamics|hyam|vertical coordinate A value|1E-5 Pa|1|F||||pfull|
|dynamics|hybm|vertical coordinate B value|none|1|F||||pfull|
|dynamics|grid_lon|longitude|degrees_E|2|F||||grid_x,grid_y|
|dynamics|grid_lat|latitude|degrees_N|2|F||||grid_x,grid_y|
|dynamics|grid_lont|longitude|degrees_E|2|F||||grid_xt,grid_yt|
|dynamics|grid_latt|latitude|degrees_N|2|F||||grid_xt,grid_yt|
|dynamics|area|cell area|m**2|2|F||||grid_xt,grid_yt|
|dynamics|dx|dx|m|2|F||||grid_xt,grid_y|
|dynamics|dy|dy|m|2|F||||grid_x,grid_yt|
|dynamics|zsurf|surface height|m|2|F||||grid_xt,grid_yt|
|dynamics|zs|Original Mean Terrain|m|2|F||||grid_xt,grid_yt|
|dynamics|ze|Hybrid_Z_surface|m|3|F||||grid_xt,grid_yt,pfull|
|dynamics|oro|Land/Water Mask|none|2|F||||grid_xt,grid_yt|
|dynamics|sgh|Terrain Standard deviation|m|2|F||||grid_xt,grid_yt|
|dynamics|ps_ic|initial surface pressure|Pa|2|F||||grid_xt,grid_yt|
|dynamics|ua_ic|initial zonal wind|m/s|3|F||||grid_xt,grid_yt,pfull|
|dynamics|va_ic|initial meridional wind|m/s|3|F||||grid_xt,grid_yt,pfull|
|dynamics|ppt_ic|initial potential temperature|K|3|F||||grid_xt,grid_yt,pfull|
|dynamics|sphum_ic|initial surface pressure|Pa|2|F||||grid_xt,grid_yt|
|dynamics|ps|surface pressure|Pa|2|T| -1.0000000E+10|   40000.00|   110000.0|grid_xt,grid_yt|
|dynamics|mq|mountain torque|Hadleys per unit area|2|T| -1.0000000E+10|||grid_xt,grid_yt|
|dynamics|aam|angular momentum|kg*m^2/s|2|T| -1.0000000E+10|||grid_xt,grid_yt|
|dynamics|amdt|angular momentum error|kg*m^2/s^2|2|T| -1.0000000E+10|||grid_xt,grid_yt|
|dynamics|pret|total precipitation|mm/day|2|T| -1.0000000E+10|||grid_xt,grid_yt|
|dynamics|prew|water precipitation|mm/day|2|T| -1.0000000E+10|||grid_xt,grid_yt|
|dynamics|prer|rain precipitation|mm/day|2|T| -1.0000000E+10|||grid_xt,grid_yt|
|dynamics|prei|ice precipitation|mm/day|2|T| -1.0000000E+10|||grid_xt,grid_yt|
|dynamics|pres|snow precipitation|mm/day|2|T| -1.0000000E+10|||grid_xt,grid_yt|
|dynamics|preg|graupel precipitation|mm/day|2|T| -1.0000000E+10|||grid_xt,grid_yt|
|dynamics|qv_dt_gfdlmp|water vapor specific humidity tendency from GFDL MP|kg/kg/s|3|T| -1.0000000E+10|||grid_xt,grid_yt,pfull|
|dynamics|ql_dt_gfdlmp|total liquid water tendency from GFDL MP|kg/kg/s|3|T| -1.0000000E+10|||grid_xt,grid_yt,pfull|
|dynamics|qi_dt_gfdlmp|total ice water tendency from GFDL MP|kg/kg/s|3|T| -1.0000000E+10|||grid_xt,grid_yt,pfull|
|dynamics|liq_wat_dt_gfdlmp|liquid water tracer tendency from GFDL MP|kg/kg/s|3|T| -1.0000000E+10|||grid_xt,grid_yt,pfull|
|dynamics|ice_wat_dt_gfdlmp|ice water tracer tendency from GFDL MP|kg/kg/s|3|T| -1.0000000E+10|||grid_xt,grid_yt,pfull|
|dynamics|qr_dt_gfdlmp|rain water tendency from GFDL MP|kg/kg/s|3|T| -1.0000000E+10|||grid_xt,grid_yt,pfull|
|dynamics|qg_dt_gfdlmp|graupel tendency from GFDL MP|kg/kg/s|3|T| -1.0000000E+10|||grid_xt,grid_yt,pfull|
|dynamics|qs_dt_gfdlmp|snow water tendency from GFDL MP|kg/kg/s|3|T| -1.0000000E+10|||grid_xt,grid_yt,pfull|
|dynamics|T_dt_gfdlmp|temperature tendency from GFDL MP|K/s|3|T| -1.0000000E+10|||grid_xt,grid_yt,pfull|
|dynamics|u_dt_gfdlmp|zonal wind tendency from GFDL MP|m/s/s|3|T| -1.0000000E+10|||grid_xt,grid_yt,pfull|
|dynamics|v_dt_gfdlmp|meridional wind tendency from GFDL MP|m/s/s|3|T| -1.0000000E+10|||grid_xt,grid_yt,pfull|
|dynamics|T_dt_phys|temperature tendency from physics|K/s|3|T| -1.0000000E+10|||grid_xt,grid_yt,pfull|
|dynamics|u_dt_phys|zonal wind tendency from physics|m/s/s|3|T| -1.0000000E+10|||grid_xt,grid_yt,pfull|
|dynamics|v_dt_phys|meridional wind tendency from physics|m/s/s|3|T| -1.0000000E+10|||grid_xt,grid_yt,pfull|
|dynamics|qv_dt_phys|water vapor specific humidity tendency from physics|kg/kg/s|3|T| -1.0000000E+10|||grid_xt,grid_yt,pfull|
|dynamics|ql_dt_phys|total liquid water tendency from physics|kg/kg/s|3|T| -1.0000000E+10|||grid_xt,grid_yt,pfull|
|dynamics|qi_dt_phys|total ice water tendency from physics|kg/kg/s|3|T| -1.0000000E+10|||grid_xt,grid_yt,pfull|
|dynamics|liq_wat_dt_phys|liquid water tracer tendency from physics|kg/kg/s|3|T| -1.0000000E+10|||grid_xt,grid_yt,pfull|
|dynamics|ice_wat_dt_phys|ice water tracer tendency from physics|kg/kg/s|3|T| -1.0000000E+10|||grid_xt,grid_yt,pfull|
|dynamics|qr_dt_phys|rain water tendency from physics|kg/kg/s|3|T| -1.0000000E+10|||grid_xt,grid_yt,pfull|
|dynamics|qg_dt_phys|graupel tendency from physics|kg/kg/s|3|T| -1.0000000E+10|||grid_xt,grid_yt,pfull|
|dynamics|qs_dt_phys|snow water tendency from physics|kg/kg/s|3|T| -1.0000000E+10|||grid_xt,grid_yt,pfull|
|dynamics|T_dt_sg|temperature tendency from 2dz subgrid mixing|K/s|3|T| -1.0000000E+10|||grid_xt,grid_yt,pfull|
|dynamics|u_dt_sg|zonal wind tendency from 2dz subgrid mixing|m/s/s|3|T| -1.0000000E+10|||grid_xt,grid_yt,pfull|
|dynamics|v_dt_sg|meridional wind tendency from 2dz subgrid mixing|m/s/s|3|T| -1.0000000E+10|||grid_xt,grid_yt,pfull|
|dynamics|qv_dt_sg|water vapor tendency from 2dz subgrid mixing|kg/kg/s|3|T| -1.0000000E+10|||grid_xt,grid_yt,pfull|
|dynamics|t_dt_nudge|temperature tendency from nudging|K/s|3|T| -1.0000000E+10|||grid_xt,grid_yt,pfull|
|dynamics|ps_dt_nudge|surface pressure tendency from nudging|Pa/s|2|T| -1.0000000E+10|||grid_xt,grid_yt|
|dynamics|delp_dt_nudge|pressure thickness tendency from nudging|Pa/s|3|T| -1.0000000E+10|||grid_xt,grid_yt,pfull|
|dynamics|u_dt_nudge|zonal wind tendency from nudging|m/s/s|3|T| -1.0000000E+10|||grid_xt,grid_yt,pfull|
|dynamics|v_dt_nudge|meridional wind tendency from nudging|m/s/s|3|T| -1.0000000E+10|||grid_xt,grid_yt,pfull|
|dynamics|qv_dt_nudge|specific humidity tendency from nudging|kg/kg/s|3|T| -1.0000000E+10|||grid_xt,grid_yt,pfull|
|dynamics|z50|50-mb height|m|2|T| -1.0000000E+10|||grid_xt,grid_yt|
|dynamics|u50|50-mb u|m/s|2|T| -1.0000000E+10|||grid_xt,grid_yt|
|dynamics|v50|50-mb v|m/s|2|T| -1.0000000E+10|||grid_xt,grid_yt|
|dynamics|t50|50-mb temperature|K|2|T| -1.0000000E+10|||grid_xt,grid_yt|
|dynamics|q50|50-mb specific humidity|kg/kg|2|T| -1.0000000E+10|||grid_xt,grid_yt|
|dynamics|ql50|50-mb cloud water mass mixing ratio|kg/kg|2|T| -1.0000000E+10|||grid_xt,grid_yt|
|dynamics|qi50|50-mb cloud ice mass mixing ratio|kg/kg|2|T| -1.0000000E+10|||grid_xt,grid_yt|
|dynamics|qr50|50-mb rain mass mixing ratio|kg/kg|2|T| -1.0000000E+10|||grid_xt,grid_yt|
|dynamics|qs50|50-mb snow mass mixing ratio|kg/kg|2|T| -1.0000000E+10|||grid_xt,grid_yt|
|dynamics|qg50|50-mb graupel mass mixing ratio|kg/kg|2|T| -1.0000000E+10|||grid_xt,grid_yt|
|dynamics|cf50|50-mb cloud fraction|1|2|T| -1.0000000E+10|||grid_xt,grid_yt|
|dynamics|omg50|50-mb omega|Pa/s|2|T| -1.0000000E+10|||grid_xt,grid_yt|
|dynamics|z200|200-mb height|m|2|T| -1.0000000E+10|||grid_xt,grid_yt|
|dynamics|u200|200-mb u|m/s|2|T| -1.0000000E+10|||grid_xt,grid_yt|
|dynamics|v200|200-mb v|m/s|2|T| -1.0000000E+10|||grid_xt,grid_yt|
|dynamics|t200|200-mb temperature|K|2|T| -1.0000000E+10|||grid_xt,grid_yt|
|dynamics|q200|200-mb specific humidity|kg/kg|2|T| -1.0000000E+10|||grid_xt,grid_yt|
|dynamics|ql200|200-mb cloud water mass mixing ratio|kg/kg|2|T| -1.0000000E+10|||grid_xt,grid_yt|
|dynamics|qi200|200-mb cloud ice mass mixing ratio|kg/kg|2|T| -1.0000000E+10|||grid_xt,grid_yt|
|dynamics|qr200|200-mb rain mass mixing ratio|kg/kg|2|T| -1.0000000E+10|||grid_xt,grid_yt|
|dynamics|qs200|200-mb snow mass mixing ratio|kg/kg|2|T| -1.0000000E+10|||grid_xt,grid_yt|
|dynamics|qg200|200-mb graupel mass mixing ratio|kg/kg|2|T| -1.0000000E+10|||grid_xt,grid_yt|
|dynamics|cf200|200-mb cloud fraction|1|2|T| -1.0000000E+10|||grid_xt,grid_yt|
|dynamics|omg200|200-mb omega|Pa/s|2|T| -1.0000000E+10|||grid_xt,grid_yt|
|dynamics|z500|500-mb height|m|2|T| -1.0000000E+10|||grid_xt,grid_yt|
|dynamics|u500|500-mb u|m/s|2|T| -1.0000000E+10|||grid_xt,grid_yt|
|dynamics|v500|500-mb v|m/s|2|T| -1.0000000E+10|||grid_xt,grid_yt|
|dynamics|t500|500-mb temperature|K|2|T| -1.0000000E+10|||grid_xt,grid_yt|
|dynamics|q500|500-mb specific humidity|kg/kg|2|T| -1.0000000E+10|||grid_xt,grid_yt|
|dynamics|ql500|500-mb cloud water mass mixing ratio|kg/kg|2|T| -1.0000000E+10|||grid_xt,grid_yt|
|dynamics|qi500|500-mb cloud ice mass mixing ratio|kg/kg|2|T| -1.0000000E+10|||grid_xt,grid_yt|
|dynamics|qr500|500-mb rain mass mixing ratio|kg/kg|2|T| -1.0000000E+10|||grid_xt,grid_yt|
|dynamics|qs500|500-mb snow mass mixing ratio|kg/kg|2|T| -1.0000000E+10|||grid_xt,grid_yt|
|dynamics|qg500|500-mb graupel mass mixing ratio|kg/kg|2|T| -1.0000000E+10|||grid_xt,grid_yt|
|dynamics|cf500|500-mb cloud fraction|1|2|T| -1.0000000E+10|||grid_xt,grid_yt|
|dynamics|omg500|500-mb omega|Pa/s|2|T| -1.0000000E+10|||grid_xt,grid_yt|
|dynamics|z700|700-mb height|m|2|T| -1.0000000E+10|||grid_xt,grid_yt|
|dynamics|u700|700-mb u|m/s|2|T| -1.0000000E+10|||grid_xt,grid_yt|
|dynamics|v700|700-mb v|m/s|2|T| -1.0000000E+10|||grid_xt,grid_yt|
|dynamics|t700|700-mb temperature|K|2|T| -1.0000000E+10|||grid_xt,grid_yt|
|dynamics|q700|700-mb specific humidity|kg/kg|2|T| -1.0000000E+10|||grid_xt,grid_yt|
|dynamics|ql700|700-mb cloud water mass mixing ratio|kg/kg|2|T| -1.0000000E+10|||grid_xt,grid_yt|
|dynamics|qi700|700-mb cloud ice mass mixing ratio|kg/kg|2|T| -1.0000000E+10|||grid_xt,grid_yt|
|dynamics|qr700|700-mb rain mass mixing ratio|kg/kg|2|T| -1.0000000E+10|||grid_xt,grid_yt|
|dynamics|qs700|700-mb snow mass mixing ratio|kg/kg|2|T| -1.0000000E+10|||grid_xt,grid_yt|
|dynamics|qg700|700-mb graupel mass mixing ratio|kg/kg|2|T| -1.0000000E+10|||grid_xt,grid_yt|
|dynamics|cf700|700-mb cloud fraction|1|2|T| -1.0000000E+10|||grid_xt,grid_yt|
|dynamics|omg700|700-mb omega|Pa/s|2|T| -1.0000000E+10|||grid_xt,grid_yt|
|dynamics|z850|850-mb height|m|2|T| -1.0000000E+10|||grid_xt,grid_yt|
|dynamics|u850|850-mb u|m/s|2|T| -1.0000000E+10|||grid_xt,grid_yt|
|dynamics|v850|850-mb v|m/s|2|T| -1.0000000E+10|||grid_xt,grid_yt|
|dynamics|t850|850-mb temperature|K|2|T| -1.0000000E+10|||grid_xt,grid_yt|
|dynamics|q850|850-mb specific humidity|kg/kg|2|T| -1.0000000E+10|||grid_xt,grid_yt|
|dynamics|ql850|850-mb cloud water mass mixing ratio|kg/kg|2|T| -1.0000000E+10|||grid_xt,grid_yt|
|dynamics|qi850|850-mb cloud ice mass mixing ratio|kg/kg|2|T| -1.0000000E+10|||grid_xt,grid_yt|
|dynamics|qr850|850-mb rain mass mixing ratio|kg/kg|2|T| -1.0000000E+10|||grid_xt,grid_yt|
|dynamics|qs850|850-mb snow mass mixing ratio|kg/kg|2|T| -1.0000000E+10|||grid_xt,grid_yt|
|dynamics|qg850|850-mb graupel mass mixing ratio|kg/kg|2|T| -1.0000000E+10|||grid_xt,grid_yt|
|dynamics|cf850|850-mb cloud fraction|1|2|T| -1.0000000E+10|||grid_xt,grid_yt|
|dynamics|omg850|850-mb omega|Pa/s|2|T| -1.0000000E+10|||grid_xt,grid_yt|
|dynamics|z925|925-mb height|m|2|T| -1.0000000E+10|||grid_xt,grid_yt|
|dynamics|u925|925-mb u|m/s|2|T| -1.0000000E+10|||grid_xt,grid_yt|
|dynamics|v925|925-mb v|m/s|2|T| -1.0000000E+10|||grid_xt,grid_yt|
|dynamics|t925|925-mb temperature|K|2|T| -1.0000000E+10|||grid_xt,grid_yt|
|dynamics|q925|925-mb specific humidity|kg/kg|2|T| -1.0000000E+10|||grid_xt,grid_yt|
|dynamics|ql925|925-mb cloud water mass mixing ratio|kg/kg|2|T| -1.0000000E+10|||grid_xt,grid_yt|
|dynamics|qi925|925-mb cloud ice mass mixing ratio|kg/kg|2|T| -1.0000000E+10|||grid_xt,grid_yt|
|dynamics|qr925|925-mb rain mass mixing ratio|kg/kg|2|T| -1.0000000E+10|||grid_xt,grid_yt|
|dynamics|qs925|925-mb snow mass mixing ratio|kg/kg|2|T| -1.0000000E+10|||grid_xt,grid_yt|
|dynamics|qg925|925-mb graupel mass mixing ratio|kg/kg|2|T| -1.0000000E+10|||grid_xt,grid_yt|
|dynamics|cf925|925-mb cloud fraction|1|2|T| -1.0000000E+10|||grid_xt,grid_yt|
|dynamics|omg925|925-mb omega|Pa/s|2|T| -1.0000000E+10|||grid_xt,grid_yt|
|dynamics|z1000|1000-mb height|m|2|T| -1.0000000E+10|||grid_xt,grid_yt|
|dynamics|u1000|1000-mb u|m/s|2|T| -1.0000000E+10|||grid_xt,grid_yt|
|dynamics|v1000|1000-mb v|m/s|2|T| -1.0000000E+10|||grid_xt,grid_yt|
|dynamics|t1000|1000-mb temperature|K|2|T| -1.0000000E+10|||grid_xt,grid_yt|
|dynamics|q1000|1000-mb specific humidity|kg/kg|2|T| -1.0000000E+10|||grid_xt,grid_yt|
|dynamics|ql1000|1000-mb cloud water mass mixing ratio|kg/kg|2|T| -1.0000000E+10|||grid_xt,grid_yt|
|dynamics|qi1000|1000-mb cloud ice mass mixing ratio|kg/kg|2|T| -1.0000000E+10|||grid_xt,grid_yt|
|dynamics|qr1000|1000-mb rain mass mixing ratio|kg/kg|2|T| -1.0000000E+10|||grid_xt,grid_yt|
|dynamics|qs1000|1000-mb snow mass mixing ratio|kg/kg|2|T| -1.0000000E+10|||grid_xt,grid_yt|
|dynamics|qg1000|1000-mb graupel mass mixing ratio|kg/kg|2|T| -1.0000000E+10|||grid_xt,grid_yt|
|dynamics|cf1000|1000-mb cloud fraction|1|2|T| -1.0000000E+10|||grid_xt,grid_yt|
|dynamics|omg1000|1000-mb omega|Pa/s|2|T| -1.0000000E+10|||grid_xt,grid_yt|
|dynamics|u_plev|zonal wind|m/sec|3|T| -1.0000000E+10|  -330.0000|   330.0000|grid_xt,grid_yt,plev|
|dynamics|v_plev|meridional wind|m/sec|3|T| -1.0000000E+10|  -330.0000|   330.0000|grid_xt,grid_yt,plev|
|dynamics|t_plev|temperature|K|3|T| -1.0000000E+10|   100.0000|   350.0000|grid_xt,grid_yt,plev|
|dynamics|h_plev|height|m|3|T| -1.0000000E+10|||grid_xt,grid_yt,plev|
|dynamics|q_plev|specific humidity|kg/kg|3|T| -1.0000000E+10|||grid_xt,grid_yt,plev|
|dynamics|ql_plev|cloud water mass mixing ratio|kg/kg|3|T| -1.0000000E+10|||grid_xt,grid_yt,plev|
|dynamics|qi_plev|cloud ice mass mixing ratio|kg/kg|3|T| -1.0000000E+10|||grid_xt,grid_yt,plev|
|dynamics|qr_plev|rain mass mixing ratio|kg/kg|3|T| -1.0000000E+10|||grid_xt,grid_yt,plev|
|dynamics|qs_plev|snow mass mixing ratio|kg/kg|3|T| -1.0000000E+10|||grid_xt,grid_yt,plev|
|dynamics|qg_plev|graupel mass mixing ratio|kg/kg|3|T| -1.0000000E+10|||grid_xt,grid_yt,plev|
|dynamics|cf_plev|cloud fraction|1|3|T| -1.0000000E+10|||grid_xt,grid_yt,plev|
|dynamics|omg_plev|omega|Pa/s|3|T| -1.0000000E+10|||grid_xt,grid_yt,plev|
|dynamics|t_plev_ave|layer-averaged temperature|K|3|T| -1.0000000E+10|||grid_xt,grid_yt,plev_ave|
|dynamics|q_plev_ave|layer-averaged specific humidity|kg/kg|3|T| -1.0000000E+10|||grid_xt,grid_yt,plev_ave|
|dynamics|qv_dt_gfdlmp_plev_ave|layer-averaged water vapor specific humidity tendency from GFDL MP|kg/kg/s|3|T| -1.0000000E+10|||grid_xt,grid_yt,plev_ave|
|dynamics|t_dt_gfdlmp_plev_ave|layer-averaged temperature tendency from GFDL MP|K/s|3|T| -1.0000000E+10|||grid_xt,grid_yt,plev_ave|
|dynamics|qv_dt_phys_plev_ave|layer-averaged water vapor specific humidity tendency from physics|kg/kg/s|3|T| -1.0000000E+10|||grid_xt,grid_yt,plev_ave|
|dynamics|t_dt_phys_plev_ave|layer-averaged temperature tendency from physics|K/s|3|T| -1.0000000E+10|||grid_xt,grid_yt,plev_ave|
|dynamics|tm|mean 300-500 mb temp|K|2|T| -1.0000000E+10|   140.0000|   400.0000|grid_xt,grid_yt|
|dynamics|slp|sea-level pressure|mb|2|T| -1.0000000E+10|   800.0000|   1200.000|grid_xt,grid_yt|
|dynamics|pmask|masking pressure at lowest level|mb|2|T| -1.0000000E+10|||grid_xt,grid_yt|
|dynamics|pmaskv2|masking pressure at lowest level|mb|2|T| -1.0000000E+10|||grid_xt,grid_yt|
|dynamics|cat15|de-pression < 1000|mb|2|T| -1.0000000E+10|||grid_xt,grid_yt|
|dynamics|cat25|de-pression < 980|mb|2|T| -1.0000000E+10|||grid_xt,grid_yt|
|dynamics|cat35|de-pression < 964|mb|2|T| -1.0000000E+10|||grid_xt,grid_yt|
|dynamics|cat45|de-pression < 944|mb|2|T| -1.0000000E+10|||grid_xt,grid_yt|
|dynamics|f15|Cat15 frequency|none|2|T| -1.0000000E+10|||grid_xt,grid_yt|
|dynamics|f25|Cat25 frequency|none|2|T| -1.0000000E+10|||grid_xt,grid_yt|
|dynamics|f35|Cat35 frequency|none|2|T| -1.0000000E+10|||grid_xt,grid_yt|
|dynamics|f45|Cat45 frequency|none|2|T| -1.0000000E+10|||grid_xt,grid_yt|
|dynamics|ucomp|zonal wind|m/sec|3|T| -1.0000000E+10|  -330.0000|   330.0000|grid_xt,grid_yt,pfull|
|dynamics|vcomp|meridional wind|m/sec|3|T| -1.0000000E+10|  -330.0000|   330.0000|grid_xt,grid_yt,pfull|
|dynamics|w|vertical wind|m/sec|3|T| -1.0000000E+10|  -100.0000|   100.0000|grid_xt,grid_yt,pfull|
|dynamics|temp|temperature|K|3|T| -1.0000000E+10|   100.0000|   350.0000|grid_xt,grid_yt,pfull|
|dynamics|theta|potential temperature|K|3|T| -1.0000000E+10|||grid_xt,grid_yt,pfull|
|dynamics|theta_e|equivalent potential temperature|K|3|T| -1.0000000E+10|||grid_xt,grid_yt,pfull|
|dynamics|omega|omega|Pa/s|3|T| -1.0000000E+10|||grid_xt,grid_yt,pfull|
|dynamics|divg|instantaneous divergence|1/s|3|T| -1.0000000E+10|||grid_xt,grid_yt,pfull|
|dynamics|divg_mean|timestep-mean divergence|1/s|3|T| -1.0000000E+10|||grid_xt,grid_yt,pfull|
|dynamics|diss_est|Dissipation estimate|J/kg/s|3|T| -1.0000000E+10| -1.0000000E+07|  1.0000000E+07|grid_xt,grid_yt,pfull|
|dynamics|diss_heat|Dissipative heating rate|K/s|3|T| -1.0000000E+10|||grid_xt,grid_yt,pfull|
|dynamics|hght|height|m|3|T| -1.0000000E+10|||grid_xt,grid_yt,pfull|
|dynamics|rh|Relative Humidity|%|3|T| -1.0000000E+10|||grid_xt,grid_yt,pfull|
|dynamics|delp|pressure thickness (GFS moist-mass)|pa|3|T| -1.0000000E+10|||grid_xt,grid_yt,pfull|
|dynamics|delp_dycore|pressure thickness as seen by the dynamical core|pa|3|T| -1.0000000E+10|||grid_xt,grid_yt,pfull|
|dynamics|delz|height thickness|m|3|T| -1.0000000E+10|||grid_xt,grid_yt,pfull|
|dynamics|pfnh|non-hydrostatic pressure (GFS moist-mass)|pa|3|T| -1.0000000E+10|||grid_xt,grid_yt,pfull|
|dynamics|ppnh|non-hydrostatic pressure perturbation|pa|3|T| -1.0000000E+10|||grid_xt,grid_yt,pfull|
|dynamics|qn|cloud condensate|kg/m/s^2|3|T| -1.0000000E+10|||grid_xt,grid_yt,pfull|
|dynamics|qp|precip condensate|kg/m/s^2|3|T| -1.0000000E+10|||grid_xt,grid_yt,pfull|
|dynamics|qdt|Dqv/Dt: fast moist phys|kg/kg/sec|3|T| -1.0000000E+10|||grid_xt,grid_yt,pfull|
|dynamics|reflectivity|Stoelinga simulated reflectivity|dBz|3|T| -1.0000000E+10|||grid_xt,grid_yt,pfull|
|dynamics|vort|vorticity|1/s|3|T| -1.0000000E+10|||grid_xt,grid_yt,pfull|
|dynamics|pv|potential vorticity|1/s|3|T| -1.0000000E+10|||grid_xt,grid_yt,pfull|
|dynamics|pv350K|350-K potential vorticity; needs x350 scaling|(K m**2) / (kg s)|2|T| -1.0000000E+10|||grid_xt,grid_yt|
|dynamics|pv550K|550-K potential vorticity; needs x550 scaling|(K m**2) / (kg s)|2|T| -1.0000000E+10|||grid_xt,grid_yt|
|dynamics|uq|zonal moisture flux|Kg/Kg*m/sec|3|T| -1.0000000E+10|||grid_xt,grid_yt,pfull|
|dynamics|vq|meridional moisture flux|Kg/Kg*m/sec|3|T| -1.0000000E+10|||grid_xt,grid_yt,pfull|
|dynamics|ut|zonal heat flux|K*m/sec|3|T| -1.0000000E+10|||grid_xt,grid_yt,pfull|
|dynamics|vt|meridional heat flux|K*m/sec|3|T| -1.0000000E+10|||grid_xt,grid_yt,pfull|
|dynamics|uu|zonal flux of zonal wind|(m/sec)^2|3|T| -1.0000000E+10|||grid_xt,grid_yt,pfull|
|dynamics|uv|zonal flux of meridional wind|(m/sec)^2|3|T| -1.0000000E+10|||grid_xt,grid_yt,pfull|
|dynamics|vv|meridional flux of meridional wind|(m/sec)^2|3|T| -1.0000000E+10|||grid_xt,grid_yt,pfull|
|dynamics|uw|vertical zonal momentum flux|N/m**2|3|T| -1.0000000E+10|||grid_xt,grid_yt,pfull|
|dynamics|vw|vertical meridional momentum flux|N/m**|3|T| -1.0000000E+10|||grid_xt,grid_yt,pfull|
|dynamics|wq|vertical moisture flux|Kg/Kg*m/sec|3|T| -1.0000000E+10|||grid_xt,grid_yt,pfull|
|dynamics|wt|vertical heat flux|K*m/sec|3|T| -1.0000000E+10|||grid_xt,grid_yt,pfull|
|dynamics|ww|vertical flux of vertical wind|(m/sec)^2|3|T| -1.0000000E+10|||grid_xt,grid_yt,pfull|
|dynamics|uq_vi|vertical integral of uq|Kg/Kg*m/sec*Pa|2|T| -1.0000000E+10|||grid_xt,grid_yt|
|dynamics|vq_vi|vertical integral of vq|Kg/Kg*m/sec*Pa|2|T| -1.0000000E+10|||grid_xt,grid_yt|
|dynamics|ut_vi|vertical integral of ut|K*m/sec*Pa|2|T| -1.0000000E+10|||grid_xt,grid_yt|
|dynamics|vt_vi|vertical integral of vt|K*m/sec*Pa|2|T| -1.0000000E+10|||grid_xt,grid_yt|
|dynamics|uu_vi|vertical integral of uu|(m/sec)^2*Pa|2|T| -1.0000000E+10|||grid_xt,grid_yt|
|dynamics|uv_vi|vertical integral of uv|(m/sec)^2*Pa|2|T| -1.0000000E+10|||grid_xt,grid_yt|
|dynamics|vv_vi|vertical integral of vv|(m/sec)^2*Pa|2|T| -1.0000000E+10|||grid_xt,grid_yt|
|dynamics|wq_vi|vertical integral of wq|Kg/Kg*m/sec*Pa|2|T| -1.0000000E+10|||grid_xt,grid_yt|
|dynamics|wt_vi|vertical integral of wt|K*m/sec*Pa|2|T| -1.0000000E+10|||grid_xt,grid_yt|
|dynamics|uw_vi|vertical integral of uw|(m/sec)^2*Pa|2|T| -1.0000000E+10|||grid_xt,grid_yt|
|dynamics|vw_vi|vertical integral of vw|(m/sec)^2*Pa|2|T| -1.0000000E+10|||grid_xt,grid_yt|
|dynamics|ww_vi|vertical integral of ww|(m/sec)^2*Pa|2|T| -1.0000000E+10|||grid_xt,grid_yt|
|dynamics|te|Total Energy|J/m/s^2|2|T| -1.0000000E+10|||grid_xt,grid_yt|
|dynamics|ke|Total KE|m^2/s^2|2|T| -1.0000000E+10|||grid_xt,grid_yt|
|dynamics|ws|Terrain W|m/s|2|T| -1.0000000E+10|||grid_xt,grid_yt|
|dynamics|max_reflectivity|Stoelinga simulated maximum (composite) reflectivity|dBz|2|T| -1.0000000E+10|||grid_xt,grid_yt|
|dynamics|base_reflectivity|Stoelinga simulated base (1 km AGL) reflectivity|dBz|2|T| -1.0000000E+10|||grid_xt,grid_yt|
|dynamics|4km_reflectivity|Stoelinga simulated base reflectivity|dBz|2|T| -1.0000000E+10|||grid_xt,grid_yt|
|dynamics|echo_top|Echo top ( <= 18.5 dBz )|m|2|T|  -1000.000|||grid_xt,grid_yt|
|dynamics|m10C_reflectivity|Reflectivity at -10C level|m|2|T| -1.0000000E+10|||grid_xt,grid_yt|
|dynamics|40dBz_height|Height of 40 dBz reflectivity|m|2|T| -1.0000000E+10|||grid_xt,grid_yt|
|dynamics|vorts|surface vorticity|1/s|2|T| -1.0000000E+10|||grid_xt,grid_yt|
|dynamics|us|surface u-wind|m/sec|2|T| -1.0000000E+10|  -200.0000|   200.0000|grid_xt,grid_yt|
|dynamics|vs|surface v-wind|m/sec|2|T| -1.0000000E+10|  -200.0000|   200.0000|grid_xt,grid_yt|
|dynamics|tq|Total water path|kg/m**2|2|T| -1.0000000E+10|||grid_xt,grid_yt|
|dynamics|iw|Ice water path|kg/m**2|2|T| -1.0000000E+10|||grid_xt,grid_yt|
|dynamics|lw|Liquid water path|kg/m**2|2|T| -1.0000000E+10|||grid_xt,grid_yt|
|dynamics|ts|Skin temperature|K|2|T||||grid_xt,grid_yt|
|dynamics|tb|lowest layer temperature|K|2|T||||grid_xt,grid_yt|
|dynamics|ctt|cloud_top temperature|K|2|T|  1.0000000E+10|||grid_xt,grid_yt|
|dynamics|ctp|cloud_top pressure|hPa|2|T|  1.0000000E+10|||grid_xt,grid_yt|
|dynamics|ctz|cloud_top height|hPa|2|T|  -1000.000|||grid_xt,grid_yt|
|dynamics|cape|Convective available potential energy (surface-based)|J/kg|2|T| -1.0000000E+10|||grid_xt,grid_yt|
|dynamics|cin|Convective inhibition (surface-based)|J/kg|2|T| -1.0000000E+10|||grid_xt,grid_yt|
|dynamics|BRN|Bulk Richardson Number|nondim|2|T| -1.0000000E+10|||grid_xt,grid_yt|
|dynamics|shear06|0--6 km shear|m/s|2|T| -1.0000000E+10|||grid_xt,grid_yt|
|dynamics|intqv|Vertically Integrated Water Vapor|kg/m**2|2|T| -1.0000000E+10|||grid_xt,grid_yt|
|dynamics|intql|Vertically Integrated Cloud Water|kg/m**2|2|T| -1.0000000E+10|||grid_xt,grid_yt|
|dynamics|intqi|Vertically Integrated Cloud Ice|kg/m**2|2|T| -1.0000000E+10|||grid_xt,grid_yt|
|dynamics|intqr|Vertically Integrated Rain|kg/m**2|2|T| -1.0000000E+10|||grid_xt,grid_yt|
|dynamics|intqs|Vertically Integrated Snow|kg/m**2|2|T| -1.0000000E+10|||grid_xt,grid_yt|
|dynamics|intqg|Vertically Integrated Graupel|kg/m**2|2|T| -1.0000000E+10|||grid_xt,grid_yt|
|dynamics|acl|Column-averaged Cl mixing ratio|kg/kg|2|T| -1.0000000E+10|||grid_xt,grid_yt|
|dynamics|acl2|Column-averaged Cl2 mixing ratio|kg/kg|2|T| -1.0000000E+10|||grid_xt,grid_yt|
|dynamics|acly|Column-averaged total chlorine mixing ratio|kg/kg|2|T| -1.0000000E+10|||grid_xt,grid_yt|
|dynamics|vort850|850-mb vorticity|1/s|2|T| -1.0000000E+10|||grid_xt,grid_yt|
|dynamics|vort200|200-mb vorticity|1/s|2|T| -1.0000000E+10|||grid_xt,grid_yt|
|dynamics|s200|200-mb wind_speed|m/s|2|T| -1.0000000E+10|||grid_xt,grid_yt|
|dynamics|sl12|12th L wind_speed|m/s|2|T| -1.0000000E+10|||grid_xt,grid_yt|
|dynamics|sl13|13th L wind_speed|m/s|2|T| -1.0000000E+10|||grid_xt,grid_yt|
|dynamics|qn200|200mb condensate|kg/m/s^2|2|T| -1.0000000E+10|||grid_xt,grid_yt|
|dynamics|qn500|500mb condensate|kg/m/s^2|2|T| -1.0000000E+10|||grid_xt,grid_yt|
|dynamics|qn850|850mb condensate|kg/m/s^2|2|T| -1.0000000E+10|||grid_xt,grid_yt|
|dynamics|vort500|500-mb vorticity|1/s|2|T| -1.0000000E+10|||grid_xt,grid_yt|
|dynamics|rain5km|5-km AGL liquid water|kg/kg|2|T| -1.0000000E+10|||grid_xt,grid_yt|
|dynamics|w200|200-mb w-wind|m/s|2|T| -1.0000000E+10|||grid_xt,grid_yt|
|dynamics|w500|500-mb w-wind|m/s|2|T| -1.0000000E+10|||grid_xt,grid_yt|
|dynamics|w700|700-mb w-wind|m/s|2|T| -1.0000000E+10|||grid_xt,grid_yt|
|dynamics|w850|850-mb w-wind|m/s|2|T| -1.0000000E+10|||grid_xt,grid_yt|
|dynamics|w5km|5-km AGL w-wind|m/s|2|T| -1.0000000E+10|||grid_xt,grid_yt|
|dynamics|w2500m|2.5-km AGL w-wind|m/s|2|T| -1.0000000E+10|||grid_xt,grid_yt|
|dynamics|w1km|1-km AGL w-wind|m/s|2|T| -1.0000000E+10|||grid_xt,grid_yt|
|dynamics|wmaxup|column-maximum updraft (below 100 mb)|m/s|2|T| -1.0000000E+10|||grid_xt,grid_yt|
|dynamics|wmaxdn|column-maximum downdraft (below 100 mb)|m/s|2|T| -1.0000000E+10|||grid_xt,grid_yt|
|dynamics|x850|850-mb vertical comp. of helicity|m/s**2|2|T| -1.0000000E+10|||grid_xt,grid_yt|
|dynamics|srh01|0-1 km Storm Relative Helicity|m/s**2|2|T| -1.0000000E+10|||grid_xt,grid_yt|
|dynamics|srh03|0-3 km Storm Relative Helicity|m/s**2|2|T| -1.0000000E+10|||grid_xt,grid_yt|
|dynamics|ustm|u Component of Storm Motion|m/s|2|T| -1.0000000E+10|||grid_xt,grid_yt|
|dynamics|vstm|v Component of Storm Motion|m/s|2|T| -1.0000000E+10|||grid_xt,grid_yt|
|dynamics|srh25|2-5 km Storm Relative Helicity|m/s**2|2|T| -1.0000000E+10|||grid_xt,grid_yt|
|dynamics|uh03|0-3 km Updraft Helicity|m/s**2|2|T| -1.0000000E+10|||grid_xt,grid_yt|
|dynamics|uh25|2-5 km Updraft Helicity|m/s**2|2|T| -1.0000000E+10|||grid_xt,grid_yt|
|dynamics|w100m|100-m AGL w-wind|m/s|2|T| -1.0000000E+10|||grid_xt,grid_yt|
|dynamics|u100m|100-m AGL u-wind|m/s|2|T| -1.0000000E+10|||grid_xt,grid_yt|
|dynamics|v100m|100-m AGL v-wind|m/s|2|T| -1.0000000E+10|||grid_xt,grid_yt|
|dynamics|rh10|10-mb relative humidity|%|2|T| -1.0000000E+10|||grid_xt,grid_yt|
|dynamics|rh50|50-mb relative humidity|%|2|T| -1.0000000E+10|||grid_xt,grid_yt|
|dynamics|rh100|100-mb relative humidity|%|2|T| -1.0000000E+10|||grid_xt,grid_yt|
|dynamics|rh200|200-mb relative humidity|%|2|T| -1.0000000E+10|||grid_xt,grid_yt|
|dynamics|rh250|250-mb relative humidity|%|2|T| -1.0000000E+10|||grid_xt,grid_yt|
|dynamics|rh300|300-mb relative humidity|%|2|T| -1.0000000E+10|||grid_xt,grid_yt|
|dynamics|rh500|500-mb relative humidity|%|2|T| -1.0000000E+10|||grid_xt,grid_yt|
|dynamics|rh700|700-mb relative humidity|%|2|T| -1.0000000E+10|||grid_xt,grid_yt|
|dynamics|rh850|850-mb relative humidity|%|2|T| -1.0000000E+10|||grid_xt,grid_yt|
|dynamics|rh925|925-mb relative humidity|%|2|T| -1.0000000E+10|||grid_xt,grid_yt|
|dynamics|rh1000|1000-mb relative humidity|%|2|T| -1.0000000E+10|||grid_xt,grid_yt|
|dynamics|dp10|10-mb dew point|K|2|T| -1.0000000E+10|||grid_xt,grid_yt|
|dynamics|dp50|50-mb dew point|K|2|T| -1.0000000E+10|||grid_xt,grid_yt|
|dynamics|dp100|100-mb dew point|K|2|T| -1.0000000E+10|||grid_xt,grid_yt|
|dynamics|dp200|200-mb dew point|K|2|T| -1.0000000E+10|||grid_xt,grid_yt|
|dynamics|dp250|250-mb dew point|K|2|T| -1.0000000E+10|||grid_xt,grid_yt|
|dynamics|dp300|300-mb dew point|K|2|T| -1.0000000E+10|||grid_xt,grid_yt|
|dynamics|dp500|500-mb dew point|K|2|T| -1.0000000E+10|||grid_xt,grid_yt|
|dynamics|dp700|700-mb dew point|K|2|T| -1.0000000E+10|||grid_xt,grid_yt|
|dynamics|dp850|850-mb dew point|K|2|T| -1.0000000E+10|||grid_xt,grid_yt|
|dynamics|dp925|925-mb dew point|K|2|T| -1.0000000E+10|||grid_xt,grid_yt|
|dynamics|dp1000|1000-mb dew point|K|2|T| -1.0000000E+10|||grid_xt,grid_yt|
|dynamics|theta_e100|100-mb equivalent potential temperature|K|2|T| -1.0000000E+10|||grid_xt,grid_yt|
|dynamics|theta_e200|200-mb equivalent potential temperature|K|2|T| -1.0000000E+10|||grid_xt,grid_yt|
|dynamics|theta_e250|250-mb equivalent potential temperature|K|2|T| -1.0000000E+10|||grid_xt,grid_yt|
|dynamics|theta_e300|300-mb equivalent potential temperature|K|2|T| -1.0000000E+10|||grid_xt,grid_yt|
|dynamics|theta_e500|500-mb equivalent potential temperature|K|2|T| -1.0000000E+10|||grid_xt,grid_yt|
|dynamics|theta_e700|700-mb equivalent potential temperature|K|2|T| -1.0000000E+10|||grid_xt,grid_yt|
|dynamics|theta_e850|850-mb equivalent potential temperature|K|2|T| -1.0000000E+10|||grid_xt,grid_yt|
|dynamics|theta_e925|925-mb equivalent potential temperature|K|2|T| -1.0000000E+10|||grid_xt,grid_yt|
|dynamics|theta_e1000|1000-mb equivalent potential temperature|K|2|T| -1.0000000E+10|||grid_xt,grid_yt|
|dynamics|rh10_cmip|10-mb relative humidity (CMIP)|%|2|T| -1.0000000E+10|||grid_xt,grid_yt|
|dynamics|rh50_cmip|50-mb relative humidity (CMIP)|%|2|T| -1.0000000E+10|||grid_xt,grid_yt|
|dynamics|rh100_cmip|100-mb relative humidity (CMIP)|%|2|T| -1.0000000E+10|||grid_xt,grid_yt|
|dynamics|rh250_cmip|250-mb relative humidity (CMIP)|%|2|T| -1.0000000E+10|||grid_xt,grid_yt|
|dynamics|rh300_cmip|300-mb relative humidity (CMIP)|%|2|T| -1.0000000E+10|||grid_xt,grid_yt|
|dynamics|rh500_cmip|500-mb relative humidity (CMIP)|%|2|T| -1.0000000E+10|||grid_xt,grid_yt|
|dynamics|rh700_cmip|700-mb relative humidity (CMIP)|%|2|T| -1.0000000E+10|||grid_xt,grid_yt|
|dynamics|rh850_cmip|850-mb relative humidity (CMIP)|%|2|T| -1.0000000E+10|||grid_xt,grid_yt|
|dynamics|rh925_cmip|925-mb relative humidity (CMIP)|%|2|T| -1.0000000E+10|||grid_xt,grid_yt|
|dynamics|rh1000_cmip|1000-mb relative humidity (CMIP)|%|2|T| -1.0000000E+10|||grid_xt,grid_yt|
|dynamics|sphum|specific humidity|kg/kg|3|T| -1.0000000E+10|||grid_xt,grid_yt,pfull|
|dynamics|liq_wat|cloud water mixing ratio|kg/kg|3|T| -1.0000000E+10|||grid_xt,grid_yt,pfull|
|dynamics|rainwat|rain mixing ratio|kg/kg|3|T| -1.0000000E+10|||grid_xt,grid_yt,pfull|
|dynamics|ice_wat|cloud ice mixing ratio|kg/kg|3|T| -1.0000000E+10|||grid_xt,grid_yt,pfull|
|dynamics|snowwat|snow mixing ratio|kg/kg|3|T| -1.0000000E+10|||grid_xt,grid_yt,pfull|
|dynamics|graupel|graupel mixing ratio|kg/kg|3|T| -1.0000000E+10|||grid_xt,grid_yt,pfull|
|dynamics|sgs_tke|tke|m**2/s**2|3|T| -1.0000000E+10|||grid_xt,grid_yt,pfull|
|dynamics|o3mr|ozone mixing ratio|kg/kg|3|T| -1.0000000E+10|||grid_xt,grid_yt,pfull|
|dynamics|cld_amt|cloud amount|1|3|T| -1.0000000E+10|||grid_xt,grid_yt,pfull|
|gfs_dyn|ucomp|zonal wind|m/sec|3|T| -1.0000000E+10|  -330.0000|   330.0000|grid_xt,grid_yt,pfull|
|gfs_dyn|vcomp|meridional wind|m/sec|3|T| -1.0000000E+10|  -330.0000|   330.0000|grid_xt,grid_yt,pfull|
|gfs_dyn|pfnh|non-hydrostatic pressure|pa|3|T| -1.0000000E+10|||grid_xt,grid_yt,pfull|
|gfs_dyn|w|vertical wind|m/sec|3|T| -1.0000000E+10|  -100.0000|   100.0000|grid_xt,grid_yt,pfull|
|gfs_dyn|delz|height thickness|m|3|T| -1.0000000E+10|||grid_xt,grid_yt,pfull|
|gfs_dyn|omga|Vertical pressure velocity|pa/sec|3|T| -1.0000000E+10|||grid_xt,grid_yt,pfull|
|gfs_dyn|temp|temperature|K|3|T| -1.0000000E+10|   100.0000|   350.0000|grid_xt,grid_yt,pfull|
|gfs_dyn|delp|pressure thickness|pa|3|T| -1.0000000E+10|||grid_xt,grid_yt,pfull|
|gfs_dyn|diss_est|dissipation estimate|none|3|T| -1.0000000E+10| -1.0000000E+07|  1.0000000E+07|grid_xt,grid_yt,pfull|
|gfs_dyn|sphum|specific humidity|kg/kg|3|T| -1.0000000E+10|||grid_xt,grid_yt,pfull|
|gfs_dyn|liq_wat|cloud water mixing ratio|kg/kg|3|T| -1.0000000E+10|||grid_xt,grid_yt,pfull|
|gfs_dyn|rainwat|rain mixing ratio|kg/kg|3|T| -1.0000000E+10|||grid_xt,grid_yt,pfull|
|gfs_dyn|ice_wat|cloud ice mixing ratio|kg/kg|3|T| -1.0000000E+10|||grid_xt,grid_yt,pfull|
|gfs_dyn|snowwat|snow mixing ratio|kg/kg|3|T| -1.0000000E+10|||grid_xt,grid_yt,pfull|
|gfs_dyn|graupel|graupel mixing ratio|kg/kg|3|T| -1.0000000E+10|||grid_xt,grid_yt,pfull|
|gfs_dyn|sgs_tke|tke|m**2/s**2|3|T| -1.0000000E+10|||grid_xt,grid_yt,pfull|
|gfs_dyn|o3mr|ozone mixing ratio|kg/kg|3|T| -1.0000000E+10|||grid_xt,grid_yt,pfull|
|gfs_dyn|cld_amt|cloud amount|1|3|T| -1.0000000E+10|||grid_xt,grid_yt,pfull|
|gfs_dyn|ps|surface pressure|pa|2|T| -1.0000000E+10|||grid_xt,grid_yt|
|gfs_dyn|hs|surface geopotential height|gpm|2|T| -1.0000000E+10|||grid_xt,grid_yt|
|gfs_dyn|reflectivity|Stoelinga simulated reflectivity|dBz|3|T| -1.0000000E+10|||grid_xt,grid_yt,pfull|
|gfs_dyn|ustm|u comp of storm motion|m/s|2|T| -1.0000000E+10|||grid_xt,grid_yt|
|gfs_dyn|vstm|v comp of storm motion|m/s|2|T| -1.0000000E+10|||grid_xt,grid_yt|
|gfs_dyn|srh01|0-1km storm rel. helicity|m/s**2|2|T| -1.0000000E+10|||grid_xt,grid_yt|
|gfs_dyn|srh03|0-3km storm rel. helicity|m/s**2|2|T| -1.0000000E+10|||grid_xt,grid_yt|
|gfs_dyn|maxvort01|Max hourly 0-1km vert vorticity|1/s|2|T| -1.0000000E+10|||grid_xt,grid_yt|
|gfs_dyn|maxvort02|Max hourly 0-2km vert vorticity|1/s|2|T| -1.0000000E+10|||grid_xt,grid_yt|
|gfs_dyn|maxvorthy1|Max hourly hybrid lev1 vert. vorticity|1/s|2|T| -1.0000000E+10|||grid_xt,grid_yt|
|gfs_dyn|wmaxup|Max hourly updraft velocity|m/s|2|T| -1.0000000E+10|||grid_xt,grid_yt|
|gfs_dyn|wmaxdn|Max hourly downdraft velocity|m/s|2|T| -1.0000000E+10|||grid_xt,grid_yt|
|gfs_dyn|uhmax03|Max hourly max 0-3km updraft helicity|m/s**2|2|T| -1.0000000E+10|||grid_xt,grid_yt|
|gfs_dyn|uhmin03|Max hourly min 0-3km updraft helicity|m/s**2|2|T| -1.0000000E+10|||grid_xt,grid_yt|
|gfs_dyn|uhmax25|Max hourly max 2-5km updraft helicity|m/s**2|2|T| -1.0000000E+10|||grid_xt,grid_yt|
|gfs_dyn|uhmin25|Max hourly min 2-5km updraft helicity|m/s**2|2|T| -1.0000000E+10|||grid_xt,grid_yt|
|gfs_phys|ALBDOsfc|surface albedo|%|2|T|  9.9900001E+20|||grid_xt,grid_yt|
|gfs_phys|USWRFsfc|Interval-averaged unadjusted upward shortwave flux at the surface|W/m**2|2|T|  9.9900001E+20|||grid_xt,grid_yt|
|gfs_phys|DSWRFsfc|Interval-averaged unadjusted downward shortwave flux at the surface|W/m**2|2|T|  9.9900001E+20|||grid_xt,grid_yt|
|gfs_phys|DLWRFsfc|Interval-averaged unadjusted downward longwave flux at the surface|W/m**2|2|T|  9.9900001E+20|||grid_xt,grid_yt|
|gfs_phys|ULWRFsfc|Interval-averaged unadjusted upward longwave flux at the surface|W/m**2|2|T|  9.9900001E+20|||grid_xt,grid_yt|
|gfs_phys|duvb_ave|UV-B Downward Solar Flux|W/m**2|2|T|  9.9900001E+20|||grid_xt,grid_yt|
|gfs_phys|cduvb_ave|Clear sky UV-B Downward Solar Flux|W/m**2|2|T|  9.9900001E+20|||grid_xt,grid_yt|
|gfs_phys|vbdsf_ave|Visible Beam Downward Solar Flux|W/m**2|2|T|  9.9900001E+20|||grid_xt,grid_yt|
|gfs_phys|vddsf_ave|Visible Diffuse Downward Solar Flux|W/m**2|2|T|  9.9900001E+20|||grid_xt,grid_yt|
|gfs_phys|nbdsf_ave|Near IR Beam Downward Solar Flux|W/m**2|2|T|  9.9900001E+20|||grid_xt,grid_yt|
|gfs_phys|nddsf_ave|Near IR Diffuse Downward Solar Flux|W/m**2|2|T|  9.9900001E+20|||grid_xt,grid_yt|
|gfs_phys|csulf_avetoa|Clear Sky Upward Long Wave Flux at toa|W/m**2|2|T|  9.9900001E+20|||grid_xt,grid_yt|
|gfs_phys|csusf_avetoa|Clear Sky Upward Short Wave Flux at toa|W/m**2|2|T|  9.9900001E+20|||grid_xt,grid_yt|
|gfs_phys|csdlf_ave|Clear Sky Downward Long Wave Flux|W/m**2|2|T|  9.9900001E+20|||grid_xt,grid_yt|
|gfs_phys|csusf_ave|Clear Sky Upward Short Wave Flux|W/m**2|2|T|  9.9900001E+20|||grid_xt,grid_yt|
|gfs_phys|csdsf_ave|Clear Sky Downward Short Wave Flux|W/m**2|2|T|  9.9900001E+20|||grid_xt,grid_yt|
|gfs_phys|csulf_ave|Clear Sky Upward Long Wave Flux|W/m**2|2|T|  9.9900001E+20|||grid_xt,grid_yt|
|gfs_phys|DSWRFtoa|top of atmos downward shortwave flux [W/m**2]|W/m**2|2|T|  9.9900001E+20|||grid_xt,grid_yt|
|gfs_phys|USWRFtoa|top of atmos upward shortwave flux [W/m**2]|W/m**2|2|T|  9.9900001E+20|||grid_xt,grid_yt|
|gfs_phys|ULWRFtoa|top of atmos upward longwave flux [W/m**2]|W/m**2|2|T|  9.9900001E+20|||grid_xt,grid_yt|
|gfs_phys|TCDCclm|atmos column total cloud cover [%]|%|2|T|  9.9900001E+20|||grid_xt,grid_yt|
|gfs_phys|TCDCbndcl|boundary layer cloud layer total cloud cover|%|2|T|  9.9900001E+20|||grid_xt,grid_yt|
|gfs_phys|TCDCcnvcl|convective cloud layer total cloud cover|%|2|T|  9.9900001E+20|||grid_xt,grid_yt|
|gfs_phys|PREScnvclt|pressure at convective cloud top level|pa|2|T|  9.9900001E+20|||grid_xt,grid_yt|
|gfs_phys|PREScnvclb|pressure at convective cloud bottom level|pa|2|T|  9.9900001E+20|||grid_xt,grid_yt|
|gfs_phys|TCDChcl|high cloud level total cloud cover [%]|%|2|T|  9.9900001E+20|||grid_xt,grid_yt|
|gfs_phys|PRES_avehct|pressure high cloud top level|%|2|T|  9.9900001E+20|||grid_xt,grid_yt|
|gfs_phys|PRES_avehcb|pressure high cloud bottom level|%|2|T|  9.9900001E+20|||grid_xt,grid_yt|
|gfs_phys|TEMP_avehct|temperature high cloud top level|K|2|T|  9.9900001E+20|||grid_xt,grid_yt|
|gfs_phys|TCDCmcl|mid cloud level total cloud cover|%|2|T|  9.9900001E+20|||grid_xt,grid_yt|
|gfs_phys|PRES_avemct|pressure middle cloud top level|%|2|T|  9.9900001E+20|||grid_xt,grid_yt|
|gfs_phys|PRES_avemcb|pressure middle cloud bottom level|%|2|T|  9.9900001E+20|||grid_xt,grid_yt|
|gfs_phys|TEMP_avemct|temperature middle cloud top level|K|2|T|  9.9900001E+20|||grid_xt,grid_yt|
|gfs_phys|TCDClcl|low cloud level total cloud cover|%|2|T|  9.9900001E+20|||grid_xt,grid_yt|
|gfs_phys|PRES_avelct|pressure low cloud top level|%|2|T|  9.9900001E+20|||grid_xt,grid_yt|
|gfs_phys|PRES_avelcb|pressure low cloud bottom level|%|2|T|  9.9900001E+20|||grid_xt,grid_yt|
|gfs_phys|TEMP_avelct|temperature low cloud top level|K|2|T|  9.9900001E+20|||grid_xt,grid_yt|
|gfs_phys|fluxr_01|fluxr diagnostic 01 - GFS radiation|XXX|2|T|  9.9900001E+20|||grid_xt,grid_yt|
|gfs_phys|fluxr_02|fluxr diagnostic 02 - GFS radiation|XXX|2|T|  9.9900001E+20|||grid_xt,grid_yt|
|gfs_phys|fluxr_03|fluxr diagnostic 03 - GFS radiation|XXX|2|T|  9.9900001E+20|||grid_xt,grid_yt|
|gfs_phys|fluxr_04|fluxr diagnostic 04 - GFS radiation|XXX|2|T|  9.9900001E+20|||grid_xt,grid_yt|
|gfs_phys|fluxr_05|fluxr diagnostic 05 - GFS radiation|XXX|2|T|  9.9900001E+20|||grid_xt,grid_yt|
|gfs_phys|fluxr_06|fluxr diagnostic 06 - GFS radiation|XXX|2|T|  9.9900001E+20|||grid_xt,grid_yt|
|gfs_phys|fluxr_07|fluxr diagnostic 07 - GFS radiation|XXX|2|T|  9.9900001E+20|||grid_xt,grid_yt|
|gfs_phys|fluxr_08|fluxr diagnostic 08 - GFS radiation|XXX|2|T|  9.9900001E+20|||grid_xt,grid_yt|
|gfs_phys|fluxr_09|fluxr diagnostic 09 - GFS radiation|XXX|2|T|  9.9900001E+20|||grid_xt,grid_yt|
|gfs_phys|fluxr_10|fluxr diagnostic 10 - GFS radiation|XXX|2|T|  9.9900001E+20|||grid_xt,grid_yt|
|gfs_phys|fluxr_11|fluxr diagnostic 11 - GFS radiation|XXX|2|T|  9.9900001E+20|||grid_xt,grid_yt|
|gfs_phys|fluxr_12|fluxr diagnostic 12 - GFS radiation|XXX|2|T|  9.9900001E+20|||grid_xt,grid_yt|
|gfs_phys|fluxr_13|fluxr diagnostic 13 - GFS radiation|XXX|2|T|  9.9900001E+20|||grid_xt,grid_yt|
|gfs_phys|fluxr_14|fluxr diagnostic 14 - GFS radiation|XXX|2|T|  9.9900001E+20|||grid_xt,grid_yt|
|gfs_phys|fluxr_15|fluxr diagnostic 15 - GFS radiation|XXX|2|T|  9.9900001E+20|||grid_xt,grid_yt|
|gfs_phys|fluxr_16|fluxr diagnostic 16 - GFS radiation|XXX|2|T|  9.9900001E+20|||grid_xt,grid_yt|
|gfs_phys|fluxr_17|fluxr diagnostic 17 - GFS radiation|XXX|2|T|  9.9900001E+20|||grid_xt,grid_yt|
|gfs_phys|fluxr_18|fluxr diagnostic 18 - GFS radiation|XXX|2|T|  9.9900001E+20|||grid_xt,grid_yt|
|gfs_phys|fluxr_19|fluxr diagnostic 19 - GFS radiation|XXX|2|T|  9.9900001E+20|||grid_xt,grid_yt|
|gfs_phys|fluxr_20|fluxr diagnostic 20 - GFS radiation|XXX|2|T|  9.9900001E+20|||grid_xt,grid_yt|
|gfs_phys|fluxr_21|fluxr diagnostic 21 - GFS radiation|XXX|2|T|  9.9900001E+20|||grid_xt,grid_yt|
|gfs_phys|fluxr_22|fluxr diagnostic 22 - GFS radiation|XXX|2|T|  9.9900001E+20|||grid_xt,grid_yt|
|gfs_phys|fluxr_23|fluxr diagnostic 23 - GFS radiation|XXX|2|T|  9.9900001E+20|||grid_xt,grid_yt|
|gfs_phys|fluxr_24|fluxr diagnostic 24 - GFS radiation|XXX|2|T|  9.9900001E+20|||grid_xt,grid_yt|
|gfs_phys|fluxr_25|fluxr diagnostic 25 - GFS radiation|XXX|2|T|  9.9900001E+20|||grid_xt,grid_yt|
|gfs_phys|fluxr_26|fluxr diagnostic 26 - GFS radiation|XXX|2|T|  9.9900001E+20|||grid_xt,grid_yt|
|gfs_phys|fluxr_27|fluxr diagnostic 27 - GFS radiation|XXX|2|T|  9.9900001E+20|||grid_xt,grid_yt|
|gfs_phys|fluxr_28|fluxr diagnostic 28 - GFS radiation|XXX|2|T|  9.9900001E+20|||grid_xt,grid_yt|
|gfs_phys|fluxr_29|fluxr diagnostic 29 - GFS radiation|XXX|2|T|  9.9900001E+20|||grid_xt,grid_yt|
|gfs_phys|fluxr_30|fluxr diagnostic 30 - GFS radiation|XXX|2|T|  9.9900001E+20|||grid_xt,grid_yt|
|gfs_phys|fluxr_31|fluxr diagnostic 31 - GFS radiation|XXX|2|T|  9.9900001E+20|||grid_xt,grid_yt|
|gfs_phys|fluxr_32|fluxr diagnostic 32 - GFS radiation|XXX|2|T|  9.9900001E+20|||grid_xt,grid_yt|
|gfs_phys|fluxr_33|fluxr diagnostic 33 - GFS radiation|XXX|2|T|  9.9900001E+20|||grid_xt,grid_yt|
|gfs_phys|fluxr_34|fluxr diagnostic 34 - GFS radiation|XXX|2|T|  9.9900001E+20|||grid_xt,grid_yt|
|gfs_phys|fluxr_35|fluxr diagnostic 35 - GFS radiation|XXX|2|T|  9.9900001E+20|||grid_xt,grid_yt|
|gfs_phys|fluxr_36|fluxr diagnostic 36 - GFS radiation|XXX|2|T|  9.9900001E+20|||grid_xt,grid_yt|
|gfs_phys|fluxr_37|fluxr diagnostic 37 - GFS radiation|XXX|2|T|  9.9900001E+20|||grid_xt,grid_yt|
|gfs_phys|fluxr_38|fluxr diagnostic 38 - GFS radiation|XXX|2|T|  9.9900001E+20|||grid_xt,grid_yt|
|gfs_phys|fluxr_39|fluxr diagnostic 39 - GFS radiation|XXX|2|T|  9.9900001E+20|||grid_xt,grid_yt|
|gfs_phys|cloud_01|cloud diagnostic 01 - GFS radiation|XXX|3|T|  9.9900001E+20|||grid_xt,grid_yt,pfull|
|gfs_phys|cloud_02|cloud diagnostic 02 - GFS radiation|XXX|3|T|  9.9900001E+20|||grid_xt,grid_yt,pfull|
|gfs_phys|cloud_03|cloud diagnostic 03 - GFS radiation|XXX|3|T|  9.9900001E+20|||grid_xt,grid_yt,pfull|
|gfs_phys|cloud_04|cloud diagnostic 04 - GFS radiation|XXX|3|T|  9.9900001E+20|||grid_xt,grid_yt,pfull|
|gfs_phys|cloud_05|cloud diagnostic 05 - GFS radiation|XXX|3|T|  9.9900001E+20|||grid_xt,grid_yt,pfull|
|gfs_phys|cloud_06|cloud diagnostic 06 - GFS radiation|XXX|3|T|  9.9900001E+20|||grid_xt,grid_yt,pfull|
|gfs_phys|cloud_07|cloud diagnostic 07 - GFS radiation|XXX|3|T|  9.9900001E+20|||grid_xt,grid_yt,pfull|
|gfs_phys|cloud_08|cloud diagnostic 08 - GFS radiation|XXX|3|T|  9.9900001E+20|||grid_xt,grid_yt,pfull|
|gfs_phys|sw_upfxc|total sky upward sw flux at toa - GFS radiation|w/m**2|2|T|  9.9900001E+20|||grid_xt,grid_yt|
|gfs_phys|sw_dnfxc|total sky downward sw flux at toa - GFS radiation|w/m**2|2|T|  9.9900001E+20|||grid_xt,grid_yt|
|gfs_phys|sw_upfx0|clear sky upward sw flux at toa - GFS radiation|w/m**2|2|T|  9.9900001E+20|||grid_xt,grid_yt|
|gfs_phys|lw_upfxc|total sky upward lw flux at toa - GFS radiation|w/m**2|2|T|  9.9900001E+20|||grid_xt,grid_yt|
|gfs_phys|lw_upfx0|clear sky upward lw flux at toa - GFS radiation|w/m**2|2|T|  9.9900001E+20|||grid_xt,grid_yt|
|gfs_phys|ssrun_acc|surface storm water runoff - GFS lsm|kg/m**2|2|T|  9.9900001E+20|||grid_xt,grid_yt|
|gfs_phys|evbs_ave|Direct Evaporation from Bare Soil - GFS lsm|W/m**2|2|T|  9.9900001E+20|||grid_xt,grid_yt|
|gfs_phys|evcw_ave|Canopy water evaporation - GFS lsm|W/m**2|2|T|  9.9900001E+20|||grid_xt,grid_yt|
|gfs_phys|snohf_ave|Snow Phase Change Heat Flux - GFS lsm|W/m**2|2|T|  9.9900001E+20|||grid_xt,grid_yt|
|gfs_phys|trans_ave|transpiration - GFS lsm|W/m**2|2|T|  9.9900001E+20|||grid_xt,grid_yt|
|gfs_phys|sbsno_ave|Sublimation (evaporation from snow) - GFS lsm|W/m**2|2|T|  9.9900001E+20|||grid_xt,grid_yt|
|gfs_phys|snowc_ave|snow cover - GFS lsm|%|2|T|  9.9900001E+20|||grid_xt,grid_yt|
|gfs_phys|soilm|total column soil moisture content [kg/m**2]|kg/m**2|2|T|  9.9900001E+20|||grid_xt,grid_yt|
|gfs_phys|tmpmin2m|min temperature at 2m height|K|2|T|  9.9900001E+20|||grid_xt,grid_yt|
|gfs_phys|tmpmax2m|max temperature at 2m height|K|2|T|  9.9900001E+20|||grid_xt,grid_yt|
|gfs_phys|dusfc|surface zonal momentum flux [N/m**2]|N/m**2|2|T|  9.9900001E+20|||grid_xt,grid_yt|
|gfs_phys|dvsfc|surface meridional momentum flux|N/m**2|2|T|  9.9900001E+20|||grid_xt,grid_yt|
|gfs_phys|shtfl_ave|surface sensible heat flux|w/m**2|2|T|  9.9900001E+20|||grid_xt,grid_yt|
|gfs_phys|lhtfl_ave|surface latent heat flux [W/m**2]|w/m**2|2|T|  9.9900001E+20|||grid_xt,grid_yt|
|gfs_phys|totprcp_ave|surface precipitation rate [kg/m**2/s]|kg/m**2/s|2|T|  9.9900001E+20|||grid_xt,grid_yt|
|gfs_phys|totprcpb_ave|bucket surface precipitation rate|kg/m**2/s|2|T|  9.9900001E+20|||grid_xt,grid_yt|
|gfs_phys|gflux_ave|surface ground heat flux [W/m**2]|W/m**2|2|T|  9.9900001E+20|||grid_xt,grid_yt|
|gfs_phys|DSWRF|Interval-averaged zenith-angle-adjusted downward shortwave flux at the surface|w/m**2|2|T|  9.9900001E+20|||grid_xt,grid_yt|
|gfs_phys|USWRF|Interval-averaged zenith-angle-adjusted upward shortwave flux at the surface|w/m**2|2|T|  9.9900001E+20|||grid_xt,grid_yt|
|gfs_phys|DLWRF|Interval-averaged surface-temperature-adjusted downward longwave flux at the surface|w/m**2|2|T|  9.9900001E+20|||grid_xt,grid_yt|
|gfs_phys|ULWRF|Interval-averaged surface-temperature-adjusted upward longwave flux at the surface|w/m**2|2|T|  9.9900001E+20|||grid_xt,grid_yt|
|gfs_phys|sunsd_acc|sunshine duration time|s|2|T|  9.9900001E+20|||grid_xt,grid_yt|
|gfs_phys|watr_acc|total water runoff|kg/m**2|2|T|  9.9900001E+20|||grid_xt,grid_yt|
|gfs_phys|pevpr_ave|averaged potential evaporation rate|W/M**2|2|T|  9.9900001E+20|||grid_xt,grid_yt|
|gfs_phys|cwork_ave|cloud work function (valid only with sas)|J/kg|2|T|  9.9900001E+20|||grid_xt,grid_yt|
|gfs_phys|u-gwd_ave|surface zonal gravity wave stress|N/m**2|2|T|  9.9900001E+20|||grid_xt,grid_yt|
|gfs_phys|v-gwd_ave|surface meridional gravity wave stress|N/m**2|2|T|  9.9900001E+20|||grid_xt,grid_yt|
|gfs_phys|psmean|surface pressure|kPa|2|T|  9.9900001E+20|||grid_xt,grid_yt|
|gfs_phys|cnvprcp_ave|averaged surface convective precipitation rate|kg/m**2/s|2|T|  9.9900001E+20|||grid_xt,grid_yt|
|gfs_phys|cnvprcpb_ave|averaged bucket surface convective precipitation rate|kg/m**2/s|2|T|  9.9900001E+20|||grid_xt,grid_yt|
|gfs_phys|cnvprcp|surface convective precipitation rate|kg/m**2/s|2|T|  9.9900001E+20|||grid_xt,grid_yt|
|gfs_phys|spfhmin2m|minimum specific humidity at 2m height|kg/kg|2|T|  9.9900001E+20|||grid_xt,grid_yt|
|gfs_phys|spfhmax2m|maximum specific humidity at 2m height|kg/kg|2|T|  9.9900001E+20|||grid_xt,grid_yt|
|gfs_phys|u10mmax|maximum (magnitude) u-wind|m/s|2|T|  9.9900001E+20|||grid_xt,grid_yt|
|gfs_phys|v10mmax|maximum (magnitude) v-wind|m/s|2|T|  9.9900001E+20|||grid_xt,grid_yt|
|gfs_phys|wind10mmax|maximum wind speed|m/s|2|T|  9.9900001E+20|||grid_xt,grid_yt|
|gfs_phys|rain|total rain at this time step|XXX|2|T|  9.9900001E+20|||grid_xt,grid_yt|
|gfs_phys|rainc|convective rain at this time step|XXX|2|T|  9.9900001E+20|||grid_xt,grid_yt|
|gfs_phys|ice|ice fall at this time step|XXX|2|T|  9.9900001E+20|||grid_xt,grid_yt|
|gfs_phys|snow|snow fall at this time step|XXX|2|T|  9.9900001E+20|||grid_xt,grid_yt|
|gfs_phys|graupel|graupel fall at this time step|XXX|2|T|  9.9900001E+20|||grid_xt,grid_yt|
|gfs_phys|totice_ave|surface ice precipitation rate|kg/m**2/s|2|T|  9.9900001E+20|||grid_xt,grid_yt|
|gfs_phys|toticeb_ave|bucket surface ice precipitation rate|kg/m**2/s|2|T|  9.9900001E+20|||grid_xt,grid_yt|
|gfs_phys|totsnw_ave|surface snow precipitation rate|kg/m**2/s|2|T|  9.9900001E+20|||grid_xt,grid_yt|
|gfs_phys|totsnwb_ave|bucket surface snow precipitation rate|kg/m**2/s|2|T|  9.9900001E+20|||grid_xt,grid_yt|
|gfs_phys|totgrp_ave|surface graupel precipitation rate|kg/m**2/s|2|T|  9.9900001E+20|||grid_xt,grid_yt|
|gfs_phys|totgrpb_ave|bucket surface graupel precipitation rate|kg/m**2/s|2|T|  9.9900001E+20|||grid_xt,grid_yt|
|gfs_phys|u10m|10 meter u wind [m/s]|m/s|2|T|  9.9900001E+20|||grid_xt,grid_yt|
|gfs_phys|v10m|10 meter v wind [m/s]|m/s|2|T|  9.9900001E+20|||grid_xt,grid_yt|
|gfs_phys|dpt2m|2 meter dew point temperature [K]|K|2|T|  9.9900001E+20|||grid_xt,grid_yt|
|gfs_phys|hgt_hyblev1|layer 1 height|m|2|T|  9.9900001E+20|||grid_xt,grid_yt|
|gfs_phys|psurf|surface pressure [Pa]|Pa|2|T|  9.9900001E+20|||grid_xt,grid_yt|
|gfs_phys|hpbl|surface planetary boundary layer height [m]|m|2|T|  9.9900001E+20|||grid_xt,grid_yt|
|gfs_phys|hgamt|ysu counter-gradient heat flux factor|XXX|2|T|  9.9900001E+20|||grid_xt,grid_yt|
|gfs_phys|hfxpbl|ysu entrainment heat flux factor|XXX|2|T|  9.9900001E+20|||grid_xt,grid_yt|
|gfs_phys|xmb_shal|cloud base mass flux from mass-flux shal cnv|XXX|2|T|  9.9900001E+20|||grid_xt,grid_yt|
|gfs_phys|tfac_shal|Tadv/Tcnv factor from  mass-flux shal cnv|XXX|2|T|  9.9900001E+20|||grid_xt,grid_yt|
|gfs_phys|sigma_shal|updraft fractional area from mass-flux shal cnv|XXX|2|T|  9.9900001E+20|||grid_xt,grid_yt|
|gfs_phys|pwat|atmos column precipitable water [kg/m**2]|kg/m**2|2|T|  9.9900001E+20|||grid_xt,grid_yt|
|gfs_phys|tmp_hyblev1|layer 1 temperature|K|2|T|  9.9900001E+20|||grid_xt,grid_yt|
|gfs_phys|spfh_hyblev1|layer 1 specific humidity|kg/kg|2|T|  9.9900001E+20|||grid_xt,grid_yt|
|gfs_phys|ugrd_hyblev1|layer 1 zonal wind|m/s|2|T|  9.9900001E+20|||grid_xt,grid_yt|
|gfs_phys|vgrd_hyblev1|layer 1 meridional wind|m/s|2|T|  9.9900001E+20|||grid_xt,grid_yt|
|gfs_phys|sfexc|Exchange Coefficient|kg/m2/s|2|T|  9.9900001E+20|||grid_xt,grid_yt|
|gfs_phys|acond|Aerodynamic conductance|m/s|2|T|  9.9900001E+20|||grid_xt,grid_yt|
|gfs_phys|DLWRFI|Instantaneous surface-temperature-adjusted downward longwave flux at the surface|w/m**2|2|T|  9.9900001E+20|||grid_xt,grid_yt|
|gfs_phys|ULWRFI|Instantaneous surface-temperature-adjusted upward longwave flux at the surface|w/m**2|2|T|  9.9900001E+20|||grid_xt,grid_yt|
|gfs_phys|DSWRFI|Instantaneous zenith-angle-adjusted downward shortwave flux at the surface|w/m**2|2|T|  9.9900001E+20|||grid_xt,grid_yt|
|gfs_phys|USWRFI|Instantaneous zenith-angle-adjusted upward shortwave flux at the surface|w/m**2|2|T|  9.9900001E+20|||grid_xt,grid_yt|
|gfs_phys|dusfci|instantaneous u component of surface stress|XXX|2|T|  9.9900001E+20|||grid_xt,grid_yt|
|gfs_phys|dvsfci|instantaneous v component of surface stress|XXX|2|T|  9.9900001E+20|||grid_xt,grid_yt|
|gfs_phys|shtfl|instantaneous surface sensible heat flux|XXX|2|T|  9.9900001E+20|||grid_xt,grid_yt|
|gfs_phys|lhtfl|instantaneous surface latent heat flux|W/m**2|2|T|  9.9900001E+20|||grid_xt,grid_yt|
|gfs_phys|gfluxi|instantaneous surface ground heat flux|W/m**2|2|T|  9.9900001E+20|||grid_xt,grid_yt|
|gfs_phys|pevpr|instantaneous surface potential evaporation|W/M**2|2|T|  9.9900001E+20|||grid_xt,grid_yt|
|gfs_phys|wilt|wiltimg point (volumetric)|Proportion|2|T|  9.9900001E+20|||grid_xt,grid_yt|
|gfs_phys|fldcp|Field Capacity (volumetric)|fraction|2|T|  9.9900001E+20|||grid_xt,grid_yt|
|gfs_phys|wet1|normalized soil wetness|XXX|2|T|  9.9900001E+20|||grid_xt,grid_yt|
|gfs_phys|cpofp|Percent frozen precipitation|%|2|T|  9.9900001E+20|||grid_xt,grid_yt|
|cosp|cltisccp|ISCCP Total Cloud Fraction / cloud_area_fraction|%|2|T| -1.0000000E+30|||grid_xt,grid_yt|
|cosp|meantbisccp|ISCCP all-sky 10.5 micron brightness temperature / toa_brightness_temperature|K|2|T| -1.0000000E+30|||grid_xt,grid_yt|
|cosp|meantbclrisccp|ISCCP clear-sky 10.5 micron brightness temperature / toa_brightness_temperature_assuming_clear_sky|K|2|T| -1.0000000E+30|||grid_xt,grid_yt|
|cosp|pctisccp|ISCCP Mean Cloud Top Pressure / air_pressure_at_cloud_top|hPa|2|T| -1.0000000E+30|||grid_xt,grid_yt|
|cosp|tauisccp|ISCCP Mean Optical Depth / atmosphere_optical_thickness_due_to_cloud|1|2|T| -1.0000000E+30|||grid_xt,grid_yt|
|cosp|albisccp|ISCCP Mean Cloud Albedo / cloud_albedo|1|2|T| -1.0000000E+30|||grid_xt,grid_yt|
|cosp|misr_meanztop|MISR Mean Cloud Top Height / cloud_top_altitude|m|2|T| -1.0000000E+30|||grid_xt,grid_yt|
|cosp|misr_cldarea|MISR cloud cover / cloud_area_fraction|1|2|T| -1.0000000E+30|||grid_xt,grid_yt|
|cosp|cltmodis|MODIS Total Cloud Fraction / cloud_area_fraction|%|2|T| -1.0000000E+30|||grid_xt,grid_yt|
|cosp|clwmodis|MODIS Liquid Cloud Fraction / cloud_area_fraction|%|2|T| -1.0000000E+30|||grid_xt,grid_yt|
|cosp|climodis|MODIS Ice Cloud Fraction / cloud_area_fraction|%|2|T| -1.0000000E+30|||grid_xt,grid_yt|
|cosp|clhmodis|MODIS High Level Cloud Fraction / cloud_area_fraction_in_atmosphere_layer|%|2|T| -1.0000000E+30|||grid_xt,grid_yt|
|cosp|clmmodis|MODIS Mid Level Cloud Fraction / cloud_area_fraction_in_atmosphere_layer|%|2|T| -1.0000000E+30|||grid_xt,grid_yt|
|cosp|cllmodis|MODIS Low Level Cloud Fraction / cloud_area_fraction_in_atmosphere_layer|%|2|T| -1.0000000E+30|||grid_xt,grid_yt|
|cosp|tautmodis|MODIS Total Cloud Optical Thickness / atmosphere_optical_thickness_due_to_cloud|1|2|T| -1.0000000E+30|||grid_xt,grid_yt|
|cosp|tauwmodis|MODIS Liquid Cloud Optical Thickness / atmosphere_optical_thickness_due_to_cloud|1|2|T| -1.0000000E+30|||grid_xt,grid_yt|
|cosp|tauimodis|MODIS Ice Cloud Optical Thickness / atmosphere_optical_thickness_due_to_cloud|1|2|T| -1.0000000E+30|||grid_xt,grid_yt|
|cosp|tautlogmodis|MODIS Total Cloud Optical Thickness (Log10 Mean) / atmosphere_optical_thickness_due_to_cloud|1|2|T| -1.0000000E+30|||grid_xt,grid_yt|
|cosp|tauwlogmodis|MODIS Liquid Cloud Optical Thickness (Log10 Mean) / atmosphere_optical_thickness_due_to_cloud|1|2|T| -1.0000000E+30|||grid_xt,grid_yt|
|cosp|tauilogmodis|MODIS Ice Cloud Optical Thickness (Log10 Mean) / atmosphere_optical_thickness_due_to_cloud|1|2|T| -1.0000000E+30|||grid_xt,grid_yt|
|cosp|reffclwmodis|MODIS Liquid Cloud Particle Size / effective_radius_of_cloud_liquid_water_particle|m|2|T| -1.0000000E+30|||grid_xt,grid_yt|
|cosp|reffclimodis|MODIS Ice Cloud Particle Size / effective_radius_of_cloud_liquid_water_particle|m|2|T| -1.0000000E+30|||grid_xt,grid_yt|
|cosp|pctmodis|MODIS Cloud Top Pressure / air_pressure_at_cloud_top|hPa|2|T| -1.0000000E+30|||grid_xt,grid_yt|
|cosp|lwpmodis|MODIS Cloud Liquid Water Path / atmosphere_cloud_liquid_water_content|kg m-2|2|T| -1.0000000E+30|||grid_xt,grid_yt|
|cosp|iwpmodis|MODIS Cloud Ice Water Path / atmosphere_mass_content_of_cloud_ice|kg m-2|2|T| -1.0000000E+30|||grid_xt,grid_yt|
|cosp|cltlidarradar|CALIPSO and CloudSat Total Cloud Fraction / cloud_area_fraction|%|2|T| -1.0000000E+30|||grid_xt,grid_yt|
|cosp|cllcalipsoice|CALIPSO Ice Low Cloud Fraction / cloud_area_fraction_in_atmosphere_layer|%|2|T| -1.0000000E+30|||grid_xt,grid_yt|
|cosp|clmcalipsoice|CALIPSO Ice Mid Cloud Fraction / cloud_area_fraction_in_atmosphere_layer|%|2|T| -1.0000000E+30|||grid_xt,grid_yt|
|cosp|clhcalipsoice|CALIPSO Ice High Cloud Fraction / cloud_area_fraction_in_atmosphere_layer|%|2|T| -1.0000000E+30|||grid_xt,grid_yt|
|cosp|cltcalipsoice|CALIPSO Ice Total Cloud Fraction / cloud_area_fraction|%|2|T| -1.0000000E+30|||grid_xt,grid_yt|
|cosp|cllcalipsoliq|CALIPSO Liquid Low Cloud Fraction / cloud_area_fraction_in_atmosphere_layer|%|2|T| -1.0000000E+30|||grid_xt,grid_yt|
|cosp|clmcalipsoliq|CALIPSO Liquid Mid Cloud Fraction / cloud_area_fraction_in_atmosphere_layer|%|2|T| -1.0000000E+30|||grid_xt,grid_yt|
|cosp|clhcalipsoliq|CALIPSO Liquid High Cloud Fraction / cloud_area_fraction_in_atmosphere_layer|%|2|T| -1.0000000E+30|||grid_xt,grid_yt|
|cosp|cltcalipsoliq|CALIPSO Liquid Total Cloud Fraction / cloud_area_fraction|%|2|T| -1.0000000E+30|||grid_xt,grid_yt|
|cosp|cllcalipsoun|CALIPSO Undefined-Phase Low Level Cloud Fraction / cloud_area_fraction_in_atmosphere_layer|%|2|T| -1.0000000E+30|||grid_xt,grid_yt|
|cosp|clmcalipsoun|CALIPSO Undefined-Phase Mid Level Cloud Fraction / cloud_area_fraction_in_atmosphere_layer|%|2|T| -1.0000000E+30|||grid_xt,grid_yt|
|cosp|clhcalipsoun|CALIPSO Undefined-Phase High Level Cloud Fraction / cloud_area_fraction_in_atmosphere_layer|%|2|T| -1.0000000E+30|||grid_xt,grid_yt|
|cosp|cltcalipsoun|CALIPSO Undefined-Phase Total Cloud Fraction / cloud_area_fraction|%|2|T| -1.0000000E+30|||grid_xt,grid_yt|
|cosp|cllcalipso|CALIPSO Low Level Cloud Fraction / cloud_area_fraction_in_atmosphere_layer|%|2|T| -1.0000000E+30|||grid_xt,grid_yt|
|cosp|clmcalipso|CALIPSO Mid Level Cloud Fraction / cloud_area_fraction_in_atmosphere_layer|%|2|T| -1.0000000E+30|||grid_xt,grid_yt|
|cosp|clhcalipso|CALIPSO High Level Cloud Fraction / cloud_area_fraction_in_atmosphere_layer|%|2|T| -1.0000000E+30|||grid_xt,grid_yt|
|cosp|cltcalipso|CALIPSO Total Cloud Fraction / cloud_area_fraction|%|2|T| -1.0000000E+30|||grid_xt,grid_yt|
|cosp|clopaquecalipso|CALIPSO Opaque Cloud Cover / opaque_cloud_cover|%|2|T| -1.0000000E+30|||grid_xt,grid_yt|
|cosp|clthincalipso|CALIPSO Thin Cloud Cover / thin_cloud_cover|%|2|T| -1.0000000E+30|||grid_xt,grid_yt|
|cosp|clzopaquecalipso|CALIPSO z_opaque Altitude / z_opaque|m|2|T| -1.0000000E+30|||grid_xt,grid_yt|
|cosp|clopaquetemp|CALIPSO Opaque Cloud Temperature / opaque_cloud_temperature|K|2|T| -1.0000000E+30|||grid_xt,grid_yt|
|cosp|clthintemp|CALIPSO Thin Cloud Temperature / thin_cloud_temperature|K|2|T| -1.0000000E+30|||grid_xt,grid_yt|
|cosp|clzopaquetemp|CALIPSO z_opaque Temperature / z_opaque_temperature|K|2|T| -1.0000000E+30|||grid_xt,grid_yt|
|cosp|clopaquemeanz|CALIPSO Opaque Cloud Altitude / opaque_cloud_altitude|m|2|T| -1.0000000E+30|||grid_xt,grid_yt|
|cosp|clthinmeanz|CALIPSO Thin Cloud Altitude / thin_cloud_altitude|m|2|T| -1.0000000E+30|||grid_xt,grid_yt|
|cosp|clthinemis|CALIPSO Thin Cloud Emissivity / thin_cloud_emissivity|1|2|T| -1.0000000E+30|||grid_xt,grid_yt|
|cosp|clopaquemeanzse|CALIPSO Opaque Cloud Altitude with respect to SE / opaque_cloud_altitude_se|m|2|T| -1.0000000E+30|||grid_xt,grid_yt|
|cosp|clthinmeanzse|CALIPSO Thin Cloud Altitude with respect to SE / thin_cloud_altitude_se|m|2|T| -1.0000000E+30|||grid_xt,grid_yt|
|cosp|clzopaquecalipsose|CALIPSO z_opaque Altitude with respect to SE / z_opaque_se|m|2|T| -1.0000000E+30|||grid_xt,grid_yt|
|cosp|cllgrLidar532|GROUND LIDAR Low Level Cloud Cover / grLidar532_low_cloud_cover|%|2|T| -1.0000000E+30|||grid_xt,grid_yt|
|cosp|clmgrLidar532|GROUND LIDAR Mid Level Cloud Cover / grLidar532_mid_cloud_cover|m|2|T| -1.0000000E+30|||grid_xt,grid_yt|
|cosp|clhgrLidar532|GROUND LIDAR High Level Cloud Cover / grLidar532_high_cloud_cover|%|2|T| -1.0000000E+30|||grid_xt,grid_yt|
|cosp|cltgrLidar532|GROUND LIDAR Total Cloud Cover / grLidar532_total_cloud_cover|%|2|T| -1.0000000E+30|||grid_xt,grid_yt|
|cosp|cllatlid|ATLID Low Level Cloud Cover / atlid_low_cloud_cover|%|2|T| -1.0000000E+30|||grid_xt,grid_yt|
|cosp|clmatlid|ATLID Mid Level Cloud Cover / atlid_mid_cloud_cover|%|2|T| -1.0000000E+30|||grid_xt,grid_yt|
|cosp|clhatlid|ATLID High Level Cloud Cover /  atlid_high_cloud_cover|%|2|T| -1.0000000E+30|||grid_xt,grid_yt|
|cosp|cltatlid|ATLID Total Cloud Cover / atlid_total_cloud_cover|%|2|T| -1.0000000E+30|||grid_xt,grid_yt|
|cosp|ptcloudsatflag0|Cloudsat precipitation cover for flag0|1|2|T| -1.0000000E+30|||grid_xt,grid_yt|
|cosp|ptcloudsatflag1|Cloudsat precipitation cover for flag1|1|2|T| -1.0000000E+30|||grid_xt,grid_yt|
|cosp|ptcloudsatflag2|Cloudsat precipitation cover for flag2|1|2|T| -1.0000000E+30|||grid_xt,grid_yt|
|cosp|ptcloudsatflag3|Cloudsat precipitation cover for flag3|1|2|T| -1.0000000E+30|||grid_xt,grid_yt|
|cosp|ptcloudsatflag4|Cloudsat precipitation cover for flag4|1|2|T| -1.0000000E+30|||grid_xt,grid_yt|
|cosp|ptcloudsatflag5|Cloudsat precipitation cover for flag5|1|2|T| -1.0000000E+30|||grid_xt,grid_yt|
|cosp|ptcloudsatflag6|Cloudsat precipitation cover for flag6|1|2|T| -1.0000000E+30|||grid_xt,grid_yt|
|cosp|ptcloudsatflag7|Cloudsat precipitation cover for flag7|1|2|T| -1.0000000E+30|||grid_xt,grid_yt|
|cosp|ptcloudsatflag8|Cloudsat precipitation cover for flag8|1|2|T| -1.0000000E+30|||grid_xt,grid_yt|
|cosp|ptcloudsatflag9|Cloudsat precipitation cover for flag9|1|2|T| -1.0000000E+30|||grid_xt,grid_yt|
|cosp|cloudsatpia|Cloudsat path integrated attenuation|1|2|T| -1.0000000E+30|||grid_xt,grid_yt|
|cosp|cloudsat_tcc|CloudSat Total Cloud Fraction / cloud_area_fraction|%|2|T| -1.0000000E+30|||grid_xt,grid_yt|
|cosp|cloudsat_tcc2|CloudSat Total Cloud Fraction (no 1km) / cloud_area_fraction|%|2|T| -1.0000000E+30|||grid_xt,grid_yt|
|cosp|npdfcld|# of Non-Precipitating Clouds / number_of_slwc_nonprecip|1|2|T| -1.0000000E+30|||grid_xt,grid_yt|
|cosp|npdfdrz|# of Drizzling Clouds / number_of_slwc_drizzle|1|2|T| -1.0000000E+30|||grid_xt,grid_yt|
|cosp|npdfrain|# of Precipitating Clouds / number_of_slwc_precip|1|2|T| -1.0000000E+30|||grid_xt,grid_yt|
|gfs_phys|dt3dt_1|temperature change due to physics 1|XXX|3|T|  9.9900001E+20|||grid_xt,grid_yt,pfull|
|gfs_phys|dt3dt_2|temperature change due to physics 2|XXX|3|T|  9.9900001E+20|||grid_xt,grid_yt,pfull|
|gfs_phys|dt3dt_3|temperature change due to physics 3|XXX|3|T|  9.9900001E+20|||grid_xt,grid_yt,pfull|
|gfs_phys|dt3dt_4|temperature change due to physics 4|XXX|3|T|  9.9900001E+20|||grid_xt,grid_yt,pfull|
|gfs_phys|dt3dt_5|temperature change due to physics 5|XXX|3|T|  9.9900001E+20|||grid_xt,grid_yt,pfull|
|gfs_phys|dt3dt_6|temperature change due to physics 6|XXX|3|T|  9.9900001E+20|||grid_xt,grid_yt,pfull|
|gfs_phys|dq3dt_1|moisture change due to physics 1|XXX|3|T|  9.9900001E+20|||grid_xt,grid_yt,pfull|
|gfs_phys|dq3dt_2|moisture change due to physics 2|XXX|3|T|  9.9900001E+20|||grid_xt,grid_yt,pfull|
|gfs_phys|dq3dt_3|moisture change due to physics 3|XXX|3|T|  9.9900001E+20|||grid_xt,grid_yt,pfull|
|gfs_phys|dq3dt_4|moisture change due to physics 4|XXX|3|T|  9.9900001E+20|||grid_xt,grid_yt,pfull|
|gfs_phys|dq3dt_5|moisture change due to physics 5|XXX|3|T|  9.9900001E+20|||grid_xt,grid_yt,pfull|
|gfs_phys|dq3dt_6|moisture change due to physics 6|XXX|3|T|  9.9900001E+20|||grid_xt,grid_yt,pfull|
|gfs_phys|dq3dt_7|moisture change due to physics 7|XXX|3|T|  9.9900001E+20|||grid_xt,grid_yt,pfull|
|gfs_phys|dq3dt_8|moisture change due to physics 8|XXX|3|T|  9.9900001E+20|||grid_xt,grid_yt,pfull|
|gfs_phys|dq3dt_9|moisture change due to physics 9|XXX|3|T|  9.9900001E+20|||grid_xt,grid_yt,pfull|
|gfs_phys|du3dt_1|u momentum change due to physics 1|XXX|3|T|  9.9900001E+20|||grid_xt,grid_yt,pfull|
|gfs_phys|du3dt_2|u momentum change due to physics 2|XXX|3|T|  9.9900001E+20|||grid_xt,grid_yt,pfull|
|gfs_phys|du3dt_3|u momentum change due to physics 3|XXX|3|T|  9.9900001E+20|||grid_xt,grid_yt,pfull|
|gfs_phys|du3dt_4|u momentum change due to physics 4|XXX|3|T|  9.9900001E+20|||grid_xt,grid_yt,pfull|
|gfs_phys|dv3dt_1|v momentum change due to physics 1|XXX|3|T|  9.9900001E+20|||grid_xt,grid_yt,pfull|
|gfs_phys|dv3dt_2|v momentum change due to physics 2|XXX|3|T|  9.9900001E+20|||grid_xt,grid_yt,pfull|
|gfs_phys|dv3dt_3|v momentum change due to physics 3|XXX|3|T|  9.9900001E+20|||grid_xt,grid_yt,pfull|
|gfs_phys|dv3dt_4|v momentum change due to physics 4|XXX|3|T|  9.9900001E+20|||grid_xt,grid_yt,pfull|
|gfs_phys|dkt_pbl|instantaneous heat diffusion coefficient|m**2/s|3|T|  9.9900001E+20|||grid_xt,grid_yt,pfull|
|gfs_phys|flux_cg|instantaneous counter-gradient heat flux in ysu|K*m/s|3|T|  9.9900001E+20|||grid_xt,grid_yt,pfull|
|gfs_phys|flux_en|instantaneous entrainment heat flux in ysu|K*m/s|3|T|  9.9900001E+20|||grid_xt,grid_yt,pfull|
|gfs_phys|wu2_shal|updraft velocity square from shallow convection|m**2/s**2|3|T|  9.9900001E+20|||grid_xt,grid_yt,pfull|
|gfs_phys|eta_shal|normalized mass flux from shallow convection|non-dim|3|T|  9.9900001E+20|||grid_xt,grid_yt,pfull|
|gfs_phys|diss_est|dissipation rate for skeb|none|3|T|  9.9900001E+20|||grid_xt,grid_yt,pfull|
|gfs_phys|skebu_wts|perturbation velocity|m/s|3|T|  9.9900001E+20|||grid_xt,grid_yt,pfull|
|gfs_phys|skebv_wts|perturbation velocity|m/s|3|T|  9.9900001E+20|||grid_xt,grid_yt,pfull|
|gfs_phys|zmtnblck|level of dividing streamline|m/s|2|T|  9.9900001E+20|||grid_xt,grid_yt|
|gfs_phys|sppt_wts|perturbation velocity|m/s|3|T|  9.9900001E+20|||grid_xt,grid_yt,pfull|
|gfs_phys|shum_wts|perturbation velocity|m/s|3|T|  9.9900001E+20|||grid_xt,grid_yt,pfull|
|gfs_sfc|alnsf|mean nir albedo with strong cosz dependency|%|2|T|  9.9900001E+20|||grid_xt,grid_yt|
|gfs_sfc|alnwf|mean nir albedo with weak cosz dependency|%|2|T|  9.9900001E+20|||grid_xt,grid_yt|
|gfs_sfc|alvsf|mean vis albedo with strong cosz dependency|%|2|T|  9.9900001E+20|||grid_xt,grid_yt|
|gfs_sfc|alvwf|mean vis albedo with weak cosz dependency|%|2|T|  9.9900001E+20|||grid_xt,grid_yt|
|gfs_sfc|canopy|canopy water (cnwat in gfs data)|%|2|T|  9.9900001E+20|||grid_xt,grid_yt|
|gfs_sfc|f10m|10-meter wind speed divided by lowest model wind speed|N/A|2|T|  9.9900001E+20|||grid_xt,grid_yt|
|gfs_sfc|facsf|fractional coverage with strong cosz dependency|XXX|2|T|  9.9900001E+20|||grid_xt,grid_yt|
|gfs_sfc|facwf|fractional coverage with weak cosz dependency|XXX|2|T|  9.9900001E+20|||grid_xt,grid_yt|
|gfs_sfc|ffhh|fh parameter from PBL scheme|XXX|2|T|  9.9900001E+20|||grid_xt,grid_yt|
|gfs_sfc|ffmm|fm parameter from PBL scheme|XXX|2|T|  9.9900001E+20|||grid_xt,grid_yt|
|gfs_sfc|uustar|uustar surface frictional wind|XXX|2|T|  9.9900001E+20|||grid_xt,grid_yt|
|gfs_sfc|slope|surface slope type|XXX|2|T|  9.9900001E+20|||grid_xt,grid_yt|
|gfs_sfc|fice|surface ice concentration (ice=1; no ice=0) [fraction]|fraction|2|T|  9.9900001E+20|||grid_xt,grid_yt|
|gfs_sfc|hice|sea ice thickness (icetk in gfs_data)|XXX|2|T|  9.9900001E+20|||grid_xt,grid_yt|
|gfs_sfc|snoalb|maximum snow albedo in fraction (salbd?? in gfs data)|XXX|2|T|  9.9900001E+20|||grid_xt,grid_yt|
|gfs_sfc|shdmax|maximum fractional coverage of green vegetation|XXX|2|T|  9.9900001E+20|||grid_xt,grid_yt|
|gfs_sfc|shdmin|minimum fractional coverage of green vegetation|XXX|2|T|  9.9900001E+20|||grid_xt,grid_yt|
|gfs_sfc|snowd|surface snow depth [m]|m|2|T|  9.9900001E+20|||grid_xt,grid_yt|
|gfs_sfc|crain|instantaneous categorical rain|number|2|T|  9.9900001E+20|||grid_xt,grid_yt|
|gfs_sfc|stype|soil type in integer 1-9|N/A|2|T|  9.9900001E+20|||grid_xt,grid_yt|
|gfs_sfc|q2m|2m specific humidity [kg/kg]|kg/kg|2|T|  9.9900001E+20|||grid_xt,grid_yt|
|gfs_sfc|t2m|2m temperature [K]|K|2|T|  9.9900001E+20|||grid_xt,grid_yt|
|gfs_sfc|tsfc|surface temperature [K]|K|2|T|  9.9900001E+20|||grid_xt,grid_yt|
|gfs_sfc|qsfc|surface specific humidity [kg/kg]|kg/kg|2|T|  9.9900001E+20|||grid_xt,grid_yt|
|gfs_sfc|tg3|deep soil temperature|K|2|T|  9.9900001E+20|||grid_xt,grid_yt|
|gfs_sfc|tisfc|surface temperature over ice fraction|K|2|T|  9.9900001E+20|||grid_xt,grid_yt|
|gfs_sfc|tprcp|total precipitation|kg/m**2|2|T|  9.9900001E+20|||grid_xt,grid_yt|
|gfs_sfc|vtype|vegetation type in integer 1-13|number|2|T|  9.9900001E+20|||grid_xt,grid_yt|
|gfs_sfc|weasd|surface snow water equivalent [kg/m**2]|kg/m**2|2|T|  9.9900001E+20|||grid_xt,grid_yt|
|gfs_sfc|HGTsfc|surface geopotential height [gpm]|gpm|2|T|  9.9900001E+20|||grid_xt,grid_yt|
|gfs_sfc|SLMSKsfc|sea-land-ice mask (0-sea, 1-land, 2-ice)|N/A|2|T|  9.9900001E+20|||grid_xt,grid_yt|
|gfs_sfc|ZORLsfc|surface roughness [m]|m|2|T|  9.9900001E+20|||grid_xt,grid_yt|
|gfs_sfc|VFRACsfc|vegetation fraction|N/A|2|T|  9.9900001E+20|||grid_xt,grid_yt|
|gfs_sfc|slc_1|liquid soil mositure at layer-1|xxx|2|T|  9.9900001E+20|||grid_xt,grid_yt|
|gfs_sfc|slc_2|liquid soil mositure at layer-2|xxx|2|T|  9.9900001E+20|||grid_xt,grid_yt|
|gfs_sfc|slc_3|liquid soil mositure at layer-3|xxx|2|T|  9.9900001E+20|||grid_xt,grid_yt|
|gfs_sfc|slc_4|liquid soil mositure at layer-4|xxx|2|T|  9.9900001E+20|||grid_xt,grid_yt|
|gfs_sfc|SOILW1|volumetric soil moisture 0-10cm [fraction]|fraction|2|T|  9.9900001E+20|||grid_xt,grid_yt|
|gfs_sfc|SOILW2|volumetric soil moisture 10-40cm [fraction]|fraction|2|T|  9.9900001E+20|||grid_xt,grid_yt|
|gfs_sfc|SOILW3|volumetric soil moisture 40-100cm [fraction]|fraction|2|T|  9.9900001E+20|||grid_xt,grid_yt|
|gfs_sfc|SOILW4|volumetric soil moisture 100-200cm [fraction]|fraction|2|T|  9.9900001E+20|||grid_xt,grid_yt|
|gfs_sfc|SOILT1|soil temperature 0-10cm [K]|K|2|T|  9.9900001E+20|||grid_xt,grid_yt|
|gfs_sfc|SOILT2|soil temperature 10-40cm [K]|K|2|T|  9.9900001E+20|||grid_xt,grid_yt|
|gfs_sfc|SOILT3|soil temperature 40-100cm [K]|K|2|T|  9.9900001E+20|||grid_xt,grid_yt|
|gfs_sfc|SOILT4|soil temperature 100-200cm [K]|K|2|T|  9.9900001E+20|||grid_xt,grid_yt|
|gfs_phys|netflxsfc|net surface heat flux [W/m**2]|W/m**2|2|T|  9.9900001E+20|||grid_xt,grid_yt|
|gfs_phys|qflux_restore|restoring flux|W/m**2|2|T|  9.9900001E+20|||grid_xt,grid_yt|
|gfs_phys|tclim_iano|climatological SST plus initial anomaly|degree C|2|T|  9.9900001E+20|||grid_xt,grid_yt|
|gfs_phys|MLD|ocean mixed layer depth|m|2|T|  9.9900001E+20|||grid_xt,grid_yt|
|gfs_phys|ps_dt|surface pressure tendency|Pa/3hr|2|T|  9.9900001E+20|||grid_xt,grid_yt|
|gfs_phys|tendency_of_air_temperature_due_to_longwave_heating|temperature tendency due to longwave radiation|K/s|3|T|  9.9900001E+20|||grid_xt,grid_yt,pfull|
|gfs_phys|tendency_of_air_temperature_due_to_shortwave_heating|temperature tendency due to shortwave radiation|K/s|3|T|  9.9900001E+20|||grid_xt,grid_yt,pfull|
|gfs_phys|tendency_of_air_temperature_due_to_turbulence|temperature tendency due to turbulence scheme|K/s|3|T|  9.9900001E+20|||grid_xt,grid_yt,pfull|
|gfs_phys|tendency_of_air_temperature_due_to_deep_convection|temperature tendency due to deep convection|K/s|3|T|  9.9900001E+20|||grid_xt,grid_yt,pfull|
|gfs_phys|tendency_of_air_temperature_due_to_shallow_convection|temperature tendency due to shallow convection|K/s|3|T|  9.9900001E+20|||grid_xt,grid_yt,pfull|
|gfs_phys|tendency_of_air_temperature_due_to_microphysics|temperature tendency due to micro-physics|K/s|3|T|  9.9900001E+20|||grid_xt,grid_yt,pfull|
|gfs_phys|tendency_of_air_temperature_due_to_dissipation_of_gravity_waves|temperature tendency due to gravity wave drag|K/s|3|T|  9.9900001E+20|||grid_xt,grid_yt,pfull|
|gfs_phys|tendency_of_air_temperature_due_to_longwave_heating_assuming_clear_sky|temperature tendency due to clear sky longwave radiation|K/s|3|T|  9.9900001E+20|||grid_xt,grid_yt,pfull|
|gfs_phys|tendency_of_air_temperature_due_to_shortwave_heating_assuming_clear_sky|temperature tendency due to clear sky shortwave radiation|K/s|3|T|  9.9900001E+20|||grid_xt,grid_yt,pfull|
|gfs_phys|vertically_integrated_tendency_of_air_temperature_due_to_longwave_heating|vertically integrated temperature tendency due to longwave radiation|W/m**2|2|T|  9.9900001E+20|||grid_xt,grid_yt|
|gfs_phys|vertically_integrated_tendency_of_air_temperature_due_to_shortwave_heating|vertically integrated temperature tendency due to shortwave radiation|W/m**2|2|T|  9.9900001E+20|||grid_xt,grid_yt|
|gfs_phys|vertically_integrated_tendency_of_air_temperature_due_to_turbulence|vertically integrated temperature tendency due to turbulence scheme|W/m**2|2|T|  9.9900001E+20|||grid_xt,grid_yt|
|gfs_phys|vertically_integrated_tendency_of_air_temperature_due_to_deep_convection|vertically integrated temperature tendency due to deep convection|W/m**2|2|T|  9.9900001E+20|||grid_xt,grid_yt|
|gfs_phys|vertically_integrated_tendency_of_air_temperature_due_to_shallow_convection|vertically integrated temperature tendency due to shallow convection|W/m**2|2|T|  9.9900001E+20|||grid_xt,grid_yt|
|gfs_phys|vertically_integrated_tendency_of_air_temperature_due_to_microphysics|vertically integrated temperature tendency due to micro-physics|W/m**2|2|T|  9.9900001E+20|||grid_xt,grid_yt|
|gfs_phys|vertically_integrated_tendency_of_air_temperature_due_to_dissipation_of_gravity_waves|vertically integrated temperature tendency due to gravity wave drag|W/m**2|2|T|  9.9900001E+20|||grid_xt,grid_yt|
|gfs_phys|vertically_integrated_tendency_of_air_temperature_due_to_longwave_heating_assuming_clear_sky|vertically integrated temperature tendency due to clear sky longwave radiation|W/m**2|2|T|  9.9900001E+20|||grid_xt,grid_yt|
|gfs_phys|vertically_integrated_tendency_of_air_temperature_due_to_shortwave_heating_assuming_clear_sky|vertically integrated temperature tendency due to clear sky shortwave radiation|W/m**2|2|T|  9.9900001E+20|||grid_xt,grid_yt|
|gfs_phys|tendency_of_specific_humidity_due_to_turbulence|water vapor tendency due to turbulence scheme|kg/kg/s|3|T|  9.9900001E+20|||grid_xt,grid_yt,pfull|
|gfs_phys|tendency_of_specific_humidity_due_to_deep_convection|water vapor tendency due to deep convection|kg/kg/s|3|T|  9.9900001E+20|||grid_xt,grid_yt,pfull|
|gfs_phys|tendency_of_specific_humidity_due_to_shallow_convection|water vapor tendency due to shallow convection|kg/kg/s|3|T|  9.9900001E+20|||grid_xt,grid_yt,pfull|
|gfs_phys|tendency_of_specific_humidity_due_to_microphysics|water vapor tendency due to microphysics|kg/kg/s|3|T|  9.9900001E+20|||grid_xt,grid_yt,pfull|
|gfs_phys|tendency_of_specific_humidity_due_to_change_in_atmosphere_mass|residual water vapor tendency|kg/kg/s|3|T|  9.9900001E+20|||grid_xt,grid_yt,pfull|
|gfs_phys|vertically_integrated_tendency_of_specific_humidity_due_to_turbulence|vertically integrated water vapor tendency due to turbulence scheme|kg/m**2/s|2|T|  9.9900001E+20|||grid_xt,grid_yt|
|gfs_phys|vertically_integrated_tendency_of_specific_humidity_due_to_deep_convection|vertically integrated water vapor tendency due to deep convection|kg/m**2/s|2|T|  9.9900001E+20|||grid_xt,grid_yt|
|gfs_phys|vertically_integrated_tendency_of_specific_humidity_due_to_shallow_convection|vertically integrated water vapor tendency due to shallow convection|kg/m**2/s|2|T|  9.9900001E+20|||grid_xt,grid_yt|
|gfs_phys|vertically_integrated_tendency_of_specific_humidity_due_to_microphysics|vertically integrated water vapor tendency due to microphysics|kg/m**2/s|2|T|  9.9900001E+20|||grid_xt,grid_yt|
|gfs_phys|vertically_integrated_tendency_of_specific_humidity_due_to_change_in_atmosphere_mass|vertically integrated residual water vapor tendency|kg/m**2/s|2|T|  9.9900001E+20|||grid_xt,grid_yt||
