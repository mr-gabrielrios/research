# Makefile created by mkmf $Id: mkmf,v 18.0.18.4.2.1.2.3 2013/08/21 15:04:06 Seth.Underwood Exp $ 

SRCROOT = /home/wenchang/tigress/CM2.5_tiger2_wy_killtc/exp/../src/

CPPDEFS = -Duse_libMPI -Duse_netCDF -DSPMD -DUSE_OCEAN_BGC -DENABLE_ODA -DINTERNAL_FILE_NML -DLAND_BND_TRACERS -DMAXDIMVALS_=5000000 -DMAXFIELDMETHODS_=200 -DMAXFIELDS_=200 -Duse_shared_pointers -DOLD_COS_SG -DOLD_PT_TO_T -DUSE_LOG_DIAG_FIELD_INFO -Duse_LARGEFILE


include /home/wenchang/tigress/CM2.5_tiger2_wy_killtc/exp/../bin/tigercpu_intelmpi_18.mk


.DEFAULT:
	-echo $@ does not exist.
all: CM2.5.x
CLUBB_3D_var.o: $(SRCROOT)atmos_param/clubb/CLUBB_3D_var.F90 parameters_model.o constants_clubb.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)atmos_param/clubb/CLUBB_3D_var.F90
CLUBB_driver_SCM.o: $(SRCROOT)atmos_param/clubb/CLUBB_driver_SCM.F90 constants.o mpp.o diag_manager.o time_manager.o fms.o field_manager.o tracer_manager.o rad_utilities.o aer_ccn_act.o aer_ccn_act_k.o ice_nucl.o alt_cloud.o clubb_core.o clubb_precision.o constants_clubb.o error_code.o grid_class.o parameter_indices.o parameters_model.o parameters_tunable.o stats_subs.o stats_variables.o T_in_K_module.o variables_prognostic_module.o sat_vapor_pres.o polysvp.o aerosol_params.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)atmos_param/clubb/CLUBB_driver_SCM.F90
FMS_coupler_util.o: $(SRCROOT)ocean_shared/generic_tracers/FMS_coupler_util.F90 mpp.o coupler_types.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)ocean_shared/generic_tracers/FMS_coupler_util.F90
FMS_ocmip2_co2calc.o: $(SRCROOT)ocean_shared/generic_tracers/FMS_ocmip2_co2calc.F90 mpp.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)ocean_shared/generic_tracers/FMS_ocmip2_co2calc.F90
MG_microp_3D.o: $(SRCROOT)atmos_param/clubb/MG_microp_3D.F90 fms.o constants.o cloud_rad.o time_manager.o cldwat2m_micro.o micro_mg.o strat_cloud_utilities.o strat_netcdf.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)atmos_param/clubb/MG_microp_3D.F90
MISR_simulator.o: $(SRCROOT)atmos_param/cosp/MISR_simulator/MISR_simulator.F /home/wenchang/tigress/CM2.5_tiger2_wy_killtc/exp/../src/atmos_param/cosp/cosp_defs.H
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -I$(SRCROOT)atmos_param/cosp	$(SRCROOT)atmos_param/cosp/MISR_simulator/MISR_simulator.F
MersenneTwister.o: $(SRCROOT)shared/random_numbers/MersenneTwister.F90
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)shared/random_numbers/MersenneTwister.F90
Skw_module.o: $(SRCROOT)atmos_param/clubb/CLUBB_core/Skw_module.F90 constants_clubb.o clubb_precision.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)atmos_param/clubb/CLUBB_core/Skw_module.F90
T_in_K_module.o: $(SRCROOT)atmos_param/clubb/CLUBB_core/T_in_K_module.F90 constants_clubb.o clubb_precision.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)atmos_param/clubb/CLUBB_core/T_in_K_module.F90
a2b_edge.o: $(SRCROOT)atmos_cubed_sphere/model/a2b_edge.F90 fv_grid_utils.o fv_arrays.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)atmos_cubed_sphere/model/a2b_edge.F90
adv_center.o: $(SRCROOT)mom5/ocean_param/gotm-4.0/util/adv_center.F90 /home/wenchang/tigress/CM2.5_tiger2_wy_killtc/exp/../src/mom5/ocean_param/gotm-4.0/include/cppdefs.h /home/wenchang/tigress/CM2.5_tiger2_wy_killtc/exp/../src/mom5/ocean_param/gotm-4.0/include/version.h util.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -I$(SRCROOT)mom5/ocean_param/gotm-4.0/include	$(SRCROOT)mom5/ocean_param/gotm-4.0/util/adv_center.F90
advance_helper_module.o: $(SRCROOT)atmos_param/clubb/CLUBB_core/advance_helper_module.F90 clubb_precision.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)atmos_param/clubb/CLUBB_core/advance_helper_module.F90
advance_sclrm_Nd_module.o: $(SRCROOT)atmos_param/clubb/advance_sclrm_Nd_module.F90 clubb_precision.o grid_class.o parameters_model.o lapack_wrap.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)atmos_param/clubb/advance_sclrm_Nd_module.F90
advance_windm_edsclrm_module.o: $(SRCROOT)atmos_param/clubb/CLUBB_core/advance_windm_edsclrm_module.F90 grid_class.o parameters_model.o parameters_tunable.o model_flags.o clubb_precision.o stats_type.o stats_variables.o clip_explicit.o error_code.o constants_clubb.o sponge_layer_damping.o lapack_wrap.o diffusion.o mean_adv.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)atmos_param/clubb/CLUBB_core/advance_windm_edsclrm_module.F90
advance_wp2_wp3_module.o: $(SRCROOT)atmos_param/clubb/CLUBB_core/advance_wp2_wp3_module.F90 grid_class.o parameters_tunable.o stats_type.o stats_variables.o constants_clubb.o model_flags.o clubb_precision.o error_code.o lapack_wrap.o fill_holes.o clip_explicit.o csr_matrix_class_3array.o gmres_wrap.o gmres_cache.o diffusion.o mean_adv.o hyper_diffusion_4th_ord.o advance_helper_module.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)atmos_param/clubb/CLUBB_core/advance_wp2_wp3_module.F90
advance_xm_wpxp_module.o: $(SRCROOT)atmos_param/clubb/CLUBB_core/advance_xm_wpxp_module.F90 parameters_tunable.o constants_clubb.o parameters_model.o grid_class.o model_flags.o mono_flux_limiter.o pdf_parameter_module.o clubb_precision.o error_code.o stats_type.o stats_variables.o sponge_layer_damping.o diffusion.o mean_adv.o clip_semi_implicit.o advance_helper_module.o lapack_wrap.o pos_definite_module.o clip_explicit.o fill_holes.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)atmos_param/clubb/CLUBB_core/advance_xm_wpxp_module.F90
advance_xp2_xpyp_module.o: $(SRCROOT)atmos_param/clubb/CLUBB_core/advance_xp2_xpyp_module.F90 constants_clubb.o model_flags.o parameters_tunable.o parameters_model.o grid_class.o clubb_precision.o clip_explicit.o stats_type.o error_code.o stats_variables.o array_index.o diffusion.o mean_adv.o advance_helper_module.o lapack_wrap.o interpolation.o fill_holes.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)atmos_param/clubb/CLUBB_core/advance_xp2_xpyp_module.F90
aer_ccn_act.o: $(SRCROOT)atmos_shared/tracer_driver/aer_ccn_act/aer_ccn_act.F90 fms.o mpp.o aer_ccn_act_k.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)atmos_shared/tracer_driver/aer_ccn_act/aer_ccn_act.F90
aer_ccn_act_k.o: $(SRCROOT)atmos_shared/tracer_driver/aer_ccn_act/aer_ccn_act_k.F90
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)atmos_shared/tracer_driver/aer_ccn_act/aer_ccn_act_k.F90
aer_in_act.o: $(SRCROOT)atmos_shared/tracer_driver/aer_ccn_act/aer_in_act.F90
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)atmos_shared/tracer_driver/aer_ccn_act/aer_in_act.F90
aerosol.o: $(SRCROOT)atmos_param/sea_esf_rad/aerosol.F90 time_manager.o diag_manager.o field_manager.o tracer_manager.o mpp.o fms.o interpolator.o mpp_io.o constants.o data_override.o rad_utilities.o mpp_domains.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)atmos_param/sea_esf_rad/aerosol.F90
aerosol_cloud.o: $(SRCROOT)atmos_param/strat_cloud/aerosol_cloud.F90 fms.o constants.o rad_utilities.o mpp.o polysvp.o aerosol_params.o aer_ccn_act.o ice_nucl.o strat_cloud_utilities.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)atmos_param/strat_cloud/aerosol_cloud.F90
aerosol_params.o: $(SRCROOT)atmos_shared/tracer_driver/aer_ccn_act/aerosol_params.F90 fms.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)atmos_shared/tracer_driver/aer_ccn_act/aerosol_params.F90
aerosolrad_package.o: $(SRCROOT)atmos_param/sea_esf_rad/aerosolrad_package.F90 mpp.o fms.o mpp_io.o time_manager.o diag_manager.o interpolator.o rad_utilities.o esfsw_parameters.o longwave_params.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)atmos_param/sea_esf_rad/aerosolrad_package.F90
affinity.o: $(SRCROOT)shared/mpp/affinity.c
	$(CC) $(CPPDEFS) $(CPPFLAGS) $(CFLAGS) $(OTHERFLAGS) $(OTHER_CFLAGS) -c	$(SRCROOT)shared/mpp/affinity.c
algebraiclength.o: $(SRCROOT)mom5/ocean_param/gotm-4.0/turbulence/algebraiclength.F90 /home/wenchang/tigress/CM2.5_tiger2_wy_killtc/exp/../src/mom5/ocean_param/gotm-4.0/include/cppdefs.h /home/wenchang/tigress/CM2.5_tiger2_wy_killtc/exp/../src/mom5/ocean_param/gotm-4.0/include/version.h turbulence.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -I$(SRCROOT)mom5/ocean_param/gotm-4.0/include	$(SRCROOT)mom5/ocean_param/gotm-4.0/turbulence/algebraiclength.F90
alpha_mnb.o: $(SRCROOT)mom5/ocean_param/gotm-4.0/turbulence/alpha_mnb.F90 /home/wenchang/tigress/CM2.5_tiger2_wy_killtc/exp/../src/mom5/ocean_param/gotm-4.0/include/cppdefs.h /home/wenchang/tigress/CM2.5_tiger2_wy_killtc/exp/../src/mom5/ocean_param/gotm-4.0/include/version.h turbulence.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -I$(SRCROOT)mom5/ocean_param/gotm-4.0/include	$(SRCROOT)mom5/ocean_param/gotm-4.0/turbulence/alpha_mnb.F90
alt_cloud.o: $(SRCROOT)atmos_param/clubb/alt_cloud.F90 sat_vapor_pres.o constants.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)atmos_param/clubb/alt_cloud.F90
amip_interp.o: $(SRCROOT)shared/amip_interp/amip_interp.F90 time_interp.o time_manager.o horiz_interp.o fms.o fms_io.o constants.o platform.o mpp.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)shared/amip_interp/amip_interp.F90
anl_erf.o: $(SRCROOT)atmos_param/clubb/CLUBB_core/anl_erf.F90
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)atmos_param/clubb/CLUBB_core/anl_erf.F90
array_index.o: $(SRCROOT)atmos_param/clubb/CLUBB_core/array_index.F90
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)atmos_param/clubb/CLUBB_core/array_index.F90
array_lib.o: $(SRCROOT)atmos_param/cosp/quickbeam/array_lib.F90 /home/wenchang/tigress/CM2.5_tiger2_wy_killtc/exp/../src/atmos_param/cosp/cosp_defs.H mrgrnk.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -I$(SRCROOT)atmos_param/cosp	$(SRCROOT)atmos_param/cosp/quickbeam/array_lib.F90
astronomy.o: $(SRCROOT)shared/astronomy/astronomy.F90 fms.o time_manager.o constants.o mpp.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)shared/astronomy/astronomy.F90
atmos_age_tracer.o: $(SRCROOT)atmos_shared/tracer_driver/atmos_age_tracer.F90 fms.o time_manager.o diag_manager.o tracer_manager.o field_manager.o interpolator.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)atmos_shared/tracer_driver/atmos_age_tracer.F90
atmos_carbon_aerosol.o: $(SRCROOT)atmos_shared/tracer_driver/atmos_carbon_aerosol.F90 mpp.o fms.o time_manager.o time_interp.o diag_manager.o tracer_manager.o field_manager.o atmos_tracer_utilities.o interpolator.o constants.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)atmos_shared/tracer_driver/atmos_carbon_aerosol.F90
atmos_ch3i.o: $(SRCROOT)atmos_shared/tracer_driver/atmos_ch3i.F90 mpp.o fms.o field_manager.o tracer_manager.o time_manager.o diag_manager.o interpolator.o constants.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)atmos_shared/tracer_driver/atmos_ch3i.F90
atmos_co2.o: $(SRCROOT)atmos_shared/tracer_driver/atmos_co2.F90 mpp.o fms.o tracer_manager.o field_manager.o diag_manager.o time_manager.o constants.o data_override.o coupler_types.o atmos_ocean_fluxes.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)atmos_shared/tracer_driver/atmos_co2.F90
atmos_convection_tracer.o: $(SRCROOT)atmos_shared/tracer_driver/atmos_convection_tracer.F90 fms.o time_manager.o diag_manager.o tracer_manager.o field_manager.o atmos_tracer_utilities.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)atmos_shared/tracer_driver/atmos_convection_tracer.F90
atmos_dust.o: $(SRCROOT)atmos_shared/tracer_driver/atmos_dust.F90 mpp.o fms.o time_manager.o diag_manager.o tracer_manager.o field_manager.o atmos_tracer_utilities.o interpolator.o constants.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)atmos_shared/tracer_driver/atmos_dust.F90
atmos_lib.o: $(SRCROOT)atmos_param/cosp/quickbeam/atmos_lib.F90 /home/wenchang/tigress/CM2.5_tiger2_wy_killtc/exp/../src/atmos_param/cosp/cosp_defs.H
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -I$(SRCROOT)atmos_param/cosp	$(SRCROOT)atmos_param/cosp/quickbeam/atmos_lib.F90
atmos_model.o: $(SRCROOT)atmos_coupled/atmos_model.F90 mpp.o mpp_domains.o fms.o fms_io.o time_manager.o field_manager.o tracer_manager.o diag_integral.o atmosphere.o xgrid.o coupler_types.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)atmos_coupled/atmos_model.F90
atmos_nudge.o: $(SRCROOT)atmos_shared/atmos_nudge/atmos_nudge.F90 mpp.o fms.o time_manager.o data_override.o diag_manager.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)atmos_shared/atmos_nudge/atmos_nudge.F90
atmos_ocean_fluxes.o: $(SRCROOT)shared/coupler/atmos_ocean_fluxes.F90 mpp.o coupler_types.o field_manager.o fm_util.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)shared/coupler/atmos_ocean_fluxes.F90
atmos_radon.o: $(SRCROOT)atmos_shared/tracer_driver/atmos_radon.F90 fms.o time_manager.o diag_manager.o tracer_manager.o field_manager.o atmos_tracer_utilities.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)atmos_shared/tracer_driver/atmos_radon.F90
atmos_sea_salt.o: $(SRCROOT)atmos_shared/tracer_driver/atmos_sea_salt.F90 mpp.o fms.o time_manager.o diag_manager.o tracer_manager.o field_manager.o constants.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)atmos_shared/tracer_driver/atmos_sea_salt.F90
atmos_soa.o: $(SRCROOT)atmos_shared/tracer_driver/atmos_soa.F90 mpp.o fms.o time_manager.o diag_manager.o tracer_manager.o field_manager.o constants.o interpolator.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)atmos_shared/tracer_driver/atmos_soa.F90
atmos_sulfate.o: $(SRCROOT)atmos_shared/tracer_driver/atmos_sulfate.F90 mpp.o fms.o time_manager.o time_interp.o diag_manager.o tracer_manager.o field_manager.o interpolator.o constants.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)atmos_shared/tracer_driver/atmos_sulfate.F90
atmos_sulfur_hex.o: $(SRCROOT)atmos_shared/tracer_driver/atmos_sulfur_hex.F90 fms.o mpp.o time_manager.o diag_manager.o tracer_manager.o field_manager.o atmos_tracer_utilities.o constants.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)atmos_shared/tracer_driver/atmos_sulfur_hex.F90
atmos_tracer_driver.o: $(SRCROOT)atmos_shared/tracer_driver/atmos_tracer_driver.F90 mpp.o fms.o time_manager.o diag_manager.o astronomy.o tracer_manager.o field_manager.o atmos_tracer_utilities.o constants.o atmos_radon.o atmos_carbon_aerosol.o atmos_convection_tracer.o atmos_sulfur_hex.o atmos_ch3i.o atmos_sea_salt.o atmos_dust.o atmos_sulfate.o atmos_soa.o tropchem_driver.o strat_chem_driver.o atmos_age_tracer.o atmos_co2.o interpolator.o coupler_types.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)atmos_shared/tracer_driver/atmos_tracer_driver.F90
atmos_tracer_utilities.o: $(SRCROOT)atmos_shared/tracer_driver/atmos_tracer_utilities.F90 fms.o time_manager.o diag_manager.o tracer_manager.o field_manager.o horiz_interp.o monin_obukhov.o constants.o interpolator.o astronomy.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)atmos_shared/tracer_driver/atmos_tracer_utilities.F90
atmosphere.o: $(SRCROOT)atmos_cubed_sphere/driver/coupled/atmosphere.F90 constants.o time_manager.o fms.o mpp.o mpp_domains.o xgrid.o diag_manager.o field_manager.o tracer_manager.o fv_arrays.o fv_control.o fv_io.o fv_dynamics.o fv_nesting.o fv_diagnostics.o fv_restart.o fv_timing.o fv_physics.o fv_mp_mod.o atmos_nudge.o fv_climate_nudge.o fv_ada_nudge.o fv_nudge.o boundary.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)atmos_cubed_sphere/driver/coupled/atmosphere.F90
axis_utils.o: $(SRCROOT)shared/axis_utils/axis_utils.F90 mpp_io.o mpp.o fms.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)shared/axis_utils/axis_utils.F90
betaDistribution.o: $(SRCROOT)atmos_param/cloud_generator/betaDistribution.F90 fms.o mpp.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)atmos_param/cloud_generator/betaDistribution.F90
betts_miller.o: $(SRCROOT)atmos_param/betts_miller/betts_miller.F90 mpp.o fms.o sat_vapor_pres.o constants.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)atmos_param/betts_miller/betts_miller.F90
bm_massflux.o: $(SRCROOT)atmos_param/betts_miller/bm_massflux.F90 mpp.o fms.o sat_vapor_pres.o constants.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)atmos_param/betts_miller/bm_massflux.F90
bm_omp.o: $(SRCROOT)atmos_param/betts_miller/bm_omp.F90 mpp.o fms.o sat_vapor_pres.o constants.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)atmos_param/betts_miller/bm_omp.F90
boundary.o: $(SRCROOT)atmos_cubed_sphere/model/boundary.F90 fv_mp_mod.o constants.o mpp_domains.o mpp.o fv_arrays.o fv_timing.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)atmos_cubed_sphere/model/boundary.F90
bulkphys_rad.o: $(SRCROOT)atmos_param/sea_esf_rad/bulkphys_rad.F90 mpp.o fms.o rad_utilities.o rh_based_clouds.o diag_clouds_W.o strat_clouds_W.o standalone_clouds.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)atmos_param/sea_esf_rad/bulkphys_rad.F90
calendar.o: $(SRCROOT)atmos_param/clubb/CLUBB_core/calendar.F90 clubb_precision.o constants_clubb.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)atmos_param/clubb/CLUBB_core/calendar.F90
cana_tile.o: $(SRCROOT)land_lad2/canopy_air/cana_tile.F90 land_tile_diag_sel.o constants.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)land_lad2/canopy_air/cana_tile.F90
canopy_air.o: $(SRCROOT)land_lad2/canopy_air/canopy_air.F90 /home/wenchang/tigress/CM2.5_tiger2_wy_killtc/exp/../src/land_lad2/canopy_air/../shared/debug.inc mpp.o fms.o fms_io.o time_manager.o constants.o sphum.o nf_utils.o land_constants.o cana_tile.o land_tile.o land_tile_diag.o land_data.o land_tile_io.o land_debug.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -I$(SRCROOT)land_lad2/canopy_air/../shared	$(SRCROOT)land_lad2/canopy_air/canopy_air.F90
cg_drag.o: $(SRCROOT)atmos_param/cg_drag/cg_drag.F90 mpp.o fms.o fms_io.o time_manager.o diag_manager.o constants.o column_diagnostics.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)atmos_param/cg_drag/cg_drag.F90
check_nan.o: $(SRCROOT)atmos_param/strat_cloud/check_nan.F90 fms.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)atmos_param/strat_cloud/check_nan.F90
cldwat2m_micro.o: $(SRCROOT)atmos_param/strat_cloud/cldwat2m_micro.F90 gamma_mg.o strat_cloud_utilities.o mpp.o fms.o simple_pdf.o sat_vapor_pres.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)atmos_param/strat_cloud/cldwat2m_micro.F90
clip_explicit.o: $(SRCROOT)atmos_param/clubb/CLUBB_core/clip_explicit.F90 grid_class.o parameters_model.o model_flags.o clubb_precision.o stats_type.o stats_variables.o constants_clubb.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)atmos_param/clubb/CLUBB_core/clip_explicit.F90
clip_semi_implicit.o: $(SRCROOT)atmos_param/clubb/CLUBB_core/clip_semi_implicit.F90 clubb_precision.o constants_clubb.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)atmos_param/clubb/CLUBB_core/clip_semi_implicit.F90
cloud_generator.o: $(SRCROOT)atmos_param/cloud_generator/cloud_generator.F90 sat_vapor_pres.o constants.o mpp.o fms.o random_numbers.o betaDistribution.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)atmos_param/cloud_generator/cloud_generator.F90
cloud_interpolator.o: $(SRCROOT)shared/drifters/cloud_interpolator.F90 /home/wenchang/tigress/CM2.5_tiger2_wy_killtc/exp/../src/shared/include/fms_platform.h
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -I$(SRCROOT)shared/include	$(SRCROOT)shared/drifters/cloud_interpolator.F90
cloud_obs.o: $(SRCROOT)atmos_param/cloud_obs/cloud_obs.F90 horiz_interp.o mpp.o fms.o fms_io.o time_manager.o time_interp.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)atmos_param/cloud_obs/cloud_obs.F90
cloud_rad.o: $(SRCROOT)atmos_param/cloud_rad/cloud_rad.F90 mpp.o fms.o constants.o gamma_mg.o mg_const.o diag_manager.o time_manager.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)atmos_param/cloud_rad/cloud_rad.F90
cloud_spec.o: $(SRCROOT)atmos_param/sea_esf_rad/cloud_spec.F90 time_manager.o mpp.o fms.o tracer_manager.o field_manager.o data_override.o random_numbers.o cloud_generator.o constants.o rad_utilities.o esfsw_parameters.o strat_clouds_W.o diag_clouds_W.o zetac_clouds_W.o specified_clouds_W.o rh_based_clouds.o donner_deep_clouds_W.o uw_clouds_W.o mgrp_prscr_clds.o standalone_clouds.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)atmos_param/sea_esf_rad/cloud_spec.F90
cloud_zonal.o: $(SRCROOT)atmos_param/cloud_zonal/cloud_zonal.F90 time_manager.o time_interp.o fms.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)atmos_param/cloud_zonal/cloud_zonal.F90
cloudrad_diagnostics.o: $(SRCROOT)atmos_param/sea_esf_rad/cloudrad_diagnostics.F90 mpp.o fms.o time_manager.o diag_manager.o constants.o rad_utilities.o esfsw_parameters.o microphys_rad.o isccp_clouds.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)atmos_param/sea_esf_rad/cloudrad_diagnostics.F90
cloudrad_package.o: $(SRCROOT)atmos_param/sea_esf_rad/cloudrad_package.F90 mpp.o fms.o time_manager.o rad_utilities.o esfsw_parameters.o cloudrad_diagnostics.o bulkphys_rad.o microphys_rad.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)atmos_param/sea_esf_rad/cloudrad_package.F90
clouds.o: $(SRCROOT)atmos_param/clouds/clouds.F90 cloud_rad.o cloud_zonal.o cloud_obs.o time_manager.o mpp.o fms.o rh_clouds.o strat_cloud.o diag_cloud.o diag_manager.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)atmos_param/clouds/clouds.F90
clubb_core.o: $(SRCROOT)atmos_param/clubb/CLUBB_core/clubb_core.F90 constants_clubb.o parameters_tunable.o parameters_model.o model_flags.o grid_class.o numerical_check.o variables_diagnostic_module.o pdf_parameter_module.o advance_sclrm_Nd_module.o advance_xm_wpxp_module.o advance_xp2_xpyp_module.o surface_varnce_module.o pdf_closure_module.o mixing_length.o advance_windm_edsclrm_module.o saturation.o advance_wp2_wp3_module.o clubb_precision.o error_code.o Skw_module.o clip_explicit.o T_in_K_module.o stats_subs.o stats_type.o stats_variables.o fill_holes.o sigma_sqd_w_module.o parameter_indices.o variables_prognostic_module.o csr_matrix_class_3array.o gmres_wrap.o gmres_cache.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)atmos_param/clubb/CLUBB_core/clubb_core.F90
clubb_precision.o: $(SRCROOT)atmos_param/clubb/CLUBB_core/clubb_precision.F90
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)atmos_param/clubb/CLUBB_core/clubb_precision.F90
cmue_a.o: $(SRCROOT)mom5/ocean_param/gotm-4.0/turbulence/cmue_a.F90 /home/wenchang/tigress/CM2.5_tiger2_wy_killtc/exp/../src/mom5/ocean_param/gotm-4.0/include/cppdefs.h /home/wenchang/tigress/CM2.5_tiger2_wy_killtc/exp/../src/mom5/ocean_param/gotm-4.0/include/version.h turbulence.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -I$(SRCROOT)mom5/ocean_param/gotm-4.0/include	$(SRCROOT)mom5/ocean_param/gotm-4.0/turbulence/cmue_a.F90
cmue_b.o: $(SRCROOT)mom5/ocean_param/gotm-4.0/turbulence/cmue_b.F90 /home/wenchang/tigress/CM2.5_tiger2_wy_killtc/exp/../src/mom5/ocean_param/gotm-4.0/include/cppdefs.h /home/wenchang/tigress/CM2.5_tiger2_wy_killtc/exp/../src/mom5/ocean_param/gotm-4.0/include/version.h turbulence.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -I$(SRCROOT)mom5/ocean_param/gotm-4.0/include	$(SRCROOT)mom5/ocean_param/gotm-4.0/turbulence/cmue_b.F90
cmue_c.o: $(SRCROOT)mom5/ocean_param/gotm-4.0/turbulence/cmue_c.F90 /home/wenchang/tigress/CM2.5_tiger2_wy_killtc/exp/../src/mom5/ocean_param/gotm-4.0/include/cppdefs.h /home/wenchang/tigress/CM2.5_tiger2_wy_killtc/exp/../src/mom5/ocean_param/gotm-4.0/include/version.h turbulence.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -I$(SRCROOT)mom5/ocean_param/gotm-4.0/include	$(SRCROOT)mom5/ocean_param/gotm-4.0/turbulence/cmue_c.F90
cmue_d.o: $(SRCROOT)mom5/ocean_param/gotm-4.0/turbulence/cmue_d.F90 /home/wenchang/tigress/CM2.5_tiger2_wy_killtc/exp/../src/mom5/ocean_param/gotm-4.0/include/cppdefs.h /home/wenchang/tigress/CM2.5_tiger2_wy_killtc/exp/../src/mom5/ocean_param/gotm-4.0/include/version.h turbulence.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -I$(SRCROOT)mom5/ocean_param/gotm-4.0/include	$(SRCROOT)mom5/ocean_param/gotm-4.0/turbulence/cmue_d.F90
cmue_ma.o: $(SRCROOT)mom5/ocean_param/gotm-4.0/turbulence/cmue_ma.F90 /home/wenchang/tigress/CM2.5_tiger2_wy_killtc/exp/../src/mom5/ocean_param/gotm-4.0/include/cppdefs.h /home/wenchang/tigress/CM2.5_tiger2_wy_killtc/exp/../src/mom5/ocean_param/gotm-4.0/include/version.h turbulence.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -I$(SRCROOT)mom5/ocean_param/gotm-4.0/include	$(SRCROOT)mom5/ocean_param/gotm-4.0/turbulence/cmue_ma.F90
cmue_rf.o: $(SRCROOT)mom5/ocean_param/gotm-4.0/turbulence/cmue_rf.F90 /home/wenchang/tigress/CM2.5_tiger2_wy_killtc/exp/../src/mom5/ocean_param/gotm-4.0/include/cppdefs.h /home/wenchang/tigress/CM2.5_tiger2_wy_killtc/exp/../src/mom5/ocean_param/gotm-4.0/include/version.h turbulence.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -I$(SRCROOT)mom5/ocean_param/gotm-4.0/include	$(SRCROOT)mom5/ocean_param/gotm-4.0/turbulence/cmue_rf.F90
cmue_sg.o: $(SRCROOT)mom5/ocean_param/gotm-4.0/turbulence/cmue_sg.F90 /home/wenchang/tigress/CM2.5_tiger2_wy_killtc/exp/../src/mom5/ocean_param/gotm-4.0/include/cppdefs.h /home/wenchang/tigress/CM2.5_tiger2_wy_killtc/exp/../src/mom5/ocean_param/gotm-4.0/include/version.h turbulence.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -I$(SRCROOT)mom5/ocean_param/gotm-4.0/include	$(SRCROOT)mom5/ocean_param/gotm-4.0/turbulence/cmue_sg.F90
co2_data.o: $(SRCROOT)atmos_param/fsrad/co2_data.F90 fs_profile.o co2int.o fms.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)atmos_param/fsrad/co2_data.F90
co2int.o: $(SRCROOT)atmos_param/fsrad/co2int.F90 fs_profile.o fms.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)atmos_param/fsrad/co2int.F90
column_diagnostics.o: $(SRCROOT)shared/column_diagnostics/column_diagnostics.F90 mpp_io.o fms.o time_manager.o constants.o mpp.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)shared/column_diagnostics/column_diagnostics.F90
compute_cpsi3.o: $(SRCROOT)mom5/ocean_param/gotm-4.0/turbulence/compute_cpsi3.F90 /home/wenchang/tigress/CM2.5_tiger2_wy_killtc/exp/../src/mom5/ocean_param/gotm-4.0/include/cppdefs.h /home/wenchang/tigress/CM2.5_tiger2_wy_killtc/exp/../src/mom5/ocean_param/gotm-4.0/include/version.h turbulence.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -I$(SRCROOT)mom5/ocean_param/gotm-4.0/include	$(SRCROOT)mom5/ocean_param/gotm-4.0/turbulence/compute_cpsi3.F90
compute_rist.o: $(SRCROOT)mom5/ocean_param/gotm-4.0/turbulence/compute_rist.F90 /home/wenchang/tigress/CM2.5_tiger2_wy_killtc/exp/../src/mom5/ocean_param/gotm-4.0/include/cppdefs.h /home/wenchang/tigress/CM2.5_tiger2_wy_killtc/exp/../src/mom5/ocean_param/gotm-4.0/include/version.h turbulence.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -I$(SRCROOT)mom5/ocean_param/gotm-4.0/include	$(SRCROOT)mom5/ocean_param/gotm-4.0/turbulence/compute_rist.F90
constants.o: $(SRCROOT)shared/constants/constants.F90
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)shared/constants/constants.F90
constants_clubb.o: $(SRCROOT)atmos_param/clubb/CLUBB_core/constants_clubb.F90 clubb_precision.o constants.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)atmos_param/clubb/CLUBB_core/constants_clubb.F90
conv_closures.o: $(SRCROOT)atmos_param/shallow_cu/conv_closures.F90 fms.o conv_utilities_k.o conv_plumes_k.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)atmos_param/shallow_cu/conv_closures.F90
conv_plumes_k.o: $(SRCROOT)atmos_param/shallow_cu/conv_plumes_k.F90 /home/wenchang/tigress/CM2.5_tiger2_wy_killtc/exp/../src/shared/include/fms_platform.h aer_ccn_act_k.o conv_utilities_k.o sat_vapor_pres_k.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -I$(SRCROOT)shared/include	$(SRCROOT)atmos_param/shallow_cu/conv_plumes_k.F90
conv_utilities.o: $(SRCROOT)atmos_param/shallow_cu/conv_utilities.F90 constants.o fms.o conv_utilities_k.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)atmos_param/shallow_cu/conv_utilities.F90
conv_utilities_k.o: $(SRCROOT)atmos_param/shallow_cu/conv_utilities_k.F90 /home/wenchang/tigress/CM2.5_tiger2_wy_killtc/exp/../src/shared/include/fms_platform.h sat_vapor_pres_k.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -I$(SRCROOT)shared/include	$(SRCROOT)atmos_param/shallow_cu/conv_utilities_k.F90
convert_fluxes.o: $(SRCROOT)mom5/ocean_param/gotm-4.0/util/convert_fluxes.F90 /home/wenchang/tigress/CM2.5_tiger2_wy_killtc/exp/../src/mom5/ocean_param/gotm-4.0/include/cppdefs.h /home/wenchang/tigress/CM2.5_tiger2_wy_killtc/exp/../src/mom5/ocean_param/gotm-4.0/include/version.h eqstate.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -I$(SRCROOT)mom5/ocean_param/gotm-4.0/include	$(SRCROOT)mom5/ocean_param/gotm-4.0/util/convert_fluxes.F90
corr_matrix_module.o: $(SRCROOT)atmos_param/clubb/CLUBB_core/corr_matrix_module.F90 input_reader.o matrix_operations.o constants_clubb.o clubb_precision.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)atmos_param/clubb/CLUBB_core/corr_matrix_module.F90
cosp.o: $(SRCROOT)atmos_param/cosp/cosp.F90 /home/wenchang/tigress/CM2.5_tiger2_wy_killtc/exp/../src/atmos_param/cosp/cosp_defs.H cosp_types.o cosp_simulator.o cosp_modis_simulator.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -I$(SRCROOT)atmos_param/cosp	$(SRCROOT)atmos_param/cosp/cosp.F90
cosp_constants.o: $(SRCROOT)atmos_param/cosp/cosp_constants.F90 /home/wenchang/tigress/CM2.5_tiger2_wy_killtc/exp/../src/atmos_param/cosp/cosp_defs.H
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -I$(SRCROOT)atmos_param/cosp	$(SRCROOT)atmos_param/cosp/cosp_constants.F90
cosp_diagnostics.o: $(SRCROOT)atmos_param/cosp/cosp_diagnostics.F90 /home/wenchang/tigress/CM2.5_tiger2_wy_killtc/exp/../src/atmos_param/cosp/cosp_defs.H mpp.o fms.o time_manager.o diag_grid.o diag_manager.o cosp_types.o cosp_io.o cosp_constants.o lmd_ipsl_stats.o cosp_modis_simulator.o modis_simulator.o cosp_utils.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -I$(SRCROOT)atmos_param/cosp	$(SRCROOT)atmos_param/cosp/cosp_diagnostics.F90
cosp_driver.o: $(SRCROOT)atmos_param/cosp/cosp_driver.F90 /home/wenchang/tigress/CM2.5_tiger2_wy_killtc/exp/../src/atmos_param/cosp/cosp_defs.H mpp.o fms.o sat_vapor_pres.o time_manager.o cosp_types.o cosp.o cosp_io.o cosp_constants.o cosp_modis_simulator.o cosp_diagnostics.o cosp_utils.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -I$(SRCROOT)atmos_param/cosp	$(SRCROOT)atmos_param/cosp/cosp_driver.F90
cosp_io.o: $(SRCROOT)atmos_param/cosp/cosp_io.F90 /home/wenchang/tigress/CM2.5_tiger2_wy_killtc/exp/../src/atmos_param/cosp/cosp_defs.H cosp_constants.o cosp_types.o cosp_modis_simulator.o mpp.o fms.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -I$(SRCROOT)atmos_param/cosp	$(SRCROOT)atmos_param/cosp/cosp_io.F90
cosp_isccp_simulator.o: $(SRCROOT)atmos_param/cosp/cosp_isccp_simulator.F90 /home/wenchang/tigress/CM2.5_tiger2_wy_killtc/exp/../src/atmos_param/cosp/cosp_defs.H cosp_constants.o cosp_types.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -I$(SRCROOT)atmos_param/cosp	$(SRCROOT)atmos_param/cosp/cosp_isccp_simulator.F90
cosp_lidar.o: $(SRCROOT)atmos_param/cosp/cosp_lidar.F90 /home/wenchang/tigress/CM2.5_tiger2_wy_killtc/exp/../src/atmos_param/cosp/cosp_defs.H cosp_constants.o cosp_types.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -I$(SRCROOT)atmos_param/cosp	$(SRCROOT)atmos_param/cosp/cosp_lidar.F90
cosp_misr_simulator.o: $(SRCROOT)atmos_param/cosp/cosp_misr_simulator.F90 /home/wenchang/tigress/CM2.5_tiger2_wy_killtc/exp/../src/atmos_param/cosp/cosp_defs.H cosp_constants.o cosp_types.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -I$(SRCROOT)atmos_param/cosp	$(SRCROOT)atmos_param/cosp/cosp_misr_simulator.F90
cosp_modis_simulator.o: $(SRCROOT)atmos_param/cosp/cosp_modis_simulator.F90 /home/wenchang/tigress/CM2.5_tiger2_wy_killtc/exp/../src/atmos_param/cosp/cosp_defs.H cosp_constants.o cosp_types.o modis_simulator.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -I$(SRCROOT)atmos_param/cosp	$(SRCROOT)atmos_param/cosp/cosp_modis_simulator.F90
cosp_radar.o: $(SRCROOT)atmos_param/cosp/llnl/cosp_radar.F90 /home/wenchang/tigress/CM2.5_tiger2_wy_killtc/exp/../src/atmos_param/cosp/cosp_defs.H cosp_constants.o cosp_types.o cosp_utils.o radar_simulator_types.o array_lib.o atmos_lib.o format_input.o mrgrnk.o math_lib.o optics_lib.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -I$(SRCROOT)atmos_param/cosp	$(SRCROOT)atmos_param/cosp/llnl/cosp_radar.F90
cosp_rttov_simulator.o: $(SRCROOT)atmos_param/cosp/cosp_rttov_simulator.F90 /home/wenchang/tigress/CM2.5_tiger2_wy_killtc/exp/../src/atmos_param/cosp/cosp_defs.H cosp_constants.o cosp_types.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -I$(SRCROOT)atmos_param/cosp	$(SRCROOT)atmos_param/cosp/cosp_rttov_simulator.F90
cosp_simulator.o: $(SRCROOT)atmos_param/cosp/cosp_simulator.F90 /home/wenchang/tigress/CM2.5_tiger2_wy_killtc/exp/../src/atmos_param/cosp/cosp_defs.H cosp_types.o cosp_radar.o cosp_lidar.o cosp_isccp_simulator.o cosp_modis_simulator.o cosp_misr_simulator.o cosp_rttov_simulator.o cosp_stats.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -I$(SRCROOT)atmos_param/cosp	$(SRCROOT)atmos_param/cosp/cosp_simulator.F90
cosp_stats.o: $(SRCROOT)atmos_param/cosp/cosp_stats.F90 /home/wenchang/tigress/CM2.5_tiger2_wy_killtc/exp/../src/atmos_param/cosp/cosp_defs.H cosp_constants.o cosp_types.o llnl_stats.o lmd_ipsl_stats.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -I$(SRCROOT)atmos_param/cosp	$(SRCROOT)atmos_param/cosp/cosp_stats.F90
cosp_types.o: $(SRCROOT)atmos_param/cosp/cosp_types.F90 /home/wenchang/tigress/CM2.5_tiger2_wy_killtc/exp/../src/atmos_param/cosp/cosp_defs.H cosp_constants.o cosp_utils.o radar_simulator_types.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -I$(SRCROOT)atmos_param/cosp	$(SRCROOT)atmos_param/cosp/cosp_types.F90
cosp_utils.o: $(SRCROOT)atmos_param/cosp/cosp_utils.F90 /home/wenchang/tigress/CM2.5_tiger2_wy_killtc/exp/../src/atmos_param/cosp/cosp_defs.H cosp_constants.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -I$(SRCROOT)atmos_param/cosp	$(SRCROOT)atmos_param/cosp/cosp_utils.F90
coupler_main.o: $(SRCROOT)coupler/coupler_main.F90 constants.o time_manager.o fms.o fms_io.o diag_manager.o field_manager.o tracer_manager.o coupler_types.o data_override.o atmos_model.o land_model.o ice_model.o ocean_model.o flux_exchange.o atmos_tracer_driver.o mpp.o mpp_io.o mpp_domains.o memutils.o ensemble_manager.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)coupler/coupler_main.F90
coupler_types.o: $(SRCROOT)shared/coupler/coupler_types.F90 field_manager.o mpp.o fm_util.o time_manager.o diag_manager.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)shared/coupler/coupler_types.F90
create_xgrid.o: $(SRCROOT)shared/mosaic/create_xgrid.c /home/wenchang/tigress/CM2.5_tiger2_wy_killtc/exp/../src/shared/mosaic/mosaic_util.h /home/wenchang/tigress/CM2.5_tiger2_wy_killtc/exp/../src/shared/mosaic/create_xgrid.h /home/wenchang/tigress/CM2.5_tiger2_wy_killtc/exp/../src/shared/mosaic/constant.h
	$(CC) $(CPPDEFS) $(CPPFLAGS) $(CFLAGS) $(OTHERFLAGS) $(OTHER_CFLAGS) -c -I$(SRCROOT)shared/mosaic	$(SRCROOT)shared/mosaic/create_xgrid.c
csr_matrix_class_3array.o: $(SRCROOT)atmos_param/clubb/CLUBB_core/csr_matrix_class_3array.F90 constants_clubb.o grid_class.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)atmos_param/clubb/CLUBB_core/csr_matrix_class_3array.F90
cu_mo_trans.o: $(SRCROOT)atmos_param/cu_mo_trans/cu_mo_trans.F90 constants.o mpp.o fms.o diag_manager.o time_manager.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)atmos_param/cu_mo_trans/cu_mo_trans.F90
cumulus_closure_k.o: $(SRCROOT)atmos_param/donner_deep/cumulus_closure_k.F90 donner_types.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)atmos_param/donner_deep/cumulus_closure_k.F90
damping_driver.o: $(SRCROOT)atmos_param/damping_driver/damping_driver.F90 mg_drag.o cg_drag.o topo_drag.o mpp.o fms.o diag_manager.o time_manager.o constants.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)atmos_param/damping_driver/damping_driver.F90
data_override.o: $(SRCROOT)shared/data_override/data_override.F90 /home/wenchang/tigress/CM2.5_tiger2_wy_killtc/exp/../src/shared/include/fms_platform.h platform.o constants.o mpp_io.o mpp.o horiz_interp.o time_interp_external.o fms_io.o fms.o axis_utils.o mpp_domains.o time_manager.o diag_manager.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -I$(SRCROOT)shared/include	$(SRCROOT)shared/data_override/data_override.F90
deep_conv.o: $(SRCROOT)atmos_param/shallow_cu/deep_conv.F90 fms.o constants.o conv_utilities.o conv_utilities_k.o conv_plumes_k.o conv_closures.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)atmos_param/shallow_cu/deep_conv.F90
detr_ice_num.o: $(SRCROOT)atmos_param/moist_processes/detr_ice_num.F90 mg_const.o mpp.o fms.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)atmos_param/moist_processes/detr_ice_num.F90
diag_axis.o: $(SRCROOT)shared/diag_manager/diag_axis.F90 mpp_domains.o fms.o diag_data.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)shared/diag_manager/diag_axis.F90
diag_cloud.o: $(SRCROOT)atmos_param/diag_cloud/diag_cloud.F90 mpp.o fms.o fms_io.o constants.o time_manager.o cloud_zonal.o diag_cloud_rad.o sat_vapor_pres.o shallow_conv.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)atmos_param/diag_cloud/diag_cloud.F90
diag_cloud_rad.o: $(SRCROOT)atmos_param/diag_cloud_rad/diag_cloud_rad.F90 mpp.o fms.o cloud_rad.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)atmos_param/diag_cloud_rad/diag_cloud_rad.F90
diag_clouds_W.o: $(SRCROOT)atmos_param/sea_esf_rad/diag_clouds_W.F90 time_manager.o diag_cloud.o diag_cloud_rad.o mpp.o fms.o rad_utilities.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)atmos_param/sea_esf_rad/diag_clouds_W.F90
diag_data.o: $(SRCROOT)shared/diag_manager/diag_data.F90 /home/wenchang/tigress/CM2.5_tiger2_wy_killtc/exp/../src/shared/include/fms_platform.h time_manager.o mpp_domains.o mpp_io.o fms.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -I$(SRCROOT)shared/include	$(SRCROOT)shared/diag_manager/diag_data.F90
diag_grid.o: $(SRCROOT)shared/diag_manager/diag_grid.F90 /home/wenchang/tigress/CM2.5_tiger2_wy_killtc/exp/../src/shared/include/fms_platform.h constants.o fms.o mpp.o mpp_domains.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -I$(SRCROOT)shared/include	$(SRCROOT)shared/diag_manager/diag_grid.F90
diag_integral.o: $(SRCROOT)atmos_param/diag_integral/diag_integral.F90 /home/wenchang/tigress/CM2.5_tiger2_wy_killtc/exp/../src/shared/include/fms_platform.h time_manager.o mpp.o fms.o constants.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -I$(SRCROOT)shared/include	$(SRCROOT)atmos_param/diag_integral/diag_integral.F90
diag_manager.o: $(SRCROOT)shared/diag_manager/diag_manager.F90 /home/wenchang/tigress/CM2.5_tiger2_wy_killtc/exp/../src/shared/include/fms_platform.h time_manager.o mpp_io.o mpp.o fms.o diag_axis.o diag_util.o diag_data.o diag_table.o diag_output.o diag_grid.o constants.o mpp_domains.o fms_io.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -I$(SRCROOT)shared/include	$(SRCROOT)shared/diag_manager/diag_manager.F90
diag_output.o: $(SRCROOT)shared/diag_manager/diag_output.F90 /home/wenchang/tigress/CM2.5_tiger2_wy_killtc/exp/../src/shared/include/fms_platform.h mpp_io.o mpp_domains.o mpp.o diag_axis.o diag_data.o time_manager.o fms.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -I$(SRCROOT)shared/include	$(SRCROOT)shared/diag_manager/diag_output.F90
diag_table.o: $(SRCROOT)shared/diag_manager/diag_table.F90 mpp_io.o mpp.o fms.o time_manager.o constants.o diag_data.o diag_util.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)shared/diag_manager/diag_table.F90
diag_util.o: $(SRCROOT)shared/diag_manager/diag_util.F90 /home/wenchang/tigress/CM2.5_tiger2_wy_killtc/exp/../src/shared/include/fms_platform.h diag_data.o diag_axis.o diag_output.o diag_grid.o fms.o fms_io.o mpp_domains.o time_manager.o mpp_io.o mpp.o constants.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -I$(SRCROOT)shared/include	$(SRCROOT)shared/diag_manager/diag_util.F90
diff_center.o: $(SRCROOT)mom5/ocean_param/gotm-4.0/util/diff_center.F90 /home/wenchang/tigress/CM2.5_tiger2_wy_killtc/exp/../src/mom5/ocean_param/gotm-4.0/include/cppdefs.h /home/wenchang/tigress/CM2.5_tiger2_wy_killtc/exp/../src/mom5/ocean_param/gotm-4.0/include/version.h util.o tridiagonal_gotm.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -I$(SRCROOT)mom5/ocean_param/gotm-4.0/include	$(SRCROOT)mom5/ocean_param/gotm-4.0/util/diff_center.F90
diff_face.o: $(SRCROOT)mom5/ocean_param/gotm-4.0/util/diff_face.F90 /home/wenchang/tigress/CM2.5_tiger2_wy_killtc/exp/../src/mom5/ocean_param/gotm-4.0/include/cppdefs.h /home/wenchang/tigress/CM2.5_tiger2_wy_killtc/exp/../src/mom5/ocean_param/gotm-4.0/include/version.h util.o tridiagonal_gotm.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -I$(SRCROOT)mom5/ocean_param/gotm-4.0/include	$(SRCROOT)mom5/ocean_param/gotm-4.0/util/diff_face.F90
diffusion.o: $(SRCROOT)atmos_param/clubb/CLUBB_core/diffusion.F90 grid_class.o clubb_precision.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)atmos_param/clubb/CLUBB_core/diffusion.F90
diffusivity.o: $(SRCROOT)atmos_param/diffusivity/diffusivity.F90 constants.o mpp.o fms.o monin_obukhov.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)atmos_param/diffusivity/diffusivity.F90
dissipationeq.o: $(SRCROOT)mom5/ocean_param/gotm-4.0/turbulence/dissipationeq.F90 /home/wenchang/tigress/CM2.5_tiger2_wy_killtc/exp/../src/mom5/ocean_param/gotm-4.0/include/cppdefs.h /home/wenchang/tigress/CM2.5_tiger2_wy_killtc/exp/../src/mom5/ocean_param/gotm-4.0/include/version.h turbulence.o util.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -I$(SRCROOT)mom5/ocean_param/gotm-4.0/include	$(SRCROOT)mom5/ocean_param/gotm-4.0/turbulence/dissipationeq.F90
donner_cape_k.o: $(SRCROOT)atmos_param/donner_deep/donner_cape_k.F90 donner_types.o sat_vapor_pres_k.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)atmos_param/donner_deep/donner_cape_k.F90
donner_cloud_model_k.o: $(SRCROOT)atmos_param/donner_deep/donner_cloud_model_k.F90 donner_types.o sat_vapor_pres_k.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)atmos_param/donner_deep/donner_cloud_model_k.F90
donner_deep.o: $(SRCROOT)atmos_param/donner_deep/donner_deep.F90 donner_types.o conv_utilities_k.o conv_plumes_k.o fms_donner.o nonfms_donner.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)atmos_param/donner_deep/donner_deep.F90
donner_deep_clouds_W.o: $(SRCROOT)atmos_param/sea_esf_rad/donner_deep_clouds_W.F90 time_manager.o mpp.o fms.o rad_utilities.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)atmos_param/sea_esf_rad/donner_deep_clouds_W.F90
donner_deep_k.o: $(SRCROOT)atmos_param/donner_deep/donner_deep_k.F90 donner_types.o conv_utilities_k.o conv_plumes_k.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)atmos_param/donner_deep/donner_deep_k.F90
donner_deep_miz.o: $(SRCROOT)atmos_param/donner_deep/donner_deep_miz.F90
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)atmos_param/donner_deep/donner_deep_miz.F90
donner_lite_k.o: $(SRCROOT)atmos_param/donner_deep/donner_lite_k.F90 donner_types.o conv_utilities_k.o conv_plumes_k.o sat_vapor_pres_k.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)atmos_param/donner_deep/donner_lite_k.F90
donner_lscloud_k.o: $(SRCROOT)atmos_param/donner_deep/donner_lscloud_k.F90 donner_types.o sat_vapor_pres_k.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)atmos_param/donner_deep/donner_lscloud_k.F90
donner_meso_k.o: $(SRCROOT)atmos_param/donner_deep/donner_meso_k.F90 donner_types.o sat_vapor_pres_k.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)atmos_param/donner_deep/donner_meso_k.F90
donner_rad_k.o: $(SRCROOT)atmos_param/donner_deep/donner_rad_k.F90 donner_types.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)atmos_param/donner_deep/donner_rad_k.F90
donner_types.o: $(SRCROOT)atmos_param/donner_deep/donner_types.F90 /home/wenchang/tigress/CM2.5_tiger2_wy_killtc/exp/../src/atmos_param/donner_deep/donner_types.h
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -I$(SRCROOT)atmos_param/donner_deep	$(SRCROOT)atmos_param/donner_deep/donner_types.F90
donner_utilities_k.o: $(SRCROOT)atmos_param/donner_deep/donner_utilities_k.F90 donner_types.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)atmos_param/donner_deep/donner_utilities_k.F90
drifters.o: $(SRCROOT)shared/drifters/drifters.F90 /home/wenchang/tigress/CM2.5_tiger2_wy_killtc/exp/../src/shared/include/fms_platform.h /home/wenchang/tigress/CM2.5_tiger2_wy_killtc/exp/../src/shared/drifters/fms_switches.h /home/wenchang/tigress/CM2.5_tiger2_wy_killtc/exp/../src/shared/drifters/drifters_push.h /home/wenchang/tigress/CM2.5_tiger2_wy_killtc/exp/../src/shared/drifters/drifters_set_field.h /home/wenchang/tigress/CM2.5_tiger2_wy_killtc/exp/../src/shared/drifters/drifters_compute_k.h mpp.o mpp_domains.o drifters_core.o drifters_input.o drifters_io.o drifters_comm.o cloud_interpolator.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -I$(SRCROOT)shared/include -I$(SRCROOT)shared/drifters	$(SRCROOT)shared/drifters/drifters.F90
drifters_comm.o: $(SRCROOT)shared/drifters/drifters_comm.F90 /home/wenchang/tigress/CM2.5_tiger2_wy_killtc/exp/../src/shared/include/fms_platform.h /home/wenchang/tigress/CM2.5_tiger2_wy_killtc/exp/../src/shared/drifters/fms_switches.h mpp.o mpp_domains.o drifters_core.o drifters_input.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -I$(SRCROOT)shared/include -I$(SRCROOT)shared/drifters	$(SRCROOT)shared/drifters/drifters_comm.F90
drifters_core.o: $(SRCROOT)shared/drifters/drifters_core.F90 /home/wenchang/tigress/CM2.5_tiger2_wy_killtc/exp/../src/shared/include/fms_platform.h
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -I$(SRCROOT)shared/include	$(SRCROOT)shared/drifters/drifters_core.F90
drifters_input.o: $(SRCROOT)shared/drifters/drifters_input.F90 /home/wenchang/tigress/CM2.5_tiger2_wy_killtc/exp/../src/shared/include/fms_platform.h
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -I$(SRCROOT)shared/include	$(SRCROOT)shared/drifters/drifters_input.F90
drifters_io.o: $(SRCROOT)shared/drifters/drifters_io.F90
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)shared/drifters/drifters_io.F90
dry_adj.o: $(SRCROOT)atmos_param/dry_adj/dry_adj.F90 mpp.o fms.o constants.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)atmos_param/dry_adj/dry_adj.F90
dsd.o: $(SRCROOT)atmos_param/cosp/quickbeam/dsd.F90 /home/wenchang/tigress/CM2.5_tiger2_wy_killtc/exp/../src/atmos_param/cosp/cosp_defs.H array_lib.o math_lib.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -I$(SRCROOT)atmos_param/cosp	$(SRCROOT)atmos_param/cosp/quickbeam/dsd.F90
dyn_core.o: $(SRCROOT)atmos_cubed_sphere/model/dyn_core.F90 constants.o mpp.o mpp_domains.o mpp_parameter.o fv_mp_mod.o sw_core.o a2b_edge.o nh_core.o tp_core.o fv_timing.o fv_diagnostics.o fv_ada_nudge.o fv_nudge.o diag_manager.o fv_arrays.o boundary.o test_cases.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)atmos_cubed_sphere/model/dyn_core.F90
edt.o: $(SRCROOT)atmos_param/edt/edt.F90 constants.o mpp.o fms.o fms_io.o diag_manager.o time_manager.o monin_obukhov.o sat_vapor_pres.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)atmos_param/edt/edt.F90
endian.o: $(SRCROOT)atmos_param/clubb/CLUBB_core/endian.F90
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)atmos_param/clubb/CLUBB_core/endian.F90
ensemble_manager.o: $(SRCROOT)shared/coupler/ensemble_manager.F90 fms.o mpp.o fms_io.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)shared/coupler/ensemble_manager.F90
entrain.o: $(SRCROOT)atmos_param/entrain/entrain.F90 constants.o mpp.o fms.o diag_manager.o time_manager.o sat_vapor_pres.o monin_obukhov.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)atmos_param/entrain/entrain.F90
epsbalgebraic.o: $(SRCROOT)mom5/ocean_param/gotm-4.0/turbulence/epsbalgebraic.F90 /home/wenchang/tigress/CM2.5_tiger2_wy_killtc/exp/../src/mom5/ocean_param/gotm-4.0/include/cppdefs.h /home/wenchang/tigress/CM2.5_tiger2_wy_killtc/exp/../src/mom5/ocean_param/gotm-4.0/include/version.h turbulence.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -I$(SRCROOT)mom5/ocean_param/gotm-4.0/include	$(SRCROOT)mom5/ocean_param/gotm-4.0/turbulence/epsbalgebraic.F90
eqstate.o: $(SRCROOT)mom5/ocean_param/gotm-4.0/util/eqstate.F90 /home/wenchang/tigress/CM2.5_tiger2_wy_killtc/exp/../src/mom5/ocean_param/gotm-4.0/include/cppdefs.h /home/wenchang/tigress/CM2.5_tiger2_wy_killtc/exp/../src/mom5/ocean_param/gotm-4.0/include/version.h
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -I$(SRCROOT)mom5/ocean_param/gotm-4.0/include	$(SRCROOT)mom5/ocean_param/gotm-4.0/util/eqstate.F90
error_code.o: $(SRCROOT)atmos_param/clubb/CLUBB_core/error_code.F90 constants_clubb.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)atmos_param/clubb/CLUBB_core/error_code.F90
esfsw_driver.o: $(SRCROOT)atmos_param/sea_esf_rad/esfsw_driver.F90 mpp.o fms.o constants.o esfsw_parameters.o rad_utilities.o tracer_manager.o field_manager.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)atmos_param/sea_esf_rad/esfsw_driver.F90
esfsw_parameters.o: $(SRCROOT)atmos_param/sea_esf_rad/esfsw_parameters.F90 mpp.o fms.o rad_utilities.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)atmos_param/sea_esf_rad/esfsw_parameters.F90
external_ic.o: $(SRCROOT)atmos_cubed_sphere/tools/external_ic.F90 external_sst.o fms.o fms_io.o mpp.o mpp_parameter.o mpp_domains.o tracer_manager.o field_manager.o constants.o fv_arrays.o fv_diagnostics.o fv_grid_utils.o fv_io.o fv_mapz.o fv_mp_mod.o fv_surf_map.o fv_timing.o init_hydro.o sim_nc_mod.o fv_nudge.o boundary.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)atmos_cubed_sphere/tools/external_ic.F90
external_sst.o: $(SRCROOT)atmos_cubed_sphere/tools/external_sst.F90 amip_interp.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)atmos_cubed_sphere/tools/external_sst.F90
extrapolation.o: $(SRCROOT)atmos_param/clubb/CLUBB_core/extrapolation.F90 clubb_precision.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)atmos_param/clubb/CLUBB_core/extrapolation.F90
fft.o: $(SRCROOT)shared/fft/fft.F90 platform.o fms.o fft99.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)shared/fft/fft.F90
fft99.o: $(SRCROOT)shared/fft/fft99.F90 constants.o mpp.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)shared/fft/fft99.F90
field_manager.o: $(SRCROOT)shared/field_manager/field_manager.F90 /home/wenchang/tigress/CM2.5_tiger2_wy_killtc/exp/../src/shared/field_manager/parse.inc mpp.o mpp_io.o fms.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -I$(SRCROOT)shared/field_manager	$(SRCROOT)shared/field_manager/field_manager.F90
file_functions.o: $(SRCROOT)atmos_param/clubb/CLUBB_core/file_functions.F90 clubb_precision.o constants_clubb.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)atmos_param/clubb/CLUBB_core/file_functions.F90
fill_holes.o: $(SRCROOT)atmos_param/clubb/CLUBB_core/fill_holes.F90 grid_class.o clubb_precision.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)atmos_param/clubb/CLUBB_core/fill_holes.F90
fk_craig.o: $(SRCROOT)mom5/ocean_param/gotm-4.0/turbulence/fk_craig.F90 /home/wenchang/tigress/CM2.5_tiger2_wy_killtc/exp/../src/mom5/ocean_param/gotm-4.0/include/cppdefs.h /home/wenchang/tigress/CM2.5_tiger2_wy_killtc/exp/../src/mom5/ocean_param/gotm-4.0/include/version.h
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -I$(SRCROOT)mom5/ocean_param/gotm-4.0/include	$(SRCROOT)mom5/ocean_param/gotm-4.0/turbulence/fk_craig.F90
flux_exchange.o: $(SRCROOT)coupler/flux_exchange.F90 mpp.o mpp_domains.o mpp_io.o atmos_model.o ocean_model.o ice_model.o land_model.o surface_flux.o monin_obukhov.o xgrid.o diag_integral.o diag_manager.o time_manager.o sat_vapor_pres.o constants.o fms.o data_override.o coupler_types.o atmos_ocean_fluxes.o atmos_tracer_driver.o field_manager.o tracer_manager.o stock_constants.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)coupler/flux_exchange.F90
fm_util.o: $(SRCROOT)shared/field_manager/fm_util.F90 field_manager.o fms.o mpp.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)shared/field_manager/fm_util.F90
fms.o: $(SRCROOT)shared/fms/fms.F90 mpp.o mpp_domains.o mpp_io.o fms_io.o memutils.o constants.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)shared/fms/fms.F90
fms_donner.o: $(SRCROOT)atmos_param/donner_deep/fms_donner.F90 /home/wenchang/tigress/CM2.5_tiger2_wy_killtc/exp/../src/atmos_param/donner_deep/donner_nml.h time_manager.o diag_manager.o field_manager.o tracer_manager.o atmos_tracer_utilities.o sat_vapor_pres.o fms.o fms_io.o mpp.o mpp_io.o constants.o column_diagnostics.o donner_types.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -I$(SRCROOT)atmos_param/donner_deep	$(SRCROOT)atmos_param/donner_deep/fms_donner.F90
fms_io.o: $(SRCROOT)shared/fms/fms_io.F90 /home/wenchang/tigress/CM2.5_tiger2_wy_killtc/exp/../src/shared/include/fms_platform.h /home/wenchang/tigress/CM2.5_tiger2_wy_killtc/exp/../src/shared/fms/read_data_2d.inc /home/wenchang/tigress/CM2.5_tiger2_wy_killtc/exp/../src/shared/fms/read_data_3d.inc /home/wenchang/tigress/CM2.5_tiger2_wy_killtc/exp/../src/shared/fms/read_data_4d.inc /home/wenchang/tigress/CM2.5_tiger2_wy_killtc/exp/../src/shared/fms/write_data.inc mpp_io.o mpp_domains.o mpp.o platform.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -I$(SRCROOT)shared/include -I$(SRCROOT)shared/fms	$(SRCROOT)shared/fms/fms_io.F90
format_input.o: $(SRCROOT)atmos_param/cosp/quickbeam/format_input.F90 /home/wenchang/tigress/CM2.5_tiger2_wy_killtc/exp/../src/atmos_param/cosp/cosp_defs.H array_lib.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -I$(SRCROOT)atmos_param/cosp	$(SRCROOT)atmos_param/cosp/quickbeam/format_input.F90
fs_profile.o: $(SRCROOT)atmos_param/fsrad/fs_profile.F90 fms.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)atmos_param/fsrad/fs_profile.F90
fsrad.o: $(SRCROOT)atmos_param/fsrad/fsrad.F90 mcm_lw.o mcm_sw_driver.o shortwave.o longwave.o rdparm.o rad_diag.o co2_data.o fms.o constants.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)atmos_param/fsrad/fsrad.F90
fv_ada_nudge.o: $(SRCROOT)atmos_cubed_sphere/tools/fv_ada_nudge.F90 external_sst.o diag_manager.o constants.o fms.o mpp.o mpp_domains.o time_manager.o grid.o fv_grid_utils.o fv_diagnostics.o tp_core.o fv_mapz.o fv_mp_mod.o fv_timing.o sim_nc_mod.o fv_arrays.o fms_io.o axis_utils.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)atmos_cubed_sphere/tools/fv_ada_nudge.F90
fv_arrays.o: $(SRCROOT)atmos_cubed_sphere/model/fv_arrays.F90 /home/wenchang/tigress/CM2.5_tiger2_wy_killtc/exp/../src/shared/include/fms_platform.h mpp_domains.o fms_io.o time_manager.o horiz_interp_type.o mpp.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -I$(SRCROOT)shared/include	$(SRCROOT)atmos_cubed_sphere/model/fv_arrays.F90
fv_climate_nudge.o: $(SRCROOT)atmos_cubed_sphere/tools/fv_climate_nudge.F90 fms.o mpp.o diag_manager.o read_climate_nudge_data.o time_manager.o time_interp.o get_cal_time.o constants.o fv_mapz.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)atmos_cubed_sphere/tools/fv_climate_nudge.F90
fv_control.o: $(SRCROOT)atmos_cubed_sphere/model/fv_control.F90 constants.o field_manager.o fms.o mpp.o mpp_domains.o tracer_manager.o fv_io.o fv_restart.o fv_arrays.o fv_grid_utils.o fv_eta.o fv_grid_tools.o fv_mp_mod.o test_cases.o fv_timing.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)atmos_cubed_sphere/model/fv_control.F90
fv_diagnostics.o: $(SRCROOT)atmos_cubed_sphere/tools/fv_diagnostics.F90 constants.o fms_io.o time_manager.o mpp_domains.o diag_manager.o fv_arrays.o fv_mapz.o fv_mp_mod.o fv_eta.o fv_grid_utils.o a2b_edge.o fv_surf_map.o fv_sg.o tracer_manager.o field_manager.o mpp.o sat_vapor_pres.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)atmos_cubed_sphere/tools/fv_diagnostics.F90
fv_dynamics.o: $(SRCROOT)atmos_cubed_sphere/model/fv_dynamics.F90 constants.o dyn_core.o fv_mapz.o fv_tracer2d.o fv_grid_utils.o fv_mp_mod.o fv_timing.o diag_manager.o fv_diagnostics.o mpp_domains.o field_manager.o tracer_manager.o fv_sg.o fv_nesting.o fv_arrays.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)atmos_cubed_sphere/model/fv_dynamics.F90
fv_eta.o: $(SRCROOT)atmos_cubed_sphere/tools/fv_eta.F90 constants.o fv_mp_mod.o mpp.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)atmos_cubed_sphere/tools/fv_eta.F90
fv_fill.o: $(SRCROOT)atmos_cubed_sphere/model/fv_fill.F90
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)atmos_cubed_sphere/model/fv_fill.F90
fv_grid_tools.o: $(SRCROOT)atmos_cubed_sphere/tools/fv_grid_tools.F90 constants.o fv_arrays.o fv_grid_utils.o fv_timing.o fv_mp_mod.o sorted_index.o mpp.o mpp_domains.o mpp_io.o mpp_parameter.o fms.o fms_io.o mosaic.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)atmos_cubed_sphere/tools/fv_grid_tools.F90
fv_grid_utils.o: $(SRCROOT)atmos_cubed_sphere/model/fv_grid_utils.F90 /home/wenchang/tigress/CM2.5_tiger2_wy_killtc/exp/../src/shared/include/fms_platform.h constants.o mpp.o external_sst.o mpp_domains.o mpp_parameter.o fv_arrays.o fv_eta.o fv_mp_mod.o fv_timing.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -I$(SRCROOT)shared/include	$(SRCROOT)atmos_cubed_sphere/model/fv_grid_utils.F90
fv_io.o: $(SRCROOT)atmos_cubed_sphere/tools/fv_io.F90 fms.o fms_io.o mpp.o mpp_domains.o tracer_manager.o field_manager.o external_sst.o fv_arrays.o fv_eta.o fv_mp_mod.o fv_mapz.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)atmos_cubed_sphere/tools/fv_io.F90
fv_mapz.o: $(SRCROOT)atmos_cubed_sphere/model/fv_mapz.F90 constants.o tracer_manager.o field_manager.o fv_grid_utils.o fv_fill.o fv_eta.o mpp_domains.o mpp.o fv_arrays.o lin_cloud_microphys.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)atmos_cubed_sphere/model/fv_mapz.F90
fv_mp_mod.o: $(SRCROOT)atmos_cubed_sphere/tools/fv_mp_mod.F90 fms.o mpp.o mpp_domains.o mpp_parameter.o fv_arrays.o fms_io.o ensemble_manager.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)atmos_cubed_sphere/tools/fv_mp_mod.F90
fv_nesting.o: $(SRCROOT)atmos_cubed_sphere/model/fv_nesting.F90 mpp_domains.o field_manager.o tracer_manager.o fv_sg.o fv_restart.o mpp.o boundary.o fv_mp_mod.o fv_arrays.o fv_grid_utils.o init_hydro.o constants.o fv_mapz.o fv_timing.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)atmos_cubed_sphere/model/fv_nesting.F90
fv_nudge.o: $(SRCROOT)atmos_cubed_sphere/tools/fv_nudge.F90 external_sst.o diag_manager.o constants.o fms.o mpp.o mpp_domains.o time_manager.o fv_grid_utils.o fv_diagnostics.o tp_core.o fv_mapz.o fv_mp_mod.o fv_timing.o sim_nc_mod.o fv_arrays.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)atmos_cubed_sphere/tools/fv_nudge.F90
fv_physics.o: $(SRCROOT)atmos_cubed_sphere/driver/coupled/fv_physics.F90 atmos_co2.o constants.o time_manager.o fms.o physics_driver.o field_manager.o tracer_manager.o mpp_domains.o mpp.o fv_grid_utils.o fv_arrays.o fv_eta.o fv_update_phys.o fv_sg.o fv_timing.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)atmos_cubed_sphere/driver/coupled/fv_physics.F90
fv_restart.o: $(SRCROOT)atmos_cubed_sphere/tools/fv_restart.F90 constants.o fv_arrays.o fv_io.o fv_grid_utils.o fv_diagnostics.o init_hydro.o mpp_domains.o mpp.o test_cases.o fv_mp_mod.o fv_surf_map.o tracer_manager.o field_manager.o external_ic.o fv_eta.o boundary.o fv_timing.o fms.o external_sst.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)atmos_cubed_sphere/tools/fv_restart.F90
fv_sg.o: $(SRCROOT)atmos_cubed_sphere/model/fv_sg.F90 constants.o tracer_manager.o field_manager.o lin_cloud_microphys.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)atmos_cubed_sphere/model/fv_sg.F90
fv_surf_map.o: $(SRCROOT)atmos_cubed_sphere/tools/fv_surf_map.F90 fms.o mpp.o mpp_domains.o constants.o fms_io.o fv_grid_utils.o fv_mp_mod.o fv_timing.o fv_arrays.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)atmos_cubed_sphere/tools/fv_surf_map.F90
fv_timing.o: $(SRCROOT)atmos_cubed_sphere/tools/fv_timing.F90 mpp.o fv_mp_mod.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)atmos_cubed_sphere/tools/fv_timing.F90
fv_tracer2d.o: $(SRCROOT)atmos_cubed_sphere/model/fv_tracer2d.F90 tp_core.o fv_mp_mod.o mpp_domains.o fv_timing.o boundary.o fv_arrays.o mpp.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)atmos_cubed_sphere/model/fv_tracer2d.F90
fv_update_phys.o: $(SRCROOT)atmos_cubed_sphere/model/fv_update_phys.F90 constants.o field_manager.o mpp_domains.o mpp_parameter.o mpp.o time_manager.o tracer_manager.o fv_arrays.o boundary.o fv_eta.o fv_timing.o fv_diagnostics.o atmos_nudge.o fv_climate_nudge.o fv_ada_nudge.o fv_nudge.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)atmos_cubed_sphere/model/fv_update_phys.F90
gamma_mg.o: $(SRCROOT)atmos_param/strat_cloud/gamma_mg.F90 fms.o mg_const.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)atmos_param/strat_cloud/gamma_mg.F90
gas_tf.o: $(SRCROOT)atmos_param/sea_esf_rad/gas_tf.F90 mpp.o fms.o constants.o rad_utilities.o longwave_params.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)atmos_param/sea_esf_rad/gas_tf.F90
gases.o: $(SRCROOT)atmos_param/cosp/quickbeam/gases.F90 /home/wenchang/tigress/CM2.5_tiger2_wy_killtc/exp/../src/atmos_param/cosp/cosp_defs.H
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -I$(SRCROOT)atmos_param/cosp	$(SRCROOT)atmos_param/cosp/quickbeam/gases.F90
gaussian_topog.o: $(SRCROOT)shared/topography/gaussian_topog.F90 fms.o constants.o mpp.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)shared/topography/gaussian_topog.F90
generic_BLING.o: $(SRCROOT)ocean_shared/generic_tracers/generic_BLING.F90 coupler_types.o field_manager.o mpp.o fms.o time_manager.o fm_util.o diag_manager.o constants.o generic_tracer_utils.o FMS_ocmip2_co2calc.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)ocean_shared/generic_tracers/generic_BLING.F90
generic_CFC.o: $(SRCROOT)ocean_shared/generic_tracers/generic_CFC.F90 coupler_types.o field_manager.o mpp.o time_manager.o fm_util.o generic_tracer_utils.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)ocean_shared/generic_tracers/generic_CFC.F90
generic_COBALT.o: $(SRCROOT)ocean_shared/generic_tracers/generic_COBALT.F90 coupler_types.o field_manager.o mpp.o time_manager.o fm_util.o diag_manager.o constants.o fms.o generic_tracer_utils.o FMS_ocmip2_co2calc.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)ocean_shared/generic_tracers/generic_COBALT.F90
generic_ERGOM.o: $(SRCROOT)ocean_shared/generic_tracers/generic_ERGOM.F90 coupler_types.o field_manager.o fms.o mpp.o time_manager.o fm_util.o diag_manager.o constants.o generic_tracer_utils.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)ocean_shared/generic_tracers/generic_ERGOM.F90
generic_TOPAZ.o: $(SRCROOT)ocean_shared/generic_tracers/generic_TOPAZ.F90 coupler_types.o field_manager.o mpp.o time_manager.o fm_util.o diag_manager.o constants.o generic_tracer_utils.o FMS_ocmip2_co2calc.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)ocean_shared/generic_tracers/generic_TOPAZ.F90
generic_miniBLING.o: $(SRCROOT)ocean_shared/generic_tracers/generic_miniBLING.F90 /home/wenchang/tigress/CM2.5_tiger2_wy_killtc/exp/../src/shared/include/fms_platform.h coupler_types.o field_manager.o mpp.o time_manager.o diag_manager.o constants.o data_override.o generic_tracer_utils.o FMS_ocmip2_co2calc.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -I$(SRCROOT)shared/include	$(SRCROOT)ocean_shared/generic_tracers/generic_miniBLING.F90
generic_tracer.o: $(SRCROOT)ocean_shared/generic_tracers/generic_tracer.F90 fms.o field_manager.o mpp.o time_manager.o coupler_types.o generic_tracer_utils.o generic_CFC.o generic_ERGOM.o generic_TOPAZ.o generic_BLING.o generic_miniBLING.o generic_COBALT.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)ocean_shared/generic_tracers/generic_tracer.F90
generic_tracer_utils.o: $(SRCROOT)ocean_shared/generic_tracers/generic_tracer_utils.F90 /home/wenchang/tigress/CM2.5_tiger2_wy_killtc/exp/../src/shared/include/fms_platform.h coupler_types.o FMS_coupler_util.o atmos_ocean_fluxes.o mpp.o time_manager.o diag_manager.o field_manager.o fms.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -I$(SRCROOT)shared/include	$(SRCROOT)ocean_shared/generic_tracers/generic_tracer_utils.F90
genericeq.o: $(SRCROOT)mom5/ocean_param/gotm-4.0/turbulence/genericeq.F90 /home/wenchang/tigress/CM2.5_tiger2_wy_killtc/exp/../src/mom5/ocean_param/gotm-4.0/include/cppdefs.h /home/wenchang/tigress/CM2.5_tiger2_wy_killtc/exp/../src/mom5/ocean_param/gotm-4.0/include/version.h turbulence.o util.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -I$(SRCROOT)mom5/ocean_param/gotm-4.0/include	$(SRCROOT)mom5/ocean_param/gotm-4.0/turbulence/genericeq.F90
get_cal_time.o: $(SRCROOT)shared/time_manager/get_cal_time.F90 fms.o time_manager.o mpp.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)shared/time_manager/get_cal_time.F90
glac_tile.o: $(SRCROOT)land_lad2/glacier/glac_tile.F90 /home/wenchang/tigress/CM2.5_tiger2_wy_killtc/exp/../src/shared/include/fms_platform.h mpp.o fms.o constants.o land_constants.o land_io.o land_tile_diag_sel.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -I$(SRCROOT)shared/include	$(SRCROOT)land_lad2/glacier/glac_tile.F90
glacier.o: $(SRCROOT)land_lad2/glacier/glacier.F90 mpp.o fms.o fms_io.o time_manager.o diag_manager.o constants.o glac_tile.o land_constants.o land_tile.o land_tile_diag.o land_data.o land_io.o land_tile_io.o nf_utils.o land_debug.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)land_lad2/glacier/glacier.F90
gmres_cache.o: $(SRCROOT)atmos_param/clubb/CLUBB_core/gmres_cache.F90 clubb_precision.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)atmos_param/clubb/CLUBB_core/gmres_cache.F90
gmres_wrap.o: $(SRCROOT)atmos_param/clubb/CLUBB_core/gmres_wrap.F90 gmres_cache.o clubb_precision.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)atmos_param/clubb/CLUBB_core/gmres_wrap.F90
gotm_lib_version.o: $(SRCROOT)mom5/ocean_param/gotm-4.0/turbulence/gotm_lib_version.F90 /home/wenchang/tigress/CM2.5_tiger2_wy_killtc/exp/../src/mom5/ocean_param/gotm-4.0/include/cppdefs.h /home/wenchang/tigress/CM2.5_tiger2_wy_killtc/exp/../src/mom5/ocean_param/gotm-4.0/include/version.h
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -I$(SRCROOT)mom5/ocean_param/gotm-4.0/include	$(SRCROOT)mom5/ocean_param/gotm-4.0/turbulence/gotm_lib_version.F90
gradient.o: $(SRCROOT)shared/mosaic/gradient.F90 mpp.o constants.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)shared/mosaic/gradient.F90
gradient_c2l.o: $(SRCROOT)shared/mosaic/gradient_c2l.c /home/wenchang/tigress/CM2.5_tiger2_wy_killtc/exp/../src/shared/mosaic/constant.h /home/wenchang/tigress/CM2.5_tiger2_wy_killtc/exp/../src/shared/mosaic/mosaic_util.h /home/wenchang/tigress/CM2.5_tiger2_wy_killtc/exp/../src/shared/mosaic/gradient_c2l.h
	$(CC) $(CPPDEFS) $(CPPFLAGS) $(CFLAGS) $(OTHERFLAGS) $(OTHER_CFLAGS) -c -I$(SRCROOT)shared/mosaic	$(SRCROOT)shared/mosaic/gradient_c2l.c
grey_radiation.o: $(SRCROOT)atmos_param/grey_radiation/grey_radiation.F90 mpp.o fms.o constants.o astronomy.o diag_manager.o time_manager.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)atmos_param/grey_radiation/grey_radiation.F90
grid.o: $(SRCROOT)shared/mosaic/grid.F90 mpp.o constants.o fms.o fms_io.o mosaic.o mpp_domains.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)shared/mosaic/grid.F90
grid_class.o: $(SRCROOT)atmos_param/clubb/CLUBB_core/grid_class.F90 clubb_precision.o constants_clubb.o error_code.o file_functions.o interpolation.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)atmos_param/clubb/CLUBB_core/grid_class.F90
gridinterpol.o: $(SRCROOT)mom5/ocean_param/gotm-4.0/util/gridinterpol.F90 /home/wenchang/tigress/CM2.5_tiger2_wy_killtc/exp/../src/mom5/ocean_param/gotm-4.0/include/cppdefs.h /home/wenchang/tigress/CM2.5_tiger2_wy_killtc/exp/../src/mom5/ocean_param/gotm-4.0/include/version.h
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -I$(SRCROOT)mom5/ocean_param/gotm-4.0/include	$(SRCROOT)mom5/ocean_param/gotm-4.0/util/gridinterpol.F90
hconst.o: $(SRCROOT)atmos_param/fsrad/hconst.F90
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)atmos_param/fsrad/hconst.F90
horiz_interp.o: $(SRCROOT)shared/horiz_interp/horiz_interp.F90 fms.o mpp.o constants.o horiz_interp_type.o horiz_interp_conserve.o horiz_interp_bilinear.o horiz_interp_bicubic.o horiz_interp_spherical.o mpp_io.o mpp_domains.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)shared/horiz_interp/horiz_interp.F90
horiz_interp_bicubic.o: $(SRCROOT)shared/horiz_interp/horiz_interp_bicubic.F90 mpp.o fms.o horiz_interp_type.o constants.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)shared/horiz_interp/horiz_interp_bicubic.F90
horiz_interp_bilinear.o: $(SRCROOT)shared/horiz_interp/horiz_interp_bilinear.F90 mpp.o fms.o constants.o horiz_interp_type.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)shared/horiz_interp/horiz_interp_bilinear.F90
horiz_interp_conserve.o: $(SRCROOT)shared/horiz_interp/horiz_interp_conserve.F90 mpp.o fms.o fms_io.o constants.o horiz_interp_type.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)shared/horiz_interp/horiz_interp_conserve.F90
horiz_interp_spherical.o: $(SRCROOT)shared/horiz_interp/horiz_interp_spherical.F90 mpp.o fms.o constants.o horiz_interp_type.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)shared/horiz_interp/horiz_interp_spherical.F90
horiz_interp_type.o: $(SRCROOT)shared/horiz_interp/horiz_interp_type.F90 mpp.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)shared/horiz_interp/horiz_interp_type.F90
hydrostatic_module.o: $(SRCROOT)atmos_param/clubb/CLUBB_core/hydrostatic_module.F90 constants_clubb.o grid_class.o clubb_precision.o interpolation.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)atmos_param/clubb/CLUBB_core/hydrostatic_module.F90
hyper_diffusion_4th_ord.o: $(SRCROOT)atmos_param/clubb/CLUBB_core/hyper_diffusion_4th_ord.F90 clubb_precision.o grid_class.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)atmos_param/clubb/CLUBB_core/hyper_diffusion_4th_ord.F90
icarus.o: $(SRCROOT)atmos_param/cosp/icarus-scops-3.7/icarus.F90 /home/wenchang/tigress/CM2.5_tiger2_wy_killtc/exp/../src/atmos_param/cosp/cosp_defs.H mpp.o fms.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -I$(SRCROOT)atmos_param/cosp	$(SRCROOT)atmos_param/cosp/icarus-scops-3.7/icarus.F90
ice_albedo.o: $(SRCROOT)ice_param/ice_albedo.F90 mpp.o fms.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)ice_param/ice_albedo.F90
ice_bergs.o: $(SRCROOT)ice_sis/ice_bergs.F90 constants.o fms.o mosaic.o mpp.o fms_io.o mpp_domains.o mpp_parameter.o time_manager.o diag_manager.o random_numbers.o stock_constants.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)ice_sis/ice_bergs.F90
ice_dyn.o: $(SRCROOT)ice_sis/ice_dyn.F90 mpp_domains.o constants.o ice_grid.o ice_thm.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)ice_sis/ice_dyn.F90
ice_grid.o: $(SRCROOT)ice_sis/ice_grid.F90 constants.o mpp.o mpp_domains.o fms.o mosaic.o grid.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)ice_sis/ice_grid.F90
ice_model.o: $(SRCROOT)ice_sis/ice_model.F90 mpp.o mpp_domains.o fms.o diag_manager.o time_manager.o astronomy.o coupler_types.o constants.o ocean_albedo.o ocean_rough.o ice_type.o ice_grid.o ice_spec.o ice_thm.o ice_dyn.o ice_bergs.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)ice_sis/ice_model.F90
ice_nucl.o: $(SRCROOT)atmos_shared/tracer_driver/aer_ccn_act/ice_nucl.F90 mpp.o fms.o aer_ccn_act_k.o aerosol_params.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)atmos_shared/tracer_driver/aer_ccn_act/ice_nucl.F90
ice_spec.o: $(SRCROOT)ice_sis/ice_spec.F90 fms.o mpp.o time_manager.o data_override.o constants.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)ice_sis/ice_spec.F90
ice_thm.o: $(SRCROOT)ice_sis/ice_thm.F90 constants.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)ice_sis/ice_thm.F90
ice_type.o: $(SRCROOT)ice_sis/ice_type.F90 mpp.o mpp_domains.o fms.o fms_io.o diag_manager.o time_manager.o coupler_types.o constants.o ice_grid.o ice_thm.o ice_dyn.o ice_bergs.o astronomy.o stock_constants.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)ice_sis/ice_type.F90
init_hydro.o: $(SRCROOT)atmos_cubed_sphere/tools/init_hydro.F90 constants.o fv_grid_utils.o fv_mp_mod.o field_manager.o tracer_manager.o mpp_domains.o mpp.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)atmos_cubed_sphere/tools/init_hydro.F90
input_names.o: $(SRCROOT)atmos_param/clubb/CLUBB_core/input_names.F90
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)atmos_param/clubb/CLUBB_core/input_names.F90
input_reader.o: $(SRCROOT)atmos_param/clubb/CLUBB_core/input_reader.F90 clubb_precision.o constants_clubb.o input_names.o interpolation.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)atmos_param/clubb/CLUBB_core/input_reader.F90
internal_wave.o: $(SRCROOT)mom5/ocean_param/gotm-4.0/turbulence/internal_wave.F90 /home/wenchang/tigress/CM2.5_tiger2_wy_killtc/exp/../src/mom5/ocean_param/gotm-4.0/include/cppdefs.h /home/wenchang/tigress/CM2.5_tiger2_wy_killtc/exp/../src/mom5/ocean_param/gotm-4.0/include/version.h turbulence.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -I$(SRCROOT)mom5/ocean_param/gotm-4.0/include	$(SRCROOT)mom5/ocean_param/gotm-4.0/turbulence/internal_wave.F90
interp.o: $(SRCROOT)shared/mosaic/interp.c /home/wenchang/tigress/CM2.5_tiger2_wy_killtc/exp/../src/shared/mosaic/mosaic_util.h /home/wenchang/tigress/CM2.5_tiger2_wy_killtc/exp/../src/shared/mosaic/interp.h /home/wenchang/tigress/CM2.5_tiger2_wy_killtc/exp/../src/shared/mosaic/create_xgrid.h
	$(CC) $(CPPDEFS) $(CPPFLAGS) $(CFLAGS) $(OTHERFLAGS) $(OTHER_CFLAGS) -c -I$(SRCROOT)shared/mosaic	$(SRCROOT)shared/mosaic/interp.c
interpolation.o: $(SRCROOT)atmos_param/clubb/CLUBB_core/interpolation.F90 clubb_precision.o constants_clubb.o model_flags.o error_code.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)atmos_param/clubb/CLUBB_core/interpolation.F90
interpolator.o: $(SRCROOT)atmos_shared/interpolator/interpolator.F90 /home/wenchang/tigress/CM2.5_tiger2_wy_killtc/exp/../src/shared/include/fms_platform.h mpp.o mpp_io.o mpp_domains.o diag_manager.o fms.o horiz_interp.o time_manager.o time_interp.o constants.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -I$(SRCROOT)shared/include	$(SRCROOT)atmos_shared/interpolator/interpolator.F90
isccp_cloud_types.o: $(SRCROOT)atmos_param/cosp/icarus-scops-3.7/isccp_cloud_types.F /home/wenchang/tigress/CM2.5_tiger2_wy_killtc/exp/../src/atmos_param/cosp/cosp_defs.H
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -I$(SRCROOT)atmos_param/cosp	$(SRCROOT)atmos_param/cosp/icarus-scops-3.7/isccp_cloud_types.F
isccp_clouds.o: $(SRCROOT)atmos_param/sea_esf_rad/isccp_clouds.F90 mpp.o fms.o time_manager.o diag_manager.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)atmos_param/sea_esf_rad/isccp_clouds.F90
ispralength.o: $(SRCROOT)mom5/ocean_param/gotm-4.0/turbulence/ispralength.F90 /home/wenchang/tigress/CM2.5_tiger2_wy_killtc/exp/../src/mom5/ocean_param/gotm-4.0/include/cppdefs.h /home/wenchang/tigress/CM2.5_tiger2_wy_killtc/exp/../src/mom5/ocean_param/gotm-4.0/include/version.h turbulence.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -I$(SRCROOT)mom5/ocean_param/gotm-4.0/include	$(SRCROOT)mom5/ocean_param/gotm-4.0/turbulence/ispralength.F90
kbalgebraic.o: $(SRCROOT)mom5/ocean_param/gotm-4.0/turbulence/kbalgebraic.F90 /home/wenchang/tigress/CM2.5_tiger2_wy_killtc/exp/../src/mom5/ocean_param/gotm-4.0/include/cppdefs.h /home/wenchang/tigress/CM2.5_tiger2_wy_killtc/exp/../src/mom5/ocean_param/gotm-4.0/include/version.h turbulence.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -I$(SRCROOT)mom5/ocean_param/gotm-4.0/include	$(SRCROOT)mom5/ocean_param/gotm-4.0/turbulence/kbalgebraic.F90
kbeq.o: $(SRCROOT)mom5/ocean_param/gotm-4.0/turbulence/kbeq.F90 /home/wenchang/tigress/CM2.5_tiger2_wy_killtc/exp/../src/mom5/ocean_param/gotm-4.0/include/cppdefs.h /home/wenchang/tigress/CM2.5_tiger2_wy_killtc/exp/../src/mom5/ocean_param/gotm-4.0/include/version.h turbulence.o util.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -I$(SRCROOT)mom5/ocean_param/gotm-4.0/include	$(SRCROOT)mom5/ocean_param/gotm-4.0/turbulence/kbeq.F90
kpp.o: $(SRCROOT)mom5/ocean_param/gotm-4.0/turbulence/kpp.F90 /home/wenchang/tigress/CM2.5_tiger2_wy_killtc/exp/../src/mom5/ocean_param/gotm-4.0/include/cppdefs.h /home/wenchang/tigress/CM2.5_tiger2_wy_killtc/exp/../src/mom5/ocean_param/gotm-4.0/include/version.h turbulence.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -I$(SRCROOT)mom5/ocean_param/gotm-4.0/include	$(SRCROOT)mom5/ocean_param/gotm-4.0/turbulence/kpp.F90
lagrange.o: $(SRCROOT)mom5/ocean_param/gotm-4.0/util/lagrange.F90 /home/wenchang/tigress/CM2.5_tiger2_wy_killtc/exp/../src/mom5/ocean_param/gotm-4.0/include/cppdefs.h /home/wenchang/tigress/CM2.5_tiger2_wy_killtc/exp/../src/mom5/ocean_param/gotm-4.0/include/version.h
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -I$(SRCROOT)mom5/ocean_param/gotm-4.0/include	$(SRCROOT)mom5/ocean_param/gotm-4.0/util/lagrange.F90
lake.o: $(SRCROOT)land_lad2/lake/lake.F90 mpp.o fms.o fms_io.o time_manager.o diag_manager.o constants.o land_constants.o land_io.o lake_tile.o land_tile.o land_tile_diag.o land_data.o land_tile_io.o nf_utils.o land_debug.o land_utils.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)land_lad2/lake/lake.F90
lake_tile.o: $(SRCROOT)land_lad2/lake/lake_tile.F90 /home/wenchang/tigress/CM2.5_tiger2_wy_killtc/exp/../src/shared/include/fms_platform.h mpp_domains.o mpp.o fms.o constants.o land_constants.o land_io.o land_tile_diag_sel.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -I$(SRCROOT)shared/include	$(SRCROOT)land_lad2/lake/lake_tile.F90
land_constants.o: $(SRCROOT)land_lad2/land_constants.F90 constants.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)land_lad2/land_constants.F90
land_data.o: $(SRCROOT)land_lad2/land_data.F90 mpp.o constants.o mpp_domains.o fms.o time_manager.o tracer_manager.o field_manager.o grid.o land_tile.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)land_lad2/land_data.F90
land_debug.o: $(SRCROOT)land_lad2/shared/land_debug.F90 mpp.o fms.o time_manager.o grid.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)land_lad2/shared/land_debug.F90
land_io.o: $(SRCROOT)land_lad2/shared/land_io.F90 mpp_io.o axis_utils.o constants.o fms.o mpp.o horiz_interp.o land_numerics.o nf_utils.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)land_lad2/shared/land_io.F90
land_model.o: $(SRCROOT)land_lad2/land_model.F90 /home/wenchang/tigress/CM2.5_tiger2_wy_killtc/exp/../src/land_lad2/shared/debug.inc time_manager.o mpp_domains.o mpp.o fms.o fms_io.o field_manager.o data_override.o diag_manager.o constants.o astronomy.o sphum.o tracer_manager.o land_constants.o glacier.o lake.o soil.o snow.o vegetation.o cana_tile.o canopy_air.o river.o topo_rough.o soil_tile.o vegn_tile.o lake_tile.o glac_tile.o snow_tile.o land_numerics.o land_io.o land_tile.o land_data.o nf_utils.o land_utils.o land_tile_io.o land_tile_diag.o land_debug.o vegn_static_override.o transitions.o stock_constants.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -I$(SRCROOT)land_lad2/shared	$(SRCROOT)land_lad2/land_model.F90
land_numerics.o: $(SRCROOT)land_lad2/shared/land_numerics.F90 fms.o mpp.o mpp_domains.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)land_lad2/shared/land_numerics.F90
land_tile.o: $(SRCROOT)land_lad2/land_tile.F90 fms.o land_constants.o glac_tile.o lake_tile.o soil_tile.o cana_tile.o vegn_tile.o snow_tile.o land_tile_diag_sel.o land_tile_diag_buff.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)land_lad2/land_tile.F90
land_tile_diag.o: $(SRCROOT)land_lad2/shared/land_tile_diag.F90 mpp.o time_manager.o diag_axis.o diag_manager.o diag_util.o fms.o land_tile_diag_sel.o land_tile.o land_data.o land_tile_diag_buff.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)land_lad2/shared/land_tile_diag.F90
land_tile_diag_buff.o: $(SRCROOT)land_lad2/shared/land_tile_diag_buff.F90
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)land_lad2/shared/land_tile_diag_buff.F90
land_tile_diag_sel.o: $(SRCROOT)land_lad2/shared/land_tile_diag_sel.F90 fms.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)land_lad2/shared/land_tile_diag_sel.F90
land_tile_io.o: $(SRCROOT)land_lad2/shared/land_tile_io.F90 mpp.o fms.o fms_io.o time_manager.o data_override.o nf_utils.o land_io.o land_tile.o land_data.o land_utils.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)land_lad2/shared/land_tile_io.F90
land_utils.o: $(SRCROOT)land_lad2/shared/land_utils.F90 land_tile.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)land_lad2/shared/land_utils.F90
lapack_wrap.o: $(SRCROOT)atmos_param/clubb/CLUBB_core/lapack_wrap.F90 constants_clubb.o error_code.o clubb_precision.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)atmos_param/clubb/CLUBB_core/lapack_wrap.F90
lengthscaleeq.o: $(SRCROOT)mom5/ocean_param/gotm-4.0/turbulence/lengthscaleeq.F90 /home/wenchang/tigress/CM2.5_tiger2_wy_killtc/exp/../src/mom5/ocean_param/gotm-4.0/include/cppdefs.h /home/wenchang/tigress/CM2.5_tiger2_wy_killtc/exp/../src/mom5/ocean_param/gotm-4.0/include/version.h turbulence.o util.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -I$(SRCROOT)mom5/ocean_param/gotm-4.0/include	$(SRCROOT)mom5/ocean_param/gotm-4.0/turbulence/lengthscaleeq.F90
lhsw_driver.o: $(SRCROOT)atmos_param/sea_esf_rad/lhsw_driver.F90 rad_utilities.o constants.o mpp.o fms.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)atmos_param/sea_esf_rad/lhsw_driver.F90
lidar_simulator.o: $(SRCROOT)atmos_param/cosp/actsim/lidar_simulator.F90 /home/wenchang/tigress/CM2.5_tiger2_wy_killtc/exp/../src/atmos_param/cosp/cosp_defs.H
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -I$(SRCROOT)atmos_param/cosp	$(SRCROOT)atmos_param/cosp/actsim/lidar_simulator.F90
lin_cloud_microphys.o: $(SRCROOT)atmos_param/lin_cloud_microphys/lin_cloud_microphys.F90 time_manager.o fms.o mpp.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)atmos_param/lin_cloud_microphys/lin_cloud_microphys.F90
llnl_stats.o: $(SRCROOT)atmos_param/cosp/llnl/llnl_stats.F90 /home/wenchang/tigress/CM2.5_tiger2_wy_killtc/exp/../src/atmos_param/cosp/cosp_defs.H cosp_constants.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -I$(SRCROOT)atmos_param/cosp	$(SRCROOT)atmos_param/cosp/llnl/llnl_stats.F90
lmd_ipsl_stats.o: $(SRCROOT)atmos_param/cosp/actsim/lmd_ipsl_stats.F90 /home/wenchang/tigress/CM2.5_tiger2_wy_killtc/exp/../src/atmos_param/cosp/cosp_defs.H llnl_stats.o cosp_constants.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -I$(SRCROOT)atmos_param/cosp	$(SRCROOT)atmos_param/cosp/actsim/lmd_ipsl_stats.F90
load_hydrometeor_classes.o: $(SRCROOT)atmos_param/cosp/quickbeam/load_hydrometeor_classes.F90 /home/wenchang/tigress/CM2.5_tiger2_wy_killtc/exp/../src/atmos_param/cosp/cosp_defs.H radar_simulator_types.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -I$(SRCROOT)atmos_param/cosp	$(SRCROOT)atmos_param/cosp/quickbeam/load_hydrometeor_classes.F90
load_mie_table.o: $(SRCROOT)atmos_param/cosp/quickbeam/load_mie_table.F90 /home/wenchang/tigress/CM2.5_tiger2_wy_killtc/exp/../src/atmos_param/cosp/cosp_defs.H radar_simulator_types.o mpp.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -I$(SRCROOT)atmos_param/cosp	$(SRCROOT)atmos_param/cosp/quickbeam/load_mie_table.F90
longwave.o: $(SRCROOT)atmos_param/fsrad/longwave.F90 rdparm.o hconst.o fms.o co2_data.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)atmos_param/fsrad/longwave.F90
longwave_clouds.o: $(SRCROOT)atmos_param/sea_esf_rad/longwave_clouds.F90 mpp.o fms.o constants.o rad_utilities.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)atmos_param/sea_esf_rad/longwave_clouds.F90
longwave_driver.o: $(SRCROOT)atmos_param/sea_esf_rad/longwave_driver.F90 mpp.o fms.o time_manager.o rad_utilities.o sealw99.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)atmos_param/sea_esf_rad/longwave_driver.F90
longwave_fluxes.o: $(SRCROOT)atmos_param/sea_esf_rad/longwave_fluxes.F90 mpp.o fms.o rad_utilities.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)atmos_param/sea_esf_rad/longwave_fluxes.F90
longwave_params.o: $(SRCROOT)atmos_param/sea_esf_rad/longwave_params.F90 mpp.o fms.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)atmos_param/sea_esf_rad/longwave_params.F90
longwave_tables.o: $(SRCROOT)atmos_param/sea_esf_rad/longwave_tables.F90 mpp.o fms.o rad_utilities.o longwave_params.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)atmos_param/sea_esf_rad/longwave_tables.F90
lscale_cond.o: $(SRCROOT)atmos_param/lscale_cond/lscale_cond.F90 mpp.o fms.o sat_vapor_pres.o constants.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)atmos_param/lscale_cond/lscale_cond.F90
lw_gases_stdtf.o: $(SRCROOT)atmos_param/sea_esf_rad/lw_gases_stdtf.F90 mpp.o fms.o fms_io.o rad_utilities.o gas_tf.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)atmos_param/sea_esf_rad/lw_gases_stdtf.F90
m_tracname.o: $(SRCROOT)atmos_shared/tracer_driver/tropchem/m_tracname.F90 moz.mods.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)atmos_shared/tracer_driver/tropchem/m_tracname.F90
mask.o: $(SRCROOT)ice_sis/mask.F90 rot.o fms.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)ice_sis/mask.F90
math_lib.o: $(SRCROOT)atmos_param/cosp/quickbeam/math_lib.F90 /home/wenchang/tigress/CM2.5_tiger2_wy_killtc/exp/../src/atmos_param/cosp/cosp_defs.H mrgrnk.o array_lib.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -I$(SRCROOT)atmos_param/cosp	$(SRCROOT)atmos_param/cosp/quickbeam/math_lib.F90
matrix_operations.o: $(SRCROOT)atmos_param/clubb/CLUBB_core/matrix_operations.F90 clubb_precision.o error_code.o constants_clubb.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)atmos_param/clubb/CLUBB_core/matrix_operations.F90
mcm_lw.o: $(SRCROOT)atmos_param/fsrad/mcm_lw.F90 constants.o fms.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)atmos_param/fsrad/mcm_lw.F90
mcm_sw_driver.o: $(SRCROOT)atmos_param/fsrad/mcm_sw_driver.F90 fms.o mcm_swnew.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)atmos_param/fsrad/mcm_sw_driver.F90
mcm_swnew.o: $(SRCROOT)atmos_param/fsrad/mcm_swnew.F90 mcm_swtbls.o fms.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)atmos_param/fsrad/mcm_swnew.F90
mcm_swtbls.o: $(SRCROOT)atmos_param/fsrad/mcm_swtbls.F90
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)atmos_param/fsrad/mcm_swtbls.F90
mean_adv.o: $(SRCROOT)atmos_param/clubb/CLUBB_core/mean_adv.F90 grid_class.o model_flags.o clubb_precision.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)atmos_param/clubb/CLUBB_core/mean_adv.F90
memuse.o: $(SRCROOT)shared/memutils/memuse.c
	$(CC) $(CPPDEFS) $(CPPFLAGS) $(CFLAGS) $(OTHERFLAGS) $(OTHER_CFLAGS) -c	$(SRCROOT)shared/memutils/memuse.c
memutils.o: $(SRCROOT)shared/memutils/memutils.F90 mpp.o mpp_io.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)shared/memutils/memutils.F90
mg_const.o: $(SRCROOT)atmos_param/strat_cloud/mg_const.F90 fms.o constants.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)atmos_param/strat_cloud/mg_const.F90
mg_drag.o: $(SRCROOT)atmos_param/mg_drag/mg_drag.F90 topography.o mpp.o fms.o fms_io.o constants.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)atmos_param/mg_drag/mg_drag.F90
mgrp_prscr_clds.o: $(SRCROOT)atmos_param/sea_esf_rad/mgrp_prscr_clds.F90 mpp.o fms.o constants.o rad_utilities.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)atmos_param/sea_esf_rad/mgrp_prscr_clds.F90
micro_mg.o: $(SRCROOT)atmos_param/strat_cloud/micro_mg.F90 gamma_mg.o strat_cloud_utilities.o mpp.o fms.o simple_pdf.o sat_vapor_pres.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)atmos_param/strat_cloud/micro_mg.F90
microphys_cloud.o: $(SRCROOT)atmos_param/sea_esf_rad/microphys_cloud.F90 mpp.o fms.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)atmos_param/sea_esf_rad/microphys_cloud.F90
microphys_rad.o: $(SRCROOT)atmos_param/sea_esf_rad/microphys_rad.F90 mpp.o fms.o constants.o time_manager.o diag_manager.o random_numbers.o cloud_generator.o rad_utilities.o longwave_params.o esfsw_parameters.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)atmos_param/sea_esf_rad/microphys_rad.F90
microphysics.o: $(SRCROOT)atmos_param/strat_cloud/microphysics.F90 fms.o mpp.o constants.o rotstayn_klein_mp.o morrison_gettelman_microp.o strat_cloud_utilities.o cldwat2m_micro.o micro_mg.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)atmos_param/strat_cloud/microphysics.F90
mixing_length.o: $(SRCROOT)atmos_param/clubb/CLUBB_core/mixing_length.F90 constants_clubb.o parameters_tunable.o parameters_model.o grid_class.o numerical_check.o saturation.o error_code.o model_flags.o clubb_precision.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)atmos_param/clubb/CLUBB_core/mixing_length.F90
mo_chem_utls.o: $(SRCROOT)atmos_shared/tracer_driver/tropchem/mo_chem_utls.F90 moz.mods.o m_tracname.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)atmos_shared/tracer_driver/tropchem/mo_chem_utls.F90
mo_chemdr.o: $(SRCROOT)atmos_shared/tracer_driver/tropchem/mo_chemdr.F90 moz.mods.o mo_photo.o mo_exp_slv.o mo_imp_slv.o mo_rodas_slv.o mo_usrrxt.o mo_setinv.o moz.subs.o mo_setsox.o mo_fphoto.o mo_chem_utls.o time_manager.o strat_chem_utilities.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)atmos_shared/tracer_driver/tropchem/mo_chemdr.F90
mo_chemini.o: $(SRCROOT)atmos_shared/tracer_driver/tropchem/mo_chemini.F90 mo_photo.o mo_chem_utls.o mo_usrrxt.o moz.mods.o mo_exp_slv.o mo_imp_slv.o mo_rodas_slv.o mo_read_sim_chm.o mo_fphoto.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)atmos_shared/tracer_driver/tropchem/mo_chemini.F90
mo_exp_slv.o: $(SRCROOT)atmos_shared/tracer_driver/tropchem/mo_exp_slv.F90 moz.mods.o mo_chem_utls.o moz.mat.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)atmos_shared/tracer_driver/tropchem/mo_exp_slv.F90
mo_fastjx.o: $(SRCROOT)atmos_shared/tracer_driver/tropchem/mo_fastjx.F90 fms.o mpp_io.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)atmos_shared/tracer_driver/tropchem/mo_fastjx.F90
mo_fphoto.o: $(SRCROOT)atmos_shared/tracer_driver/tropchem/mo_fphoto.F90 mpp.o fms.o tracer_manager.o field_manager.o mo_fastjx.o sat_vapor_pres.o mo_chem_utls.o moz.mods.o time_manager.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)atmos_shared/tracer_driver/tropchem/mo_fphoto.F90
mo_hook.o: $(SRCROOT)atmos_shared/tracer_driver/tropchem/mo_hook.F90 diag_manager.o time_manager.o constants.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)atmos_shared/tracer_driver/tropchem/mo_hook.F90
mo_imp_slv.o: $(SRCROOT)atmos_shared/tracer_driver/tropchem/mo_imp_slv.F90 moz.mods.o constants.o mo_chem_utls.o m_tracname.o moz.mat.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)atmos_shared/tracer_driver/tropchem/mo_imp_slv.F90
mo_jpl.o: $(SRCROOT)atmos_shared/tracer_driver/tropchem/mo_jpl.F90
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)atmos_shared/tracer_driver/tropchem/mo_jpl.F90
mo_photo.o: $(SRCROOT)atmos_shared/tracer_driver/tropchem/mo_photo.F90 fms.o mpp.o mpp_io.o time_manager.o constants.o mo_chem_utls.o moz.mods.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)atmos_shared/tracer_driver/tropchem/mo_photo.F90
mo_read_sim_chm.o: $(SRCROOT)atmos_shared/tracer_driver/tropchem/mo_read_sim_chm.F90 mpp.o mpp_io.o moz.mods.o m_tracname.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)atmos_shared/tracer_driver/tropchem/mo_read_sim_chm.F90
mo_rodas_slv.o: $(SRCROOT)atmos_shared/tracer_driver/tropchem/mo_rodas_slv.F90 moz.mods.o mo_chem_utls.o moz.mat.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)atmos_shared/tracer_driver/tropchem/mo_rodas_slv.F90
mo_setinv.o: $(SRCROOT)atmos_shared/tracer_driver/tropchem/mo_setinv.F90
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)atmos_shared/tracer_driver/tropchem/mo_setinv.F90
mo_setsox.o: $(SRCROOT)atmos_shared/tracer_driver/tropchem/mo_setsox.F90 mo_chem_utls.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)atmos_shared/tracer_driver/tropchem/mo_setsox.F90
mo_usrrxt.o: $(SRCROOT)atmos_shared/tracer_driver/tropchem/mo_usrrxt.F90 sat_vapor_pres.o constants.o strat_chem_utilities.o mo_chem_utls.o moz.mods.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)atmos_shared/tracer_driver/tropchem/mo_usrrxt.F90
model_flags.o: $(SRCROOT)atmos_param/clubb/CLUBB_core/model_flags.F90 constants_clubb.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)atmos_param/clubb/CLUBB_core/model_flags.F90
modis_simulator.o: $(SRCROOT)atmos_param/cosp/MODIS_simulator/modis_simulator.F90 /home/wenchang/tigress/CM2.5_tiger2_wy_killtc/exp/../src/atmos_param/cosp/cosp_defs.H cosp_types.o fms.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -I$(SRCROOT)atmos_param/cosp	$(SRCROOT)atmos_param/cosp/MODIS_simulator/modis_simulator.F90
moist_conv.o: $(SRCROOT)atmos_param/moist_conv/moist_conv.F90 mpp.o time_manager.o diag_manager.o sat_vapor_pres.o fms.o constants.o field_manager.o tracer_manager.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)atmos_param/moist_conv/moist_conv.F90
moist_processes.o: $(SRCROOT)atmos_param/moist_processes/moist_processes.F90 sat_vapor_pres.o time_manager.o diag_manager.o diag_data.o mpp.o fms.o field_manager.o tracer_manager.o constants.o betts_miller.o bm_massflux.o bm_omp.o donner_deep.o moist_conv.o lscale_cond.o uw_conv.o lin_cloud_microphys.o ras.o dry_adj.o strat_cloud.o detr_ice_num.o MG_microp_3D.o CLUBB_driver_SCM.o rh_clouds.o diag_cloud.o diag_integral.o cu_mo_trans.o mo_hook.o rad_utilities.o moist_processes_utils.o moistproc_kernels.o atmos_tracer_utilities.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)atmos_param/moist_processes/moist_processes.F90
moist_processes_utils.o: $(SRCROOT)atmos_param/moist_processes/moist_processes_utils.F90 sat_vapor_pres.o time_manager.o diag_manager.o constants.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)atmos_param/moist_processes/moist_processes_utils.F90
moistproc_kernels.o: $(SRCROOT)atmos_param/moist_processes/moistproc_kernels.F90 sat_vapor_pres.o time_manager.o diag_manager.o constants.o field_manager.o tracer_manager.o betts_miller.o bm_massflux.o bm_omp.o diag_cloud.o donner_deep.o moist_conv.o lscale_cond.o uw_conv.o lin_cloud_microphys.o ras.o strat_cloud.o rh_clouds.o cu_mo_trans.o atmos_tracer_utilities.o mo_hook.o rad_utilities.o moist_processes_utils.o detr_ice_num.o mpp.o MG_microp_3D.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)atmos_param/moist_processes/moistproc_kernels.F90
monin_obukhov.o: $(SRCROOT)atmos_param/monin_obukhov/monin_obukhov.F90 constants.o mpp.o fms.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)atmos_param/monin_obukhov/monin_obukhov.F90
monin_obukhov_kernel.o: $(SRCROOT)atmos_param/monin_obukhov/monin_obukhov_kernel.F90 /home/wenchang/tigress/CM2.5_tiger2_wy_killtc/exp/../src/shared/include/fms_platform.h /home/wenchang/tigress/CM2.5_tiger2_wy_killtc/exp/../src/atmos_param/monin_obukhov/monin_obukhov_interfaces.h
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -I$(SRCROOT)shared/include -I$(SRCROOT)atmos_param/monin_obukhov	$(SRCROOT)atmos_param/monin_obukhov/monin_obukhov_kernel.F90
mono_flux_limiter.o: $(SRCROOT)atmos_param/clubb/CLUBB_core/mono_flux_limiter.F90 grid_class.o constants_clubb.o clubb_precision.o error_code.o fill_holes.o stats_type.o stats_variables.o mean_adv.o lapack_wrap.o anl_erf.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)atmos_param/clubb/CLUBB_core/mono_flux_limiter.F90
morrison_gettelman_microp.o: $(SRCROOT)atmos_param/strat_cloud/morrison_gettelman_microp.F90 sat_vapor_pres.o polysvp.o constants.o mpp.o fms.o strat_cloud_utilities.o gamma_mg.o mg_const.o simple_pdf.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)atmos_param/strat_cloud/morrison_gettelman_microp.F90
mosaic.o: $(SRCROOT)shared/mosaic/mosaic.F90 mpp.o fms_io.o constants.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)shared/mosaic/mosaic.F90
mosaic_util.o: $(SRCROOT)shared/mosaic/mosaic_util.c /home/wenchang/tigress/CM2.5_tiger2_wy_killtc/exp/../src/shared/mosaic/mosaic_util.h /home/wenchang/tigress/CM2.5_tiger2_wy_killtc/exp/../src/shared/mosaic/constant.h
	$(CC) $(CPPDEFS) $(CPPFLAGS) $(CFLAGS) $(OTHERFLAGS) $(OTHER_CFLAGS) -c -I$(SRCROOT)shared/mosaic	$(SRCROOT)shared/mosaic/mosaic_util.c
moz.mat.o: $(SRCROOT)atmos_shared/tracer_driver/tropchem/moz.mat.F90 moz.mods.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)atmos_shared/tracer_driver/tropchem/moz.mat.F90
moz.mods.o: $(SRCROOT)atmos_shared/tracer_driver/tropchem/moz.mods.F90 mpp.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)atmos_shared/tracer_driver/tropchem/moz.mods.F90
moz.subs.o: $(SRCROOT)atmos_shared/tracer_driver/tropchem/moz.subs.F90 moz.mods.o mo_jpl.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)atmos_shared/tracer_driver/tropchem/moz.subs.F90
mpp.o: $(SRCROOT)shared/mpp/mpp.F90 /home/wenchang/tigress/CM2.5_tiger2_wy_killtc/exp/../src/shared/include/fms_platform.h /home/wenchang/tigress/CM2.5_tiger2_wy_killtc/exp/../src/shared/mpp/include/system_clock.h /home/wenchang/tigress/CM2.5_tiger2_wy_killtc/exp/../src/shared/mpp/include/mpp_util.inc /home/wenchang/tigress/CM2.5_tiger2_wy_killtc/exp/../src/shared/mpp/include/mpp_util_sma.inc /home/wenchang/tigress/CM2.5_tiger2_wy_killtc/exp/../src/shared/mpp/include/mpp_util_mpi.inc /home/wenchang/tigress/CM2.5_tiger2_wy_killtc/exp/../src/shared/mpp/include/mpp_util_nocomm.inc /home/wenchang/tigress/CM2.5_tiger2_wy_killtc/exp/../src/shared/mpp/include/mpp_error_a_a.h /home/wenchang/tigress/CM2.5_tiger2_wy_killtc/exp/../src/shared/mpp/include/mpp_error_a_s.h /home/wenchang/tigress/CM2.5_tiger2_wy_killtc/exp/../src/shared/mpp/include/mpp_error_s_a.h /home/wenchang/tigress/CM2.5_tiger2_wy_killtc/exp/../src/shared/mpp/include/mpp_error_s_s.h /home/wenchang/tigress/CM2.5_tiger2_wy_killtc/exp/../src/shared/mpp/include/mpp_comm.inc /home/wenchang/tigress/CM2.5_tiger2_wy_killtc/exp/../src/shared/mpp/include/mpp_comm_sma.inc /home/wenchang/tigress/CM2.5_tiger2_wy_killtc/exp/../src/shared/mpp/include/mpp_transmit_sma.h /home/wenchang/tigress/CM2.5_tiger2_wy_killtc/exp/../src/shared/mpp/include/mpp_transmit.inc /home/wenchang/tigress/CM2.5_tiger2_wy_killtc/exp/../src/shared/mpp/include/mpp_reduce_sma.h /home/wenchang/tigress/CM2.5_tiger2_wy_killtc/exp/../src/shared/mpp/include/mpp_sum_sma.h /home/wenchang/tigress/CM2.5_tiger2_wy_killtc/exp/../src/shared/mpp/include/mpp_sum.inc /home/wenchang/tigress/CM2.5_tiger2_wy_killtc/exp/../src/shared/mpp/include/mpp_comm_mpi.inc /home/wenchang/tigress/CM2.5_tiger2_wy_killtc/exp/../src/shared/mpp/include/mpp_transmit_mpi.h /home/wenchang/tigress/CM2.5_tiger2_wy_killtc/exp/../src/shared/mpp/include/mpp_reduce_mpi.h /home/wenchang/tigress/CM2.5_tiger2_wy_killtc/exp/../src/shared/mpp/include/mpp_sum_mpi.h /home/wenchang/tigress/CM2.5_tiger2_wy_killtc/exp/../src/shared/mpp/include/mpp_comm_nocomm.inc /home/wenchang/tigress/CM2.5_tiger2_wy_killtc/exp/../src/shared/mpp/include/mpp_transmit_nocomm.h /home/wenchang/tigress/CM2.5_tiger2_wy_killtc/exp/../src/shared/mpp/include/mpp_reduce_nocomm.h /home/wenchang/tigress/CM2.5_tiger2_wy_killtc/exp/../src/shared/mpp/include/mpp_sum_nocomm.h /home/wenchang/tigress/CM2.5_tiger2_wy_killtc/exp/../src/shared/mpp/include/mpp_chksum_int.h /home/wenchang/tigress/CM2.5_tiger2_wy_killtc/exp/../src/shared/mpp/include/mpp_chksum_scalar.h /home/wenchang/tigress/CM2.5_tiger2_wy_killtc/exp/../src/shared/mpp/include/mpp_chksum.h /home/wenchang/tigress/CM2.5_tiger2_wy_killtc/exp/../src/shared/mpp/include/mpp_gather.h /home/wenchang/tigress/CM2.5_tiger2_wy_killtc/exp/../src/shared/mpp/include/mpp_scatter.h mpp_parameter.o mpp_data.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -I$(SRCROOT)shared/include -I$(SRCROOT)shared/mpp/include	$(SRCROOT)shared/mpp/mpp.F90
mpp_data.o: $(SRCROOT)shared/mpp/mpp_data.F90 /home/wenchang/tigress/CM2.5_tiger2_wy_killtc/exp/../src/shared/include/fms_platform.h /home/wenchang/tigress/CM2.5_tiger2_wy_killtc/exp/../src/shared/mpp/include/mpp_data_sma.inc /home/wenchang/tigress/CM2.5_tiger2_wy_killtc/exp/../src/shared/mpp/include/mpp_data_mpi.inc /home/wenchang/tigress/CM2.5_tiger2_wy_killtc/exp/../src/shared/mpp/include/mpp_data_nocomm.inc mpp_parameter.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -I$(SRCROOT)shared/include -I$(SRCROOT)shared/mpp/include	$(SRCROOT)shared/mpp/mpp_data.F90
mpp_domains.o: $(SRCROOT)shared/mpp/mpp_domains.F90 /home/wenchang/tigress/CM2.5_tiger2_wy_killtc/exp/../src/shared/include/fms_platform.h /home/wenchang/tigress/CM2.5_tiger2_wy_killtc/exp/../src/shared/mpp/include/mpp_define_nest_domains.inc /home/wenchang/tigress/CM2.5_tiger2_wy_killtc/exp/../src/shared/mpp/include/mpp_domains_util.inc /home/wenchang/tigress/CM2.5_tiger2_wy_killtc/exp/../src/shared/mpp/include/mpp_domains_comm.inc /home/wenchang/tigress/CM2.5_tiger2_wy_killtc/exp/../src/shared/mpp/include/mpp_domains_define.inc /home/wenchang/tigress/CM2.5_tiger2_wy_killtc/exp/../src/shared/mpp/include/mpp_domains_misc.inc /home/wenchang/tigress/CM2.5_tiger2_wy_killtc/exp/../src/shared/mpp/include/mpp_update_domains2D.h /home/wenchang/tigress/CM2.5_tiger2_wy_killtc/exp/../src/shared/mpp/include/mpp_update_domains2D_nonblock.h /home/wenchang/tigress/CM2.5_tiger2_wy_killtc/exp/../src/shared/mpp/include/mpp_do_update_nonblock.h /home/wenchang/tigress/CM2.5_tiger2_wy_killtc/exp/../src/shared/mpp/include/mpp_do_updateV_nonblock.h /home/wenchang/tigress/CM2.5_tiger2_wy_killtc/exp/../src/shared/mpp/include/mpp_do_update.h /home/wenchang/tigress/CM2.5_tiger2_wy_killtc/exp/../src/shared/mpp/include/mpp_do_updateV.h /home/wenchang/tigress/CM2.5_tiger2_wy_killtc/exp/../src/shared/mpp/include/mpp_do_check.h /home/wenchang/tigress/CM2.5_tiger2_wy_killtc/exp/../src/shared/mpp/include/mpp_do_checkV.h /home/wenchang/tigress/CM2.5_tiger2_wy_killtc/exp/../src/shared/mpp/include/mpp_update_nest_domains.h /home/wenchang/tigress/CM2.5_tiger2_wy_killtc/exp/../src/shared/mpp/include/mpp_do_update_nest.h /home/wenchang/tigress/CM2.5_tiger2_wy_killtc/exp/../src/shared/mpp/include/mpp_do_redistribute.h /home/wenchang/tigress/CM2.5_tiger2_wy_killtc/exp/../src/shared/mpp/include/mpp_get_boundary.h /home/wenchang/tigress/CM2.5_tiger2_wy_killtc/exp/../src/shared/mpp/include/mpp_do_get_boundary.h /home/wenchang/tigress/CM2.5_tiger2_wy_killtc/exp/../src/shared/mpp/include/mpp_group_update.h /home/wenchang/tigress/CM2.5_tiger2_wy_killtc/exp/../src/shared/mpp/include/mpp_domains_reduce.inc /home/wenchang/tigress/CM2.5_tiger2_wy_killtc/exp/../src/shared/mpp/include/mpp_global_reduce.h /home/wenchang/tigress/CM2.5_tiger2_wy_killtc/exp/../src/shared/mpp/include/mpp_global_sum.h /home/wenchang/tigress/CM2.5_tiger2_wy_killtc/exp/../src/shared/mpp/include/mpp_global_sum_tl.h /home/wenchang/tigress/CM2.5_tiger2_wy_killtc/exp/../src/shared/mpp/include/mpp_global_field.h /home/wenchang/tigress/CM2.5_tiger2_wy_killtc/exp/../src/shared/mpp/include/mpp_do_global_field.h mpp_parameter.o mpp_data.o mpp.o mpp_memutils.o mpp_pset.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -I$(SRCROOT)shared/include -I$(SRCROOT)shared/mpp/include	$(SRCROOT)shared/mpp/mpp_domains.F90
mpp_io.o: $(SRCROOT)shared/mpp/mpp_io.F90 /home/wenchang/tigress/CM2.5_tiger2_wy_killtc/exp/../src/shared/include/fms_platform.h /home/wenchang/tigress/CM2.5_tiger2_wy_killtc/exp/../src/shared/mpp/include/mpp_io_util.inc /home/wenchang/tigress/CM2.5_tiger2_wy_killtc/exp/../src/shared/mpp/include/mpp_io_misc.inc /home/wenchang/tigress/CM2.5_tiger2_wy_killtc/exp/../src/shared/mpp/include/mpp_io_connect.inc /home/wenchang/tigress/CM2.5_tiger2_wy_killtc/exp/../src/shared/mpp/include/mpp_io_read.inc /home/wenchang/tigress/CM2.5_tiger2_wy_killtc/exp/../src/shared/mpp/include/mpp_read_2Ddecomp.h /home/wenchang/tigress/CM2.5_tiger2_wy_killtc/exp/../src/shared/mpp/include/mpp_read_compressed.h /home/wenchang/tigress/CM2.5_tiger2_wy_killtc/exp/../src/shared/mpp/include/mpp_io_write.inc /home/wenchang/tigress/CM2.5_tiger2_wy_killtc/exp/../src/shared/mpp/include/mpp_write_2Ddecomp.h /home/wenchang/tigress/CM2.5_tiger2_wy_killtc/exp/../src/shared/mpp/include/mpp_write_compressed.h /home/wenchang/tigress/CM2.5_tiger2_wy_killtc/exp/../src/shared/mpp/include/mpp_write.h mpp_parameter.o mpp.o mpp_domains.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -I$(SRCROOT)shared/include -I$(SRCROOT)shared/mpp/include	$(SRCROOT)shared/mpp/mpp_io.F90
mpp_memutils.o: $(SRCROOT)shared/mpp/mpp_memutils.F90 mpp.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)shared/mpp/mpp_memutils.F90
mpp_parameter.o: $(SRCROOT)shared/mpp/mpp_parameter.F90 /home/wenchang/tigress/CM2.5_tiger2_wy_killtc/exp/../src/shared/include/fms_platform.h
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -I$(SRCROOT)shared/include	$(SRCROOT)shared/mpp/mpp_parameter.F90
mpp_pset.o: $(SRCROOT)shared/mpp/mpp_pset.F90 /home/wenchang/tigress/CM2.5_tiger2_wy_killtc/exp/../src/shared/include/fms_platform.h mpp.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -I$(SRCROOT)shared/include	$(SRCROOT)shared/mpp/mpp_pset.F90
mpp_utilities.o: $(SRCROOT)shared/mpp/mpp_utilities.F90 mpp.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)shared/mpp/mpp_utilities.F90
mrgrnk.o: $(SRCROOT)atmos_param/cosp/quickbeam/mrgrnk.F90 /home/wenchang/tigress/CM2.5_tiger2_wy_killtc/exp/../src/atmos_param/cosp/cosp_defs.H
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -I$(SRCROOT)atmos_param/cosp	$(SRCROOT)atmos_param/cosp/quickbeam/mrgrnk.F90
mt95.o: $(SRCROOT)atmos_param/clubb/CLUBB_core/mt95.f90
	$(FC) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)atmos_param/clubb/CLUBB_core/mt95.f90
my25_turb.o: $(SRCROOT)atmos_param/my25_turb/my25_turb.F90 mpp.o fms.o fms_io.o tridiagonal.o constants.o monin_obukhov.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)atmos_param/my25_turb/my25_turb.F90
nc_cond.o: $(SRCROOT)atmos_param/strat_cloud/nc_cond.F90 mpp.o fms.o constants.o betaDistribution.o strat_cloud_utilities.o polysvp.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)atmos_param/strat_cloud/nc_cond.F90
nf_utils.o: $(SRCROOT)land_lad2/shared/nf_utils/nf_utils.F90 nfu.o nfc.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)land_lad2/shared/nf_utils/nf_utils.F90
nfc.o: $(SRCROOT)land_lad2/shared/nf_utils/nfc.F90 /home/wenchang/tigress/CM2.5_tiger2_wy_killtc/exp/../src/land_lad2/shared/nf_utils/getput_compressed.inc nfu.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -I$(SRCROOT)land_lad2/shared/nf_utils	$(SRCROOT)land_lad2/shared/nf_utils/nfc.F90
nfu.o: $(SRCROOT)land_lad2/shared/nf_utils/nfu.F90 /home/wenchang/tigress/CM2.5_tiger2_wy_killtc/exp/../src/land_lad2/shared/nf_utils/getput.inc
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -I$(SRCROOT)land_lad2/shared/nf_utils	$(SRCROOT)land_lad2/shared/nf_utils/nfu.F90
nh_core.o: $(SRCROOT)atmos_cubed_sphere/model/null/nh_core.F90 fms.o fv_arrays.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)atmos_cubed_sphere/model/null/nh_core.F90
nonfms_donner.o: $(SRCROOT)atmos_param/donner_deep/nonfms_donner.F90 /home/wenchang/tigress/CM2.5_tiger2_wy_killtc/exp/../src/atmos_param/donner_deep/donner_nml.h sat_vapor_pres_k.o donner_types.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -I$(SRCROOT)atmos_param/donner_deep	$(SRCROOT)atmos_param/donner_deep/nonfms_donner.F90
nsclock.o: $(SRCROOT)shared/mpp/nsclock.c
	$(CC) $(CPPDEFS) $(CPPFLAGS) $(CFLAGS) $(OTHERFLAGS) $(OTHER_CFLAGS) -c	$(SRCROOT)shared/mpp/nsclock.c
numerical_check.o: $(SRCROOT)atmos_param/clubb/CLUBB_core/numerical_check.F90 grid_class.o clubb_precision.o parameters_model.o pdf_parameter_module.o stats_variables.o constants_clubb.o variables_diagnostic_module.o variables_prognostic_module.o parameters_microphys.o error_code.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)atmos_param/clubb/CLUBB_core/numerical_check.F90
ocean_adv_vel_diag.o: $(SRCROOT)mom5/ocean_diag/ocean_adv_vel_diag.F90 /home/wenchang/tigress/CM2.5_tiger2_wy_killtc/exp/../src/mom5/ocean_core/ocean_memory.h constants.o diag_manager.o fms.o mpp.o time_manager.o ocean_domains.o ocean_operators.o ocean_parameters.o ocean_types.o ocean_util.o ocean_workspace.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -I$(SRCROOT)mom5/ocean_core	$(SRCROOT)mom5/ocean_diag/ocean_adv_vel_diag.F90
ocean_advection_velocity.o: $(SRCROOT)mom5/ocean_core/ocean_advection_velocity.F90 /home/wenchang/tigress/CM2.5_tiger2_wy_killtc/exp/../src/mom5/ocean_core/ocean_memory.h axis_utils.o constants.o diag_manager.o fms.o fms_io.o mpp_domains.o mpp.o ocean_domains.o ocean_operators.o ocean_parameters.o ocean_types.o ocean_workspace.o ocean_obc.o ocean_util.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -I$(SRCROOT)mom5/ocean_core	$(SRCROOT)mom5/ocean_core/ocean_advection_velocity.F90
ocean_age_tracer.o: $(SRCROOT)mom5/ocean_tracers/ocean_age_tracer.F90 field_manager.o mpp.o time_manager.o ocean_tpm_util.o fm_util.o ocean_types.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)mom5/ocean_tracers/ocean_age_tracer.F90
ocean_albedo.o: $(SRCROOT)ice_param/ocean_albedo.F90 mpp.o fms.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)ice_param/ocean_albedo.F90
ocean_barotropic.o: $(SRCROOT)mom5/ocean_core/ocean_barotropic.F90 /home/wenchang/tigress/CM2.5_tiger2_wy_killtc/exp/../src/mom5/ocean_core/ocean_memory.h constants.o diag_manager.o fms.o fms_io.o mpp_domains.o mpp.o time_manager.o time_interp_external.o ocean_bih_friction.o ocean_domains.o ocean_lap_friction.o ocean_obc_barotrop.o ocean_obc.o ocean_operators.o ocean_parameters.o ocean_types.o ocean_util.o ocean_workspace.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -I$(SRCROOT)mom5/ocean_core	$(SRCROOT)mom5/ocean_core/ocean_barotropic.F90
ocean_bbc.o: $(SRCROOT)mom5/ocean_core/ocean_bbc.F90 /home/wenchang/tigress/CM2.5_tiger2_wy_killtc/exp/../src/mom5/ocean_core/ocean_memory.h constants.o diag_manager.o fms.o mpp.o mpp_domains.o ocean_domains.o ocean_parameters.o ocean_types.o ocean_workspace.o ocean_wave_type.o ocean_wave.o ocean_util.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -I$(SRCROOT)mom5/ocean_core	$(SRCROOT)mom5/ocean_core/ocean_bbc.F90
ocean_bgc_restore.o: $(SRCROOT)mom5/ocean_bgc/ocean_bgc_restore.F90 /home/wenchang/tigress/CM2.5_tiger2_wy_killtc/exp/../src/shared/include/fms_platform.h time_manager.o diag_manager.o field_manager.o fms.o fms_io.o mpp.o time_interp_external.o mpp_domains.o constants.o ocean_tpm_util.o fm_util.o coupler_types.o ocean_types.o ocmip2_co2calc.o ocean_util.o atmos_ocean_fluxes.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -I$(SRCROOT)shared/include	$(SRCROOT)mom5/ocean_bgc/ocean_bgc_restore.F90
ocean_bih_friction.o: $(SRCROOT)mom5/ocean_param/lateral/ocean_bih_friction.F90 constants.o diag_manager.o fms.o fms_io.o mpp_domains.o mpp.o ocean_bihcgrid_friction.o ocean_bihcst_friction.o ocean_bihgen_friction.o ocean_domains.o ocean_operators.o ocean_parameters.o ocean_types.o ocean_util.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)mom5/ocean_param/lateral/ocean_bih_friction.F90
ocean_bih_tracer.o: $(SRCROOT)mom5/ocean_param/lateral/ocean_bih_tracer.F90 constants.o diag_manager.o fms.o mpp_domains.o mpp.o ocean_domains.o ocean_operators.o ocean_obc.o ocean_parameters.o ocean_types.o ocean_workspace.o ocean_util.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)mom5/ocean_param/lateral/ocean_bih_tracer.F90
ocean_bihcgrid_friction.o: $(SRCROOT)mom5/ocean_param/lateral/ocean_bihcgrid_friction.F90 /home/wenchang/tigress/CM2.5_tiger2_wy_killtc/exp/../src/mom5/ocean_core/ocean_memory.h constants.o diag_manager.o fms.o mpp_domains.o mpp.o ocean_domains.o ocean_obc.o ocean_operators.o ocean_parameters.o ocean_types.o ocean_util.o ocean_workspace.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -I$(SRCROOT)mom5/ocean_core	$(SRCROOT)mom5/ocean_param/lateral/ocean_bihcgrid_friction.F90
ocean_bihcst_friction.o: $(SRCROOT)mom5/ocean_param/lateral/ocean_bihcst_friction.F90 constants.o diag_manager.o fms.o mpp_domains.o mpp.o ocean_domains.o ocean_obc.o ocean_operators.o ocean_parameters.o ocean_types.o ocean_util.o ocean_workspace.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)mom5/ocean_param/lateral/ocean_bihcst_friction.F90
ocean_bihgen_friction.o: $(SRCROOT)mom5/ocean_param/lateral/ocean_bihgen_friction.F90 /home/wenchang/tigress/CM2.5_tiger2_wy_killtc/exp/../src/mom5/ocean_core/ocean_memory.h constants.o diag_manager.o fms.o mpp_domains.o mpp.o ocean_domains.o ocean_obc.o ocean_operators.o ocean_parameters.o ocean_types.o ocean_util.o ocean_workspace.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -I$(SRCROOT)mom5/ocean_core	$(SRCROOT)mom5/ocean_param/lateral/ocean_bihgen_friction.F90
ocean_blob.o: $(SRCROOT)mom5/ocean_blobs/ocean_blob.F90 constants.o diag_manager.o fms.o fms_io.o mpp_domains.o mpp.o ocean_blob_diag.o ocean_blob_dynamic_bottom.o ocean_blob_dynamic_free.o ocean_blob_static_free.o ocean_blob_static_bottom.o ocean_blob_util.o ocean_density.o ocean_domains.o ocean_parameters.o ocean_types.o ocean_util.o ocean_workspace.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)mom5/ocean_blobs/ocean_blob.F90
ocean_blob_diag.o: $(SRCROOT)mom5/ocean_blobs/ocean_blob_diag.F90 constants.o fms.o mpp.o ocean_blob_util.o ocean_parameters.o ocean_types.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)mom5/ocean_blobs/ocean_blob_diag.F90
ocean_blob_dynamic_bottom.o: $(SRCROOT)mom5/ocean_blobs/ocean_blob_dynamic_bottom.F90 constants.o diag_manager.o fms.o mpp.o mpp_domains.o ocean_blob_util.o ocean_density.o ocean_parameters.o ocean_types.o ocean_util.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)mom5/ocean_blobs/ocean_blob_dynamic_bottom.F90
ocean_blob_dynamic_free.o: $(SRCROOT)mom5/ocean_blobs/ocean_blob_dynamic_free.F90 constants.o diag_manager.o fms.o mpp.o mpp_domains.o ocean_blob_util.o ocean_density.o ocean_parameters.o ocean_types.o ocean_util.o ocean_workspace.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)mom5/ocean_blobs/ocean_blob_dynamic_free.F90
ocean_blob_static_bottom.o: $(SRCROOT)mom5/ocean_blobs/ocean_blob_static_bottom.F90 fms.o mpp_domains.o mpp.o ocean_blob_util.o ocean_density.o ocean_parameters.o ocean_types.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)mom5/ocean_blobs/ocean_blob_static_bottom.F90
ocean_blob_static_free.o: $(SRCROOT)mom5/ocean_blobs/ocean_blob_static_free.F90 fms.o ocean_blob_util.o ocean_density.o ocean_types.o ocean_workspace.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)mom5/ocean_blobs/ocean_blob_static_free.F90
ocean_blob_util.o: $(SRCROOT)mom5/ocean_blobs/ocean_blob_util.F90 constants.o fms.o mpp_domains.o mpp.o grid.o ocean_parameters.o ocean_types.o ocean_workspace.o ocean_util.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)mom5/ocean_blobs/ocean_blob_util.F90
ocean_convect.o: $(SRCROOT)mom5/ocean_param/vertical/ocean_convect.F90 constants.o diag_manager.o fms.o mpp.o ocean_density.o ocean_domains.o ocean_parameters.o ocean_types.o ocean_workspace.o ocean_util.o ocean_tracer_util.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)mom5/ocean_param/vertical/ocean_convect.F90
ocean_coriolis.o: $(SRCROOT)mom5/ocean_core/ocean_coriolis.F90 /home/wenchang/tigress/CM2.5_tiger2_wy_killtc/exp/../src/mom5/ocean_core/ocean_memory.h constants.o diag_manager.o fms.o mpp.o ocean_domains.o ocean_parameters.o ocean_types.o ocean_util.o ocean_workspace.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -I$(SRCROOT)mom5/ocean_core	$(SRCROOT)mom5/ocean_core/ocean_coriolis.F90
ocean_density.o: $(SRCROOT)mom5/ocean_core/ocean_density.F90 /home/wenchang/tigress/CM2.5_tiger2_wy_killtc/exp/../src/shared/include/fms_platform.h /home/wenchang/tigress/CM2.5_tiger2_wy_killtc/exp/../src/mom5/ocean_core/ocean_memory.h constants.o diag_manager.o fms.o fms_io.o mpp_domains.o mpp.o platform.o time_manager.o field_manager.o ocean_domains.o ocean_operators.o ocean_parameters.o ocean_pressure.o ocean_types.o ocean_util.o ocean_workspace.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -I$(SRCROOT)shared/include -I$(SRCROOT)mom5/ocean_core	$(SRCROOT)mom5/ocean_core/ocean_density.F90
ocean_diagnostics.o: $(SRCROOT)mom5/ocean_diag/ocean_diagnostics.F90 /home/wenchang/tigress/CM2.5_tiger2_wy_killtc/exp/../src/mom5/ocean_core/ocean_memory.h diag_manager.o fms.o mpp.o time_manager.o ocean_adv_vel_diag.o ocean_domains.o ocean_tracer_diag.o ocean_types.o ocean_velocity_diag.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -I$(SRCROOT)mom5/ocean_core	$(SRCROOT)mom5/ocean_diag/ocean_diagnostics.F90
ocean_domains.o: $(SRCROOT)mom5/ocean_core/ocean_domains.F90 fms.o mpp_domains.o mpp.o ocean_types.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)mom5/ocean_core/ocean_domains.F90
ocean_drifters.o: $(SRCROOT)mom5/ocean_diag/ocean_drifters.F90 drifters.o fms.o mpp_domains.o mpp.o time_manager.o ocean_domains.o ocean_parameters.o ocean_types.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)mom5/ocean_diag/ocean_drifters.F90
ocean_form_drag.o: $(SRCROOT)mom5/ocean_param/vertical/ocean_form_drag.F90 constants.o diag_manager.o fms.o mpp.o mpp_domains.o ocean_domains.o ocean_parameters.o ocean_types.o ocean_workspace.o ocean_util.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)mom5/ocean_param/vertical/ocean_form_drag.F90
ocean_frazil.o: $(SRCROOT)mom5/ocean_tracers/ocean_frazil.F90 /home/wenchang/tigress/CM2.5_tiger2_wy_killtc/exp/../src/mom5/ocean_core/ocean_memory.h constants.o diag_manager.o fms.o mpp.o mpp_domains.o ocean_domains.o ocean_parameters.o ocean_tpm_util.o ocean_types.o ocean_workspace.o ocean_util.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -I$(SRCROOT)mom5/ocean_core	$(SRCROOT)mom5/ocean_tracers/ocean_frazil.F90
ocean_generic_tracer.o: $(SRCROOT)mom5/ocean_bgc/ocean_generic_tracer.F90 /home/wenchang/tigress/CM2.5_tiger2_wy_killtc/exp/../src/shared/include/fms_platform.h ocean_tpm_util.o ocean_types.o ocean_parameters.o mpp.o field_manager.o generic_tracer.o generic_tracer_utils.o coupler_types.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -I$(SRCROOT)shared/include	$(SRCROOT)mom5/ocean_bgc/ocean_generic_tracer.F90
ocean_grids.o: $(SRCROOT)mom5/ocean_core/ocean_grids.F90 /home/wenchang/tigress/CM2.5_tiger2_wy_killtc/exp/../src/mom5/ocean_core/ocean_memory.h constants.o diag_manager.o fms.o mpp_domains.o mpp.o mosaic.o ocean_domains.o ocean_parameters.o ocean_types.o ocean_workspace.o ocean_util.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -I$(SRCROOT)mom5/ocean_core	$(SRCROOT)mom5/ocean_core/ocean_grids.F90
ocean_ibgc.o: $(SRCROOT)mom5/ocean_bgc/ocean_ibgc.F90 /home/wenchang/tigress/CM2.5_tiger2_wy_killtc/exp/../src/shared/include/fms_platform.h diag_manager.o field_manager.o mpp.o fms.o fms_io.o time_manager.o time_interp_external.o mpp_domains.o constants.o ocean_tpm_util.o fm_util.o ocean_types.o ocmip2_co2calc.o coupler_types.o atmos_ocean_fluxes.o ocean_util.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -I$(SRCROOT)shared/include	$(SRCROOT)mom5/ocean_bgc/ocean_ibgc.F90
ocean_increment_eta.o: $(SRCROOT)mom5/ocean_param/sources/ocean_increment_eta.F90 /home/wenchang/tigress/CM2.5_tiger2_wy_killtc/exp/../src/mom5/ocean_core/ocean_memory.h diag_manager.o fms.o mpp.o time_interp_external.o time_manager.o ocean_domains.o ocean_parameters.o ocean_types.o ocean_workspace.o ocean_util.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -I$(SRCROOT)mom5/ocean_core	$(SRCROOT)mom5/ocean_param/sources/ocean_increment_eta.F90
ocean_increment_tracer.o: $(SRCROOT)mom5/ocean_param/sources/ocean_increment_tracer.F90 /home/wenchang/tigress/CM2.5_tiger2_wy_killtc/exp/../src/mom5/ocean_core/ocean_memory.h diag_manager.o fms.o mpp.o time_interp_external.o time_manager.o ocean_domains.o ocean_parameters.o ocean_types.o ocean_workspace.o ocean_util.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -I$(SRCROOT)mom5/ocean_core	$(SRCROOT)mom5/ocean_param/sources/ocean_increment_tracer.F90
ocean_increment_velocity.o: $(SRCROOT)mom5/ocean_param/sources/ocean_increment_velocity.F90 /home/wenchang/tigress/CM2.5_tiger2_wy_killtc/exp/../src/mom5/ocean_core/ocean_memory.h diag_manager.o fms.o mpp.o time_interp_external.o time_manager.o ocean_domains.o ocean_parameters.o ocean_types.o ocean_workspace.o ocean_util.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -I$(SRCROOT)mom5/ocean_core	$(SRCROOT)mom5/ocean_param/sources/ocean_increment_velocity.F90
ocean_lap_friction.o: $(SRCROOT)mom5/ocean_param/lateral/ocean_lap_friction.F90 diag_manager.o fms.o fms_io.o mpp_domains.o mpp.o ocean_domains.o ocean_lapcst_friction.o ocean_lapgen_friction.o ocean_lapcgrid_friction.o ocean_operators.o ocean_parameters.o ocean_types.o ocean_util.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)mom5/ocean_param/lateral/ocean_lap_friction.F90
ocean_lap_tracer.o: $(SRCROOT)mom5/ocean_param/lateral/ocean_lap_tracer.F90 constants.o diag_manager.o fms.o mpp_domains.o mpp.o ocean_domains.o ocean_obc.o ocean_operators.o ocean_parameters.o ocean_types.o ocean_workspace.o ocean_util.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)mom5/ocean_param/lateral/ocean_lap_tracer.F90
ocean_lapcgrid_friction.o: $(SRCROOT)mom5/ocean_param/lateral/ocean_lapcgrid_friction.F90 /home/wenchang/tigress/CM2.5_tiger2_wy_killtc/exp/../src/mom5/ocean_core/ocean_memory.h constants.o diag_manager.o fms.o mpp_domains.o mpp.o ocean_domains.o ocean_obc.o ocean_operators.o ocean_parameters.o ocean_types.o ocean_util.o ocean_workspace.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -I$(SRCROOT)mom5/ocean_core	$(SRCROOT)mom5/ocean_param/lateral/ocean_lapcgrid_friction.F90
ocean_lapcst_friction.o: $(SRCROOT)mom5/ocean_param/lateral/ocean_lapcst_friction.F90 constants.o diag_manager.o fms.o mpp_domains.o mpp.o ocean_domains.o ocean_obc.o ocean_operators.o ocean_parameters.o ocean_types.o ocean_util.o ocean_workspace.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)mom5/ocean_param/lateral/ocean_lapcst_friction.F90
ocean_lapgen_friction.o: $(SRCROOT)mom5/ocean_param/lateral/ocean_lapgen_friction.F90 /home/wenchang/tigress/CM2.5_tiger2_wy_killtc/exp/../src/mom5/ocean_core/ocean_memory.h constants.o diag_manager.o fms.o mpp_domains.o mpp.o ocean_domains.o ocean_obc.o ocean_operators.o ocean_parameters.o ocean_types.o ocean_util.o ocean_workspace.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -I$(SRCROOT)mom5/ocean_core	$(SRCROOT)mom5/ocean_param/lateral/ocean_lapgen_friction.F90
ocean_mixdownslope.o: $(SRCROOT)mom5/ocean_param/lateral/ocean_mixdownslope.F90 /home/wenchang/tigress/CM2.5_tiger2_wy_killtc/exp/../src/mom5/ocean_core/ocean_memory.h constants.o diag_manager.o fms.o mpp_domains.o mpp.o ocean_density.o ocean_domains.o ocean_parameters.o ocean_types.o ocean_util.o ocean_tracer_util.o ocean_workspace.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -I$(SRCROOT)mom5/ocean_core	$(SRCROOT)mom5/ocean_param/lateral/ocean_mixdownslope.F90
ocean_model.o: $(SRCROOT)mom5/ocean_core/ocean_model.F90 /home/wenchang/tigress/CM2.5_tiger2_wy_killtc/exp/../src/mom5/ocean_core/ocean_memory.h fms.o fms_io.o mpp_domains.o mpp.o stock_constants.o time_interp_external.o time_manager.o ocean_advection_velocity.o ocean_barotropic.o ocean_bbc.o ocean_bih_friction.o ocean_blob.o ocean_convect.o ocean_coriolis.o ocean_density.o ocean_diagnostics.o ocean_domains.o ocean_form_drag.o ocean_grids.o ocean_increment_eta.o ocean_increment_tracer.o ocean_increment_velocity.o ocean_lap_tracer.o ocean_bih_tracer.o ocean_lap_friction.o ocean_mixdownslope.o ocean_momentum_source.o ocean_nphysics.o ocean_nphysics_new.o ocean_obc.o ocean_operators.o ocean_overexchange.o ocean_overflow.o ocean_overflow_OFP.o ocean_passive.o ocean_pressure.o ocean_rivermix.o ocean_riverspread.o ocean_parameters.o ocean_sbc.o ocean_shortwave.o ocean_sigma_transport.o ocean_sponges_eta.o ocean_sponges_tracer.o ocean_sponges_velocity.o ocean_submesoscale.o ocean_tempsalt.o ocean_thickness.o ocean_topog.o ocean_tracer_advect.o ocean_tracer.o ocean_tracer_diag.o ocean_tracer_util.o ocean_tpm.o ocean_types.o ocean_util.o ocean_velocity_advect.o ocean_velocity_diag.o ocean_velocity.o ocean_vert_mix.o ocean_vert_gotm.o ocean_workspace.o ocean_xlandinsert.o ocean_xlandmix.o ocean_drifters.o ocean_wave_type.o ocean_wave.o oda_types.o oda_driver.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -I$(SRCROOT)mom5/ocean_core	$(SRCROOT)mom5/ocean_core/ocean_model.F90
ocean_momentum_source.o: $(SRCROOT)mom5/ocean_param/sources/ocean_momentum_source.F90 constants.o diag_manager.o field_manager.o fms.o mpp_domains.o mpp.o ocean_domains.o ocean_parameters.o ocean_types.o ocean_workspace.o ocean_util.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)mom5/ocean_param/sources/ocean_momentum_source.F90
ocean_nphysics.o: $(SRCROOT)mom5/ocean_param/neutral/ocean_nphysics.F90 /home/wenchang/tigress/CM2.5_tiger2_wy_killtc/exp/../src/mom5/ocean_core/ocean_memory.h constants.o diag_manager.o fms.o mpp.o mpp_domains.o ocean_domains.o ocean_nphysics_util.o ocean_nphysicsA.o ocean_nphysicsB.o ocean_nphysicsC.o ocean_parameters.o ocean_types.o ocean_util.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -I$(SRCROOT)mom5/ocean_core	$(SRCROOT)mom5/ocean_param/neutral/ocean_nphysics.F90
ocean_nphysicsA.o: $(SRCROOT)mom5/ocean_param/neutral/ocean_nphysicsA.F90 /home/wenchang/tigress/CM2.5_tiger2_wy_killtc/exp/../src/mom5/ocean_core/ocean_memory.h constants.o diag_manager.o fms.o mpp_domains.o mpp.o time_manager.o ocean_domains.o ocean_nphysics_util.o ocean_operators.o ocean_parameters.o ocean_sigma_transport.o ocean_types.o ocean_util.o ocean_workspace.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -I$(SRCROOT)mom5/ocean_core	$(SRCROOT)mom5/ocean_param/neutral/ocean_nphysicsA.F90
ocean_nphysicsB.o: $(SRCROOT)mom5/ocean_param/neutral/ocean_nphysicsB.F90 /home/wenchang/tigress/CM2.5_tiger2_wy_killtc/exp/../src/mom5/ocean_core/ocean_memory.h constants.o diag_manager.o fms.o fms_io.o mpp_domains.o mpp.o time_manager.o ocean_domains.o ocean_nphysics_util.o ocean_operators.o ocean_parameters.o ocean_sigma_transport.o ocean_types.o ocean_util.o ocean_workspace.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -I$(SRCROOT)mom5/ocean_core	$(SRCROOT)mom5/ocean_param/neutral/ocean_nphysicsB.F90
ocean_nphysicsC.o: $(SRCROOT)mom5/ocean_param/neutral/ocean_nphysicsC.F90 /home/wenchang/tigress/CM2.5_tiger2_wy_killtc/exp/../src/mom5/ocean_core/ocean_memory.h constants.o diag_manager.o fms.o fms_io.o mpp_domains.o mpp.o time_manager.o ocean_density.o ocean_domains.o ocean_nphysics_util.o ocean_operators.o ocean_parameters.o ocean_tracer_diag.o ocean_types.o ocean_util.o ocean_workspace.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -I$(SRCROOT)mom5/ocean_core	$(SRCROOT)mom5/ocean_param/neutral/ocean_nphysicsC.F90
ocean_nphysics_diff.o: $(SRCROOT)mom5/ocean_param/neutral/ocean_nphysics_diff.F90 /home/wenchang/tigress/CM2.5_tiger2_wy_killtc/exp/../src/mom5/ocean_core/ocean_memory.h constants.o fms.o mpp_domains.o mpp.o ocean_domains.o ocean_nphysics_util_new.o ocean_operators.o ocean_parameters.o ocean_tracer_diag.o ocean_types.o ocean_util.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -I$(SRCROOT)mom5/ocean_core	$(SRCROOT)mom5/ocean_param/neutral/ocean_nphysics_diff.F90
ocean_nphysics_flux.o: $(SRCROOT)mom5/ocean_param/neutral/ocean_nphysics_flux.F90 /home/wenchang/tigress/CM2.5_tiger2_wy_killtc/exp/../src/mom5/ocean_core/ocean_memory.h constants.o fms.o mpp_domains.o mpp.o ocean_domains.o ocean_nphysics_util_new.o ocean_operators.o ocean_types.o ocean_util.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -I$(SRCROOT)mom5/ocean_core	$(SRCROOT)mom5/ocean_param/neutral/ocean_nphysics_flux.F90
ocean_nphysics_new.o: $(SRCROOT)mom5/ocean_param/neutral/ocean_nphysics_new.F90 /home/wenchang/tigress/CM2.5_tiger2_wy_killtc/exp/../src/mom5/ocean_core/ocean_memory.h constants.o fms.o fms_io.o mpp.o mpp_domains.o ocean_domains.o ocean_nphysics_util_new.o ocean_nphysics_flux.o ocean_nphysics_tensor.o ocean_nphysics_skew.o ocean_nphysics_diff.o ocean_operators.o ocean_parameters.o ocean_types.o ocean_util.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -I$(SRCROOT)mom5/ocean_core	$(SRCROOT)mom5/ocean_param/neutral/ocean_nphysics_new.F90
ocean_nphysics_skew.o: $(SRCROOT)mom5/ocean_param/neutral/ocean_nphysics_skew.F90 /home/wenchang/tigress/CM2.5_tiger2_wy_killtc/exp/../src/mom5/ocean_core/ocean_memory.h constants.o fms.o mpp.o ocean_domains.o ocean_nphysics_util_new.o ocean_parameters.o ocean_types.o ocean_util.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -I$(SRCROOT)mom5/ocean_core	$(SRCROOT)mom5/ocean_param/neutral/ocean_nphysics_skew.F90
ocean_nphysics_tensor.o: $(SRCROOT)mom5/ocean_param/neutral/ocean_nphysics_tensor.F90 /home/wenchang/tigress/CM2.5_tiger2_wy_killtc/exp/../src/mom5/ocean_core/ocean_memory.h constants.o fms.o mpp.o ocean_domains.o ocean_nphysics_skew.o ocean_sigma_transport.o ocean_types.o ocean_util.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -I$(SRCROOT)mom5/ocean_core	$(SRCROOT)mom5/ocean_param/neutral/ocean_nphysics_tensor.F90
ocean_nphysics_util.o: $(SRCROOT)mom5/ocean_param/neutral/ocean_nphysics_util.F90 /home/wenchang/tigress/CM2.5_tiger2_wy_killtc/exp/../src/mom5/ocean_core/ocean_memory.h constants.o diag_manager.o fms.o fms_io.o mpp.o mpp_domains.o time_manager.o ocean_domains.o ocean_density.o ocean_operators.o ocean_parameters.o ocean_tracer_diag.o ocean_types.o ocean_util.o ocean_tracer_util.o ocean_workspace.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -I$(SRCROOT)mom5/ocean_core	$(SRCROOT)mom5/ocean_param/neutral/ocean_nphysics_util.F90
ocean_nphysics_util_new.o: $(SRCROOT)mom5/ocean_param/neutral/ocean_nphysics_util_new.F90 /home/wenchang/tigress/CM2.5_tiger2_wy_killtc/exp/../src/mom5/ocean_core/ocean_memory.h constants.o diag_manager.o fms.o mpp.o mpp_domains.o ocean_domains.o ocean_parameters.o ocean_types.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -I$(SRCROOT)mom5/ocean_core	$(SRCROOT)mom5/ocean_param/neutral/ocean_nphysics_util_new.F90
ocean_obc.o: $(SRCROOT)mom5/ocean_core/ocean_obc.F90 /home/wenchang/tigress/CM2.5_tiger2_wy_killtc/exp/../src/shared/include/fms_platform.h constants.o data_override.o diag_manager.o fms.o fms_io.o mpp_io.o mpp_domains.o mpp.o time_interp_external.o time_manager.o tracer_manager.o ocean_util.o ocean_domains.o ocean_parameters.o ocean_types.o ocean_obc_barotrop.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -I$(SRCROOT)shared/include	$(SRCROOT)mom5/ocean_core/ocean_obc.F90
ocean_obc_barotrop.o: $(SRCROOT)mom5/ocean_core/ocean_obc_barotrop.F90 /home/wenchang/tigress/CM2.5_tiger2_wy_killtc/exp/../src/shared/include/fms_platform.h constants.o data_override.o diag_manager.o fms.o fms_io.o mpp_io.o mpp_domains.o mpp.o time_interp_external.o time_manager.o tracer_manager.o ocean_util.o ocean_domains.o ocean_parameters.o ocean_types.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -I$(SRCROOT)shared/include	$(SRCROOT)mom5/ocean_core/ocean_obc_barotrop.F90
ocean_operators.o: $(SRCROOT)mom5/ocean_core/ocean_operators.F90 /home/wenchang/tigress/CM2.5_tiger2_wy_killtc/exp/../src/mom5/ocean_core/ocean_memory.h mpp_domains.o mpp.o fms.o ocean_domains.o ocean_parameters.o ocean_types.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -I$(SRCROOT)mom5/ocean_core	$(SRCROOT)mom5/ocean_core/ocean_operators.F90
ocean_overexchange.o: $(SRCROOT)mom5/ocean_param/sources/ocean_overexchange.F90 /home/wenchang/tigress/CM2.5_tiger2_wy_killtc/exp/../src/mom5/ocean_core/ocean_memory.h constants.o diag_manager.o fms.o mpp_domains.o mpp.o ocean_density.o ocean_domains.o ocean_parameters.o ocean_types.o ocean_util.o ocean_tracer_util.o ocean_workspace.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -I$(SRCROOT)mom5/ocean_core	$(SRCROOT)mom5/ocean_param/sources/ocean_overexchange.F90
ocean_overflow.o: $(SRCROOT)mom5/ocean_param/sources/ocean_overflow.F90 /home/wenchang/tigress/CM2.5_tiger2_wy_killtc/exp/../src/mom5/ocean_core/ocean_memory.h constants.o diag_manager.o fms.o mpp_domains.o mpp.o ocean_density.o ocean_domains.o ocean_parameters.o ocean_types.o ocean_util.o ocean_tracer_util.o ocean_workspace.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -I$(SRCROOT)mom5/ocean_core	$(SRCROOT)mom5/ocean_param/sources/ocean_overflow.F90
ocean_overflow_OFP.o: $(SRCROOT)mom5/ocean_param/sources/ocean_overflow_OFP.F90 /home/wenchang/tigress/CM2.5_tiger2_wy_killtc/exp/../src/mom5/ocean_core/ocean_memory.h constants.o diag_manager.o field_manager.o fms.o mpp_domains.o mpp.o ocean_density.o ocean_domains.o ocean_parameters.o ocean_types.o ocean_util.o ocean_tracer_util.o ocean_workspace.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -I$(SRCROOT)mom5/ocean_core	$(SRCROOT)mom5/ocean_param/sources/ocean_overflow_OFP.F90
ocean_parameters.o: $(SRCROOT)mom5/ocean_core/ocean_parameters.F90 mpp.o fms.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)mom5/ocean_core/ocean_parameters.F90
ocean_passive.o: $(SRCROOT)mom5/ocean_tracers/ocean_passive.F90 /home/wenchang/tigress/CM2.5_tiger2_wy_killtc/exp/../src/mom5/ocean_core/ocean_memory.h constants.o field_manager.o fm_util.o fms.o mpp.o ocean_domains.o ocean_tpm_util.o ocean_types.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -I$(SRCROOT)mom5/ocean_core	$(SRCROOT)mom5/ocean_tracers/ocean_passive.F90
ocean_pert_co2.o: $(SRCROOT)mom5/ocean_bgc/ocean_pert_co2.F90 /home/wenchang/tigress/CM2.5_tiger2_wy_killtc/exp/../src/shared/include/fms_platform.h atmos_ocean_fluxes.o diag_manager.o field_manager.o time_manager.o fms.o mpp.o mpp_domains.o constants.o ocean_tpm_util.o fm_util.o coupler_types.o ocean_types.o ocmip2_co2calc.o ocean_util.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -I$(SRCROOT)shared/include	$(SRCROOT)mom5/ocean_bgc/ocean_pert_co2.F90
ocean_po4_pre.o: $(SRCROOT)mom5/ocean_bgc/ocean_po4_pre.F90 /home/wenchang/tigress/CM2.5_tiger2_wy_killtc/exp/../src/shared/include/fms_platform.h diag_manager.o field_manager.o mpp.o fms.o time_manager.o mpp_domains.o ocean_tpm_util.o fm_util.o ocean_types.o ocean_tracer_diag.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -I$(SRCROOT)shared/include	$(SRCROOT)mom5/ocean_bgc/ocean_po4_pre.F90
ocean_pressure.o: $(SRCROOT)mom5/ocean_core/ocean_pressure.F90 /home/wenchang/tigress/CM2.5_tiger2_wy_killtc/exp/../src/mom5/ocean_core/ocean_memory.h constants.o diag_manager.o fms.o mpp_io.o mpp.o mpp_domains.o ocean_domains.o ocean_operators.o ocean_parameters.o ocean_types.o ocean_util.o ocean_workspace.o ocean_obc.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -I$(SRCROOT)mom5/ocean_core	$(SRCROOT)mom5/ocean_core/ocean_pressure.F90
ocean_residency.o: $(SRCROOT)mom5/ocean_tracers/ocean_residency.F90 time_manager.o field_manager.o mpp.o diag_manager.o ocean_tpm_util.o fm_util.o ocean_types.o ocean_residency_meta.o ocean_residency_ml.o ocean_residency_range.o ocean_residency_integrand.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)mom5/ocean_tracers/ocean_residency.F90
ocean_residency_integrand.o: $(SRCROOT)mom5/ocean_tracers/ocean_residency_integrand.F90 field_manager.o mpp.o fms.o ocean_residency_meta.o ocean_types.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)mom5/ocean_tracers/ocean_residency_integrand.F90
ocean_residency_meta.o: $(SRCROOT)mom5/ocean_tracers/ocean_residency_meta.F90 field_manager.o fm_util.o mpp.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)mom5/ocean_tracers/ocean_residency_meta.F90
ocean_residency_ml.o: $(SRCROOT)mom5/ocean_tracers/ocean_residency_ml.F90 field_manager.o mpp.o ocean_residency_meta.o ocean_types.o ocean_tracer_diag.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)mom5/ocean_tracers/ocean_residency_ml.F90
ocean_residency_range.o: $(SRCROOT)mom5/ocean_tracers/ocean_residency_range.F90 /home/wenchang/tigress/CM2.5_tiger2_wy_killtc/exp/../src/shared/include/fms_platform.h field_manager.o mpp.o ocean_residency_meta.o ocean_types.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -I$(SRCROOT)shared/include	$(SRCROOT)mom5/ocean_tracers/ocean_residency_range.F90
ocean_rivermix.o: $(SRCROOT)mom5/ocean_param/sources/ocean_rivermix.F90 axis_utils.o constants.o diag_manager.o fms.o mpp_domains.o mpp.o ocean_domains.o ocean_parameters.o ocean_types.o ocean_workspace.o ocean_util.o ocean_tracer_util.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)mom5/ocean_param/sources/ocean_rivermix.F90
ocean_riverspread.o: $(SRCROOT)mom5/ocean_param/sources/ocean_riverspread.F90 field_manager.o fms.o mpp_domains.o mpp.o ocean_domains.o ocean_operators.o ocean_types.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)mom5/ocean_param/sources/ocean_riverspread.F90
ocean_rough.o: $(SRCROOT)ice_param/ocean_rough.F90 mpp.o fms.o constants.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)ice_param/ocean_rough.F90
ocean_sbc.o: $(SRCROOT)mom5/ocean_core/ocean_sbc.F90 /home/wenchang/tigress/CM2.5_tiger2_wy_killtc/exp/../src/shared/include/fms_platform.h /home/wenchang/tigress/CM2.5_tiger2_wy_killtc/exp/../src/mom5/ocean_core/ocean_memory.h constants.o diag_manager.o fms.o fms_io.o mpp_domains.o mpp.o time_interp_external.o time_manager.o ocean_domains.o ocean_parameters.o ocean_riverspread.o ocean_tpm.o ocean_types.o ocean_workspace.o ocean_util.o ocean_tracer_util.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -I$(SRCROOT)shared/include -I$(SRCROOT)mom5/ocean_core	$(SRCROOT)mom5/ocean_core/ocean_sbc.F90
ocean_shortwave.o: $(SRCROOT)mom5/ocean_param/sources/ocean_shortwave.F90 /home/wenchang/tigress/CM2.5_tiger2_wy_killtc/exp/../src/shared/include/fms_platform.h /home/wenchang/tigress/CM2.5_tiger2_wy_killtc/exp/../src/mom5/ocean_core/ocean_memory.h constants.o diag_manager.o fms.o mpp.o ocean_domains.o ocean_parameters.o ocean_shortwave_csiro.o ocean_shortwave_jerlov.o ocean_shortwave_gfdl.o ocean_types.o ocean_tpm_util.o ocean_workspace.o ocean_util.o ocean_tracer_util.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -I$(SRCROOT)shared/include -I$(SRCROOT)mom5/ocean_core	$(SRCROOT)mom5/ocean_param/sources/ocean_shortwave.F90
ocean_shortwave_csiro.o: $(SRCROOT)mom5/ocean_param/sources/ocean_shortwave_csiro.F90 /home/wenchang/tigress/CM2.5_tiger2_wy_killtc/exp/../src/mom5/ocean_core/ocean_memory.h axis_utils.o field_manager.o fms.o mpp.o time_interp_external.o ocean_domains.o ocean_types.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -I$(SRCROOT)mom5/ocean_core	$(SRCROOT)mom5/ocean_param/sources/ocean_shortwave_csiro.F90
ocean_shortwave_gfdl.o: $(SRCROOT)mom5/ocean_param/sources/ocean_shortwave_gfdl.F90 /home/wenchang/tigress/CM2.5_tiger2_wy_killtc/exp/../src/mom5/ocean_core/ocean_memory.h axis_utils.o constants.o diag_manager.o field_manager.o fms.o mpp.o time_interp_external.o ocean_domains.o ocean_parameters.o ocean_types.o ocean_workspace.o ocean_util.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -I$(SRCROOT)mom5/ocean_core	$(SRCROOT)mom5/ocean_param/sources/ocean_shortwave_gfdl.F90
ocean_shortwave_jerlov.o: $(SRCROOT)mom5/ocean_param/sources/ocean_shortwave_jerlov.F90 /home/wenchang/tigress/CM2.5_tiger2_wy_killtc/exp/../src/mom5/ocean_core/ocean_memory.h constants.o fms.o mpp.o time_interp_external.o ocean_domains.o ocean_parameters.o ocean_types.o ocean_workspace.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -I$(SRCROOT)mom5/ocean_core	$(SRCROOT)mom5/ocean_param/sources/ocean_shortwave_jerlov.F90
ocean_sigma_transport.o: $(SRCROOT)mom5/ocean_param/lateral/ocean_sigma_transport.F90 /home/wenchang/tigress/CM2.5_tiger2_wy_killtc/exp/../src/mom5/ocean_core/ocean_memory.h constants.o diag_manager.o fms.o fms_io.o mpp_domains.o mpp.o ocean_domains.o ocean_density.o ocean_operators.o ocean_parameters.o ocean_types.o ocean_util.o ocean_tracer_util.o ocean_workspace.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -I$(SRCROOT)mom5/ocean_core	$(SRCROOT)mom5/ocean_param/lateral/ocean_sigma_transport.F90
ocean_sponges_eta.o: $(SRCROOT)mom5/ocean_param/sources/ocean_sponges_eta.F90 /home/wenchang/tigress/CM2.5_tiger2_wy_killtc/exp/../src/mom5/ocean_core/ocean_memory.h diag_manager.o fms.o mpp.o time_interp_external.o time_manager.o ocean_domains.o ocean_parameters.o ocean_types.o ocean_workspace.o ocean_util.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -I$(SRCROOT)mom5/ocean_core	$(SRCROOT)mom5/ocean_param/sources/ocean_sponges_eta.F90
ocean_sponges_tracer.o: $(SRCROOT)mom5/ocean_param/sources/ocean_sponges_tracer.F90 /home/wenchang/tigress/CM2.5_tiger2_wy_killtc/exp/../src/mom5/ocean_core/ocean_memory.h diag_manager.o fms.o mpp.o time_interp_external.o time_manager.o ocean_domains.o ocean_parameters.o ocean_types.o ocean_workspace.o ocean_util.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -I$(SRCROOT)mom5/ocean_core	$(SRCROOT)mom5/ocean_param/sources/ocean_sponges_tracer.F90
ocean_sponges_velocity.o: $(SRCROOT)mom5/ocean_param/sources/ocean_sponges_velocity.F90 /home/wenchang/tigress/CM2.5_tiger2_wy_killtc/exp/../src/mom5/ocean_core/ocean_memory.h diag_manager.o fms.o mpp.o time_interp_external.o time_manager.o ocean_domains.o ocean_parameters.o ocean_types.o ocean_workspace.o ocean_util.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -I$(SRCROOT)mom5/ocean_core	$(SRCROOT)mom5/ocean_param/sources/ocean_sponges_velocity.F90
ocean_submesoscale.o: $(SRCROOT)mom5/ocean_param/lateral/ocean_submesoscale.F90 /home/wenchang/tigress/CM2.5_tiger2_wy_killtc/exp/../src/mom5/ocean_core/ocean_memory.h constants.o diag_manager.o fms.o mpp_domains.o mpp.o time_manager.o ocean_domains.o ocean_operators.o ocean_parameters.o ocean_tracer_diag.o ocean_types.o ocean_util.o ocean_tracer_util.o ocean_workspace.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -I$(SRCROOT)mom5/ocean_core	$(SRCROOT)mom5/ocean_param/lateral/ocean_submesoscale.F90
ocean_tempsalt.o: $(SRCROOT)mom5/ocean_tracers/ocean_tempsalt.F90 /home/wenchang/tigress/CM2.5_tiger2_wy_killtc/exp/../src/mom5/ocean_core/ocean_memory.h constants.o fms.o mpp.o ocean_domains.o ocean_parameters.o ocean_tpm_util.o ocean_types.o ocean_workspace.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -I$(SRCROOT)mom5/ocean_core	$(SRCROOT)mom5/ocean_tracers/ocean_tempsalt.F90
ocean_thickness.o: $(SRCROOT)mom5/ocean_core/ocean_thickness.F90 /home/wenchang/tigress/CM2.5_tiger2_wy_killtc/exp/../src/mom5/ocean_core/ocean_memory.h constants.o diag_manager.o fms.o fms_io.o mpp_domains.o mpp.o ocean_domains.o ocean_grids.o ocean_parameters.o ocean_tracer_util.o ocean_types.o ocean_util.o ocean_workspace.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -I$(SRCROOT)mom5/ocean_core	$(SRCROOT)mom5/ocean_core/ocean_thickness.F90
ocean_topog.o: $(SRCROOT)mom5/ocean_core/ocean_topog.F90 /home/wenchang/tigress/CM2.5_tiger2_wy_killtc/exp/../src/mom5/ocean_core/ocean_memory.h fms.o mpp_domains.o mpp.o axis_utils.o ocean_domains.o ocean_parameters.o ocean_types.o ocean_util.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -I$(SRCROOT)mom5/ocean_core	$(SRCROOT)mom5/ocean_core/ocean_topog.F90
ocean_tpm.o: $(SRCROOT)mom5/ocean_tracers/ocean_tpm.F90 mpp.o mpp_domains.o ocean_types.o ocean_age_tracer.o ocean_residency.o ocean_pert_co2.o ocmip2_abiotic.o ocmip2_cfc.o ocmip2_biotic.o ocean_bgc_restore.o ocmip2_he.o ocean_po4_pre.o ocean_ibgc.o ocean_generic_tracer.o ocean_frazil.o ocean_tempsalt.o ocean_passive.o transport_matrix.o time_manager.o coupler_types.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)mom5/ocean_tracers/ocean_tpm.F90
ocean_tpm_util.o: $(SRCROOT)mom5/ocean_tracers/ocean_tpm_util.F90 field_manager.o fms.o mpp.o fm_util.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)mom5/ocean_tracers/ocean_tpm_util.F90
ocean_tracer.o: $(SRCROOT)mom5/ocean_tracers/ocean_tracer.F90 /home/wenchang/tigress/CM2.5_tiger2_wy_killtc/exp/../src/mom5/ocean_core/ocean_memory.h constants.o diag_manager.o field_manager.o fms.o fms_io.o fm_util.o mpp_domains.o mpp_io.o mpp.o platform.o time_manager.o transport_matrix.o ocean_blob.o ocean_convect.o ocean_domains.o ocean_density.o ocean_frazil.o ocean_bih_tracer.o ocean_generic_tracer.o ocean_lap_tracer.o ocean_obc.o ocean_parameters.o ocean_passive.o ocean_shortwave.o ocean_tempsalt.o ocean_tpm.o ocean_tpm_util.o ocean_tracer_advect.o ocean_tracer_diag.o ocean_tracer_util.o ocean_thickness.o ocean_types.o ocean_util.o ocean_vert_mix.o ocean_workspace.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -I$(SRCROOT)mom5/ocean_core	$(SRCROOT)mom5/ocean_tracers/ocean_tracer.F90
ocean_tracer_advect.o: $(SRCROOT)mom5/ocean_tracers/ocean_tracer_advect.F90 /home/wenchang/tigress/CM2.5_tiger2_wy_killtc/exp/../src/mom5/ocean_core/ocean_memory.h constants.o diag_manager.o fms.o fms_io.o mpp_domains.o mpp.o ocean_domains.o ocean_obc.o ocean_parameters.o ocean_topog.o ocean_tracer_util.o ocean_types.o ocean_workspace.o ocean_util.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -I$(SRCROOT)mom5/ocean_core	$(SRCROOT)mom5/ocean_tracers/ocean_tracer_advect.F90
ocean_tracer_diag.o: $(SRCROOT)mom5/ocean_diag/ocean_tracer_diag.F90 /home/wenchang/tigress/CM2.5_tiger2_wy_killtc/exp/../src/mom5/ocean_core/ocean_memory.h constants.o diag_manager.o fms.o mpp_domains.o mpp.o time_manager.o ocean_density.o ocean_domains.o ocean_obc.o ocean_operators.o ocean_parameters.o ocean_tracer_util.o ocean_types.o ocean_util.o ocean_workspace.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -I$(SRCROOT)mom5/ocean_core	$(SRCROOT)mom5/ocean_diag/ocean_tracer_diag.F90
ocean_tracer_util.o: $(SRCROOT)mom5/ocean_diag/ocean_tracer_util.F90 /home/wenchang/tigress/CM2.5_tiger2_wy_killtc/exp/../src/mom5/ocean_core/ocean_memory.h constants.o fms.o mpp.o platform.o diag_manager.o ocean_domains.o ocean_parameters.o ocean_types.o ocean_util.o ocean_workspace.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -I$(SRCROOT)mom5/ocean_core	$(SRCROOT)mom5/ocean_diag/ocean_tracer_util.F90
ocean_types.o: $(SRCROOT)mom5/ocean_core/ocean_types.F90 /home/wenchang/tigress/CM2.5_tiger2_wy_killtc/exp/../src/shared/include/fms_platform.h /home/wenchang/tigress/CM2.5_tiger2_wy_killtc/exp/../src/mom5/ocean_core/ocean_memory.h coupler_types.o field_manager.o fms.o mpp_domains.o mpp.o time_manager.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -I$(SRCROOT)shared/include -I$(SRCROOT)mom5/ocean_core	$(SRCROOT)mom5/ocean_core/ocean_types.F90
ocean_util.o: $(SRCROOT)mom5/ocean_core/ocean_util.F90 /home/wenchang/tigress/CM2.5_tiger2_wy_killtc/exp/../src/mom5/ocean_core/ocean_memory.h platform.o constants.o diag_manager.o mpp_domains.o mpp.o time_manager.o ocean_domains.o ocean_parameters.o ocean_types.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -I$(SRCROOT)mom5/ocean_core	$(SRCROOT)mom5/ocean_core/ocean_util.F90
ocean_velocity.o: $(SRCROOT)mom5/ocean_core/ocean_velocity.F90 /home/wenchang/tigress/CM2.5_tiger2_wy_killtc/exp/../src/mom5/ocean_core/ocean_memory.h constants.o diag_manager.o fms.o fms_io.o mpp_domains.o mpp.o time_interp_external.o ocean_bih_friction.o ocean_coriolis.o ocean_domains.o ocean_form_drag.o ocean_lap_friction.o ocean_momentum_source.o ocean_obc.o ocean_operators.o ocean_parameters.o ocean_pressure.o ocean_types.o ocean_util.o ocean_velocity_advect.o ocean_velocity_diag.o ocean_vert_mix.o ocean_workspace.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -I$(SRCROOT)mom5/ocean_core	$(SRCROOT)mom5/ocean_core/ocean_velocity.F90
ocean_velocity_advect.o: $(SRCROOT)mom5/ocean_core/ocean_velocity_advect.F90 /home/wenchang/tigress/CM2.5_tiger2_wy_killtc/exp/../src/mom5/ocean_core/ocean_memory.h diag_manager.o fms.o mpp.o ocean_domains.o ocean_obc.o ocean_operators.o ocean_parameters.o ocean_types.o ocean_workspace.o ocean_util.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -I$(SRCROOT)mom5/ocean_core	$(SRCROOT)mom5/ocean_core/ocean_velocity_advect.F90
ocean_velocity_diag.o: $(SRCROOT)mom5/ocean_diag/ocean_velocity_diag.F90 /home/wenchang/tigress/CM2.5_tiger2_wy_killtc/exp/../src/shared/include/fms_platform.h /home/wenchang/tigress/CM2.5_tiger2_wy_killtc/exp/../src/mom5/ocean_core/ocean_memory.h constants.o diag_manager.o fms.o mpp_domains.o mpp.o time_manager.o ocean_bih_friction.o ocean_coriolis.o ocean_domains.o ocean_form_drag.o ocean_lap_friction.o ocean_momentum_source.o ocean_operators.o ocean_parameters.o ocean_pressure.o ocean_types.o ocean_util.o ocean_velocity_advect.o ocean_vert_mix.o ocean_workspace.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -I$(SRCROOT)shared/include -I$(SRCROOT)mom5/ocean_core	$(SRCROOT)mom5/ocean_diag/ocean_velocity_diag.F90
ocean_vert_chen.o: $(SRCROOT)mom5/ocean_param/vertical/ocean_vert_chen.F90 constants.o diag_manager.o fms.o fms_io.o mpp_domains.o mpp.o ocean_density.o ocean_domains.o ocean_parameters.o ocean_shortwave_csiro.o ocean_types.o ocean_util.o ocean_vert_util.o ocean_workspace.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)mom5/ocean_param/vertical/ocean_vert_chen.F90
ocean_vert_const.o: $(SRCROOT)mom5/ocean_param/vertical/ocean_vert_const.F90 constants.o diag_manager.o fms.o mpp_io.o mpp.o ocean_density.o ocean_domains.o ocean_parameters.o ocean_types.o ocean_workspace.o ocean_util.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)mom5/ocean_param/vertical/ocean_vert_const.F90
ocean_vert_gotm.o: $(SRCROOT)mom5/ocean_param/vertical/ocean_vert_gotm.F90 /home/wenchang/tigress/CM2.5_tiger2_wy_killtc/exp/../src/mom5/ocean_core/ocean_memory.h constants.o diag_manager.o fms.o fms_io.o mpp_domains.o mpp.o tridiagonal_gotm.o turbulence.o ocean_domains.o ocean_parameters.o ocean_types.o ocean_util.o ocean_workspace.o ocean_obc.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -I$(SRCROOT)mom5/ocean_core	$(SRCROOT)mom5/ocean_param/vertical/ocean_vert_gotm.F90
ocean_vert_kpp_mom4p0.o: $(SRCROOT)mom5/ocean_param/vertical/ocean_vert_kpp_mom4p0.F90 /home/wenchang/tigress/CM2.5_tiger2_wy_killtc/exp/../src/mom5/ocean_core/ocean_memory.h constants.o diag_manager.o fms.o mpp_domains.o mpp.o ocean_density.o ocean_domains.o ocean_parameters.o ocean_types.o ocean_workspace.o ocean_util.o ocean_tracer_util.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -I$(SRCROOT)mom5/ocean_core	$(SRCROOT)mom5/ocean_param/vertical/ocean_vert_kpp_mom4p0.F90
ocean_vert_kpp_mom4p1.o: $(SRCROOT)mom5/ocean_param/vertical/ocean_vert_kpp_mom4p1.F90 /home/wenchang/tigress/CM2.5_tiger2_wy_killtc/exp/../src/mom5/ocean_core/ocean_memory.h constants.o diag_manager.o fms.o mpp_domains.o mpp.o ocean_density.o ocean_domains.o ocean_parameters.o ocean_types.o ocean_workspace.o ocean_util.o ocean_tracer_util.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -I$(SRCROOT)mom5/ocean_core	$(SRCROOT)mom5/ocean_param/vertical/ocean_vert_kpp_mom4p1.F90
ocean_vert_kpp_test.o: $(SRCROOT)mom5/ocean_param/vertical/ocean_vert_kpp_test.F90 /home/wenchang/tigress/CM2.5_tiger2_wy_killtc/exp/../src/mom5/ocean_core/ocean_memory.h constants.o diag_manager.o fms.o mpp_domains.o mpp.o ocean_density.o ocean_domains.o ocean_parameters.o ocean_types.o ocean_vert_util.o ocean_workspace.o ocean_util.o ocean_tracer_util.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -I$(SRCROOT)mom5/ocean_core	$(SRCROOT)mom5/ocean_param/vertical/ocean_vert_kpp_test.F90
ocean_vert_mix.o: $(SRCROOT)mom5/ocean_param/vertical/ocean_vert_mix.F90 constants.o diag_manager.o field_manager.o fms.o mpp_domains.o mpp.o ocean_domains.o ocean_parameters.o ocean_tracer_util.o ocean_types.o ocean_util.o ocean_vert_const.o ocean_vert_chen.o ocean_vert_gotm.o ocean_vert_kpp_test.o ocean_vert_kpp_mom4p0.o ocean_vert_kpp_mom4p1.o ocean_vert_pp.o ocean_vert_tidal.o ocean_vert_util.o ocean_workspace.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)mom5/ocean_param/vertical/ocean_vert_mix.F90
ocean_vert_pp.o: $(SRCROOT)mom5/ocean_param/vertical/ocean_vert_pp.F90 constants.o diag_manager.o fms.o mpp.o ocean_density.o ocean_domains.o ocean_parameters.o ocean_types.o ocean_vert_util.o ocean_workspace.o ocean_util.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)mom5/ocean_param/vertical/ocean_vert_pp.F90
ocean_vert_tidal.o: $(SRCROOT)mom5/ocean_param/vertical/ocean_vert_tidal.F90 /home/wenchang/tigress/CM2.5_tiger2_wy_killtc/exp/../src/mom5/ocean_core/ocean_memory.h constants.o diag_manager.o fms.o mpp_domains.o mpp.o ocean_domains.o ocean_operators.o ocean_parameters.o ocean_types.o ocean_workspace.o ocean_util.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -I$(SRCROOT)mom5/ocean_core	$(SRCROOT)mom5/ocean_param/vertical/ocean_vert_tidal.F90
ocean_vert_tidal_test.o: $(SRCROOT)mom5/ocean_param/vertical/ocean_vert_tidal_test.F90 /home/wenchang/tigress/CM2.5_tiger2_wy_killtc/exp/../src/mom5/ocean_core/ocean_memory.h constants.o diag_manager.o fms.o mpp_domains.o mpp.o ocean_domains.o ocean_operators.o ocean_parameters.o ocean_types.o ocean_workspace.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -I$(SRCROOT)mom5/ocean_core	$(SRCROOT)mom5/ocean_param/vertical/ocean_vert_tidal_test.F90
ocean_vert_util.o: $(SRCROOT)mom5/ocean_param/vertical/ocean_vert_util.F90 /home/wenchang/tigress/CM2.5_tiger2_wy_killtc/exp/../src/mom5/ocean_core/ocean_memory.h constants.o diag_manager.o fms.o mpp.o ocean_density.o ocean_domains.o ocean_parameters.o ocean_types.o ocean_workspace.o ocean_util.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -I$(SRCROOT)mom5/ocean_core	$(SRCROOT)mom5/ocean_param/vertical/ocean_vert_util.F90
ocean_wave.o: $(SRCROOT)mom5/ocean_wave/ocean_wave.F90 /home/wenchang/tigress/CM2.5_tiger2_wy_killtc/exp/../src/mom5/ocean_core/ocean_memory.h constants.o diag_manager.o fms.o mpp.o fms_io.o ocean_domains.o mpp_domains.o ocean_types.o ocean_parameters.o ocean_workspace.o ocean_util.o ocean_wave_type.o data_override.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -I$(SRCROOT)mom5/ocean_core	$(SRCROOT)mom5/ocean_wave/ocean_wave.F90
ocean_wave_type.o: $(SRCROOT)mom5/ocean_wave/ocean_wave_type.F90 /home/wenchang/tigress/CM2.5_tiger2_wy_killtc/exp/../src/shared/include/fms_platform.h /home/wenchang/tigress/CM2.5_tiger2_wy_killtc/exp/../src/mom5/ocean_core/ocean_memory.h
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -I$(SRCROOT)shared/include -I$(SRCROOT)mom5/ocean_core	$(SRCROOT)mom5/ocean_wave/ocean_wave_type.F90
ocean_workspace.o: $(SRCROOT)mom5/ocean_core/ocean_workspace.F90 /home/wenchang/tigress/CM2.5_tiger2_wy_killtc/exp/../src/mom5/ocean_core/ocean_memory.h ocean_types.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -I$(SRCROOT)mom5/ocean_core	$(SRCROOT)mom5/ocean_core/ocean_workspace.F90
ocean_xlandinsert.o: $(SRCROOT)mom5/ocean_param/sources/ocean_xlandinsert.F90 constants.o diag_manager.o field_manager.o fms.o mpp_domains.o mpp.o ocean_domains.o ocean_parameters.o ocean_types.o ocean_workspace.o ocean_util.o ocean_tracer_util.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)mom5/ocean_param/sources/ocean_xlandinsert.F90
ocean_xlandmix.o: $(SRCROOT)mom5/ocean_param/sources/ocean_xlandmix.F90 constants.o diag_manager.o field_manager.o fms.o mpp_domains.o mpp.o ocean_domains.o ocean_parameters.o ocean_types.o ocean_workspace.o ocean_util.o ocean_tracer_util.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)mom5/ocean_param/sources/ocean_xlandmix.F90
ocmip2_abiotic.o: $(SRCROOT)mom5/ocean_bgc/ocmip2_abiotic.F90 /home/wenchang/tigress/CM2.5_tiger2_wy_killtc/exp/../src/shared/include/fms_platform.h atmos_ocean_fluxes.o diag_manager.o field_manager.o time_manager.o fms.o fms_io.o mpp.o time_interp_external.o mpp_domains.o constants.o ocean_tpm_util.o fm_util.o coupler_types.o ocean_types.o ocmip2_co2calc.o ocean_util.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -I$(SRCROOT)shared/include	$(SRCROOT)mom5/ocean_bgc/ocmip2_abiotic.F90
ocmip2_biotic.o: $(SRCROOT)mom5/ocean_bgc/ocmip2_biotic.F90 /home/wenchang/tigress/CM2.5_tiger2_wy_killtc/exp/../src/shared/include/fms_platform.h time_manager.o diag_manager.o field_manager.o fms.o fms_io.o mpp.o time_interp_external.o mpp_domains.o constants.o ocean_tpm_util.o fm_util.o coupler_types.o ocean_types.o ocmip2_co2calc.o ocean_util.o atmos_ocean_fluxes.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -I$(SRCROOT)shared/include	$(SRCROOT)mom5/ocean_bgc/ocmip2_biotic.F90
ocmip2_cfc.o: $(SRCROOT)mom5/ocean_bgc/ocmip2_cfc.F90 /home/wenchang/tigress/CM2.5_tiger2_wy_killtc/exp/../src/shared/include/fms_platform.h atmos_ocean_fluxes.o time_manager.o field_manager.o fms.o ocean_tpm_util.o constants.o fm_util.o mpp.o coupler_types.o ocean_types.o diag_manager.o ocean_util.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -I$(SRCROOT)shared/include	$(SRCROOT)mom5/ocean_bgc/ocmip2_cfc.F90
ocmip2_co2calc.o: $(SRCROOT)mom5/ocean_bgc/ocmip2_co2calc.F90
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)mom5/ocean_bgc/ocmip2_co2calc.F90
ocmip2_he.o: $(SRCROOT)mom5/ocean_bgc/ocmip2_he.F90 /home/wenchang/tigress/CM2.5_tiger2_wy_killtc/exp/../src/shared/include/fms_platform.h atmos_ocean_fluxes.o time_manager.o time_interp_external.o diag_manager.o field_manager.o fms.o fms_io.o ocean_tpm_util.o fm_util.o mpp.o coupler_types.o ocean_types.o mpp_domains.o ocean_util.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -I$(SRCROOT)shared/include	$(SRCROOT)mom5/ocean_bgc/ocmip2_he.F90
oda_core.o: $(SRCROOT)shared/oda_tools/oda_core.F90 fms.o mpp.o mpp_domains.o time_manager.o get_cal_time.o axis_utils.o constants.o oda_types.o write_ocean_data.o mpp_io.o field_manager.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)shared/oda_tools/oda_core.F90
oda_core_ecda.o: $(SRCROOT)shared/oda_tools/oda_core_ecda.F90 fms.o mpp.o mpp_io.o mpp_domains.o mpp_memutils.o time_manager.o get_cal_time.o axis_utils.o horiz_interp_type.o horiz_interp_bilinear.o constants.o oda_types.o xbt_drop_rate_adjust.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)shared/oda_tools/oda_core_ecda.F90
oda_driver.o: $(SRCROOT)mom5/ocean_core/oda_driver.F90 constants.o ensemble_manager.o mpp_domains.o mpp.o mpp_io.o oda_types.o oda_core.o time_manager.o write_ocean_data.o ocean_convect.o ocean_domains.o ocean_grids.o ocean_parameters.o ocean_topog.o ocean_types.o ocean_workspace.o fms.o diag_manager.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)mom5/ocean_core/oda_driver.F90
oda_types.o: $(SRCROOT)shared/oda_tools/oda_types.F90 time_manager.o mpp.o mpp_domains.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)shared/oda_tools/oda_types.F90
ode_solvers.o: $(SRCROOT)mom5/ocean_param/gotm-4.0/util/ode_solvers.F90 /home/wenchang/tigress/CM2.5_tiger2_wy_killtc/exp/../src/mom5/ocean_param/gotm-4.0/include/cppdefs.h /home/wenchang/tigress/CM2.5_tiger2_wy_killtc/exp/../src/mom5/ocean_param/gotm-4.0/include/version.h
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -I$(SRCROOT)mom5/ocean_param/gotm-4.0/include	$(SRCROOT)mom5/ocean_param/gotm-4.0/util/ode_solvers.F90
optical_path.o: $(SRCROOT)atmos_param/sea_esf_rad/optical_path.F90 mpp.o fms.o constants.o rad_utilities.o longwave_params.o lw_gases_stdtf.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)atmos_param/sea_esf_rad/optical_path.F90
optics_lib.o: $(SRCROOT)atmos_param/cosp/quickbeam/optics_lib.F90 /home/wenchang/tigress/CM2.5_tiger2_wy_killtc/exp/../src/atmos_param/cosp/cosp_defs.H
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -I$(SRCROOT)atmos_param/cosp	$(SRCROOT)atmos_param/cosp/quickbeam/optics_lib.F90
original_fms_rad.o: $(SRCROOT)atmos_param/sea_esf_rad/original_fms_rad.F90 mcm_lw.o mcm_sw_driver.o fsrad.o clouds.o time_manager.o mpp.o fms.o rad_utilities.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)atmos_param/sea_esf_rad/original_fms_rad.F90
output_grads.o: $(SRCROOT)atmos_param/clubb/CLUBB_core/output_grads.F90 /home/wenchang/tigress/CM2.5_tiger2_wy_killtc/exp/../src/atmos_param/clubb/CLUBB_core/recl.inc clubb_precision.o constants_clubb.o stat_file_module.o model_flags.o endian.o calendar.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -I$(SRCROOT)atmos_param/clubb/CLUBB_core	$(SRCROOT)atmos_param/clubb/CLUBB_core/output_grads.F90
output_netcdf.o: $(SRCROOT)atmos_param/clubb/CLUBB_core/output_netcdf.F90 stat_file_module.o clubb_precision.o constants_clubb.o parameters_model.o parameters_tunable.o parameter_indices.o model_flags.o parameters_microphys.o parameters_radiation.o calendar.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)atmos_param/clubb/CLUBB_core/output_netcdf.F90
ozone.o: $(SRCROOT)atmos_param/sea_esf_rad/ozone.F90 mpp.o fms.o fms_io.o time_manager.o diag_manager.o time_interp.o constants.o interpolator.o tracer_manager.o field_manager.o rad_utilities.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)atmos_param/sea_esf_rad/ozone.F90
parameter_indices.o: $(SRCROOT)atmos_param/clubb/CLUBB_core/parameter_indices.F90
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)atmos_param/clubb/CLUBB_core/parameter_indices.F90
parameters_microphys.o: $(SRCROOT)atmos_param/clubb/CLUBB_core/parameters_microphys.F90 clubb_precision.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)atmos_param/clubb/CLUBB_core/parameters_microphys.F90
parameters_model.o: $(SRCROOT)atmos_param/clubb/CLUBB_core/parameters_model.F90 clubb_precision.o constants_clubb.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)atmos_param/clubb/CLUBB_core/parameters_model.F90
parameters_radiation.o: $(SRCROOT)atmos_param/clubb/CLUBB_core/parameters_radiation.F90 clubb_precision.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)atmos_param/clubb/CLUBB_core/parameters_radiation.F90
parameters_tunable.o: $(SRCROOT)atmos_param/clubb/CLUBB_core/parameters_tunable.F90 parameter_indices.o grid_class.o clubb_precision.o constants_clubb.o error_code.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)atmos_param/clubb/CLUBB_core/parameters_tunable.F90
pdf_closure_module.o: $(SRCROOT)atmos_param/clubb/CLUBB_core/pdf_closure_module.F90 constants_clubb.o parameters_model.o parameters_tunable.o pdf_parameter_module.o anl_erf.o numerical_check.o saturation.o error_code.o stats_variables.o clubb_precision.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)atmos_param/clubb/CLUBB_core/pdf_closure_module.F90
pdf_parameter_module.o: $(SRCROOT)atmos_param/clubb/CLUBB_core/pdf_parameter_module.F90 clubb_precision.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)atmos_param/clubb/CLUBB_core/pdf_parameter_module.F90
pf_to_mr.o: $(SRCROOT)atmos_param/cosp/llnl/pf_to_mr.F /home/wenchang/tigress/CM2.5_tiger2_wy_killtc/exp/../src/atmos_param/cosp/cosp_defs.H
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -I$(SRCROOT)atmos_param/cosp	$(SRCROOT)atmos_param/cosp/llnl/pf_to_mr.F
physics_driver.o: $(SRCROOT)atmos_param/physics_driver/physics_driver.F90 time_manager.o field_manager.o tracer_manager.o atmos_tracer_driver.o mpp.o fms.o fms_io.o constants.o diag_manager.o rad_utilities.o cosp_driver.o moist_processes.o moistproc_kernels.o vert_turb_driver.o vert_diff_driver.o radiation_driver.o cloud_spec.o aerosol.o radiative_gases.o damping_driver.o grey_radiation.o CLUBB_driver_SCM.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)atmos_param/physics_driver/physics_driver.F90
platform.o: $(SRCROOT)shared/platform/platform.F90 /home/wenchang/tigress/CM2.5_tiger2_wy_killtc/exp/../src/shared/include/fms_platform.h
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -I$(SRCROOT)shared/include	$(SRCROOT)shared/platform/platform.F90
polysvp.o: $(SRCROOT)atmos_param/strat_cloud/polysvp.F90 fms.o sat_vapor_pres.o constants.o mg_const.o strat_cloud_utilities.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)atmos_param/strat_cloud/polysvp.F90
pos_definite_module.o: $(SRCROOT)atmos_param/clubb/CLUBB_core/pos_definite_module.F90 grid_class.o constants_clubb.o clubb_precision.o error_code.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)atmos_param/clubb/CLUBB_core/pos_definite_module.F90
potentialml.o: $(SRCROOT)mom5/ocean_param/gotm-4.0/turbulence/potentialml.F90 /home/wenchang/tigress/CM2.5_tiger2_wy_killtc/exp/../src/mom5/ocean_param/gotm-4.0/include/cppdefs.h /home/wenchang/tigress/CM2.5_tiger2_wy_killtc/exp/../src/mom5/ocean_param/gotm-4.0/include/version.h turbulence.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -I$(SRCROOT)mom5/ocean_param/gotm-4.0/include	$(SRCROOT)mom5/ocean_param/gotm-4.0/turbulence/potentialml.F90
prec_scops.o: $(SRCROOT)atmos_param/cosp/llnl/prec_scops.F /home/wenchang/tigress/CM2.5_tiger2_wy_killtc/exp/../src/atmos_param/cosp/cosp_defs.H
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -I$(SRCROOT)atmos_param/cosp	$(SRCROOT)atmos_param/cosp/llnl/prec_scops.F
production.o: $(SRCROOT)mom5/ocean_param/gotm-4.0/turbulence/production.F90 /home/wenchang/tigress/CM2.5_tiger2_wy_killtc/exp/../src/mom5/ocean_param/gotm-4.0/include/cppdefs.h /home/wenchang/tigress/CM2.5_tiger2_wy_killtc/exp/../src/mom5/ocean_param/gotm-4.0/include/version.h turbulence.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -I$(SRCROOT)mom5/ocean_param/gotm-4.0/include	$(SRCROOT)mom5/ocean_param/gotm-4.0/turbulence/production.F90
q2over2eq.o: $(SRCROOT)mom5/ocean_param/gotm-4.0/turbulence/q2over2eq.F90 /home/wenchang/tigress/CM2.5_tiger2_wy_killtc/exp/../src/mom5/ocean_param/gotm-4.0/include/cppdefs.h /home/wenchang/tigress/CM2.5_tiger2_wy_killtc/exp/../src/mom5/ocean_param/gotm-4.0/include/version.h turbulence.o util.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -I$(SRCROOT)mom5/ocean_param/gotm-4.0/include	$(SRCROOT)mom5/ocean_param/gotm-4.0/turbulence/q2over2eq.F90
quicksort.o: $(SRCROOT)shared/drifters/quicksort.F90
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)shared/drifters/quicksort.F90
r_ratio.o: $(SRCROOT)mom5/ocean_param/gotm-4.0/turbulence/r_ratio.F90 /home/wenchang/tigress/CM2.5_tiger2_wy_killtc/exp/../src/mom5/ocean_param/gotm-4.0/include/cppdefs.h /home/wenchang/tigress/CM2.5_tiger2_wy_killtc/exp/../src/mom5/ocean_param/gotm-4.0/include/version.h turbulence.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -I$(SRCROOT)mom5/ocean_param/gotm-4.0/include	$(SRCROOT)mom5/ocean_param/gotm-4.0/turbulence/r_ratio.F90
rad_diag.o: $(SRCROOT)atmos_param/fsrad/rad_diag.F90 rdparm.o hconst.o longwave.o fms.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)atmos_param/fsrad/rad_diag.F90
rad_output_file.o: $(SRCROOT)atmos_param/sea_esf_rad/rad_output_file.F90 mpp.o fms.o time_manager.o diag_manager.o diag_data.o constants.o rad_utilities.o esfsw_parameters.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)atmos_param/sea_esf_rad/rad_output_file.F90
rad_utilities.o: $(SRCROOT)atmos_param/sea_esf_rad/rad_utilities.F90 mpp.o fms.o field_manager.o time_manager.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)atmos_param/sea_esf_rad/rad_utilities.F90
radar_simulator.o: $(SRCROOT)atmos_param/cosp/quickbeam/radar_simulator.F90 /home/wenchang/tigress/CM2.5_tiger2_wy_killtc/exp/../src/atmos_param/cosp/cosp_defs.H mrgrnk.o array_lib.o math_lib.o optics_lib.o radar_simulator_types.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -I$(SRCROOT)atmos_param/cosp	$(SRCROOT)atmos_param/cosp/quickbeam/radar_simulator.F90
radar_simulator_types.o: $(SRCROOT)atmos_param/cosp/quickbeam/radar_simulator_types.F90 /home/wenchang/tigress/CM2.5_tiger2_wy_killtc/exp/../src/atmos_param/cosp/cosp_defs.H
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -I$(SRCROOT)atmos_param/cosp	$(SRCROOT)atmos_param/cosp/quickbeam/radar_simulator_types.F90
radiation_diag.o: $(SRCROOT)atmos_param/sea_esf_rad/radiation_diag.F90 mpp.o fms.o constants.o rad_utilities.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)atmos_param/sea_esf_rad/radiation_diag.F90
radiation_driver.o: $(SRCROOT)atmos_param/radiation_driver/radiation_driver.F90 mpp.o fms.o fms_io.o diag_manager.o diag_data.o time_manager.o sat_vapor_pres.o constants.o data_override.o rad_utilities.o esfsw_parameters.o diag_integral.o astronomy.o original_fms_rad.o sea_esf_rad.o rad_output_file.o cloudrad_package.o cloudrad_diagnostics.o microphys_rad.o aerosolrad_package.o field_manager.o tracer_manager.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)atmos_param/radiation_driver/radiation_driver.F90
radiative_gases.o: $(SRCROOT)atmos_param/sea_esf_rad/radiative_gases.F90 time_manager.o diag_manager.o mpp.o fms.o fms_io.o time_interp.o tracer_manager.o field_manager.o rad_utilities.o ozone.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)atmos_param/sea_esf_rad/radiative_gases.F90
random_numbers.o: $(SRCROOT)shared/random_numbers/random_numbers.F90 MersenneTwister.o time_manager.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)shared/random_numbers/random_numbers.F90
ras.o: $(SRCROOT)atmos_param/ras/ras.F90 mpp.o sat_vapor_pres.o constants.o diag_manager.o time_manager.o fms.o field_manager.o tracer_manager.o rad_utilities.o aer_ccn_act.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)atmos_param/ras/ras.F90
rdparm.o: $(SRCROOT)atmos_param/fsrad/rdparm.F90 fms.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)atmos_param/fsrad/rdparm.F90
read_climate_nudge_data.o: $(SRCROOT)atmos_cubed_sphere/tools/read_climate_nudge_data.F90 fms.o mpp.o mpp_io.o constants.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)atmos_cubed_sphere/tools/read_climate_nudge_data.F90
read_mosaic.o: $(SRCROOT)shared/mosaic/read_mosaic.c /home/wenchang/tigress/CM2.5_tiger2_wy_killtc/exp/../src/shared/mosaic/read_mosaic.h /home/wenchang/tigress/CM2.5_tiger2_wy_killtc/exp/../src/shared/mosaic/constant.h /home/wenchang/tigress/CM2.5_tiger2_wy_killtc/exp/../src/shared/mosaic/mosaic_util.h
	$(CC) $(CPPDEFS) $(CPPFLAGS) $(CFLAGS) $(OTHERFLAGS) $(OTHER_CFLAGS) -c -I$(SRCROOT)shared/mosaic	$(SRCROOT)shared/mosaic/read_mosaic.c
rh_based_clouds.o: $(SRCROOT)atmos_param/sea_esf_rad/rh_based_clouds.F90 mpp.o fms.o rh_clouds.o rad_utilities.o constants.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)atmos_param/sea_esf_rad/rh_based_clouds.F90
rh_clouds.o: $(SRCROOT)atmos_param/rh_clouds/rh_clouds.F90 mpp.o fms.o fms_io.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)atmos_param/rh_clouds/rh_clouds.F90
river.o: $(SRCROOT)land_lad2/river/river.F90 mpp.o fms.o mpp_domains.o fms_io.o diag_manager.o time_manager.o river_type.o river_physics.o constants.o stock_constants.o land_tile.o land_data.o lake_tile.o mpp_io.o grid.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)land_lad2/river/river.F90
river_physics.o: $(SRCROOT)land_lad2/river/river_physics.F90 mpp.o fms.o mpp_domains.o diag_manager.o river_type.o lake.o lake_tile.o constants.o land_debug.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)land_lad2/river/river_physics.F90
river_type.o: $(SRCROOT)land_lad2/river/river_type.F90 time_manager.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)land_lad2/river/river_type.F90
rot.o: $(SRCROOT)ice_sis/rot.F90 mpp.o axis_utils.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)ice_sis/rot.F90
rotstayn_klein_mp.o: $(SRCROOT)atmos_param/strat_cloud/rotstayn_klein_mp.F90 mpp.o fms.o cloud_generator.o constants.o aer_in_act.o sat_vapor_pres.o polysvp.o strat_cloud_utilities.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)atmos_param/strat_cloud/rotstayn_klein_mp.F90
sat_vapor_pres.o: $(SRCROOT)shared/sat_vapor_pres/sat_vapor_pres.F90 constants.o fms.o mpp_io.o mpp.o sat_vapor_pres_k.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)shared/sat_vapor_pres/sat_vapor_pres.F90
sat_vapor_pres_k.o: $(SRCROOT)shared/sat_vapor_pres/sat_vapor_pres_k.F90
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)shared/sat_vapor_pres/sat_vapor_pres_k.F90
saturation.o: $(SRCROOT)atmos_param/clubb/CLUBB_core/saturation.F90 model_flags.o clubb_precision.o constants_clubb.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)atmos_param/clubb/CLUBB_core/saturation.F90
scops.o: $(SRCROOT)atmos_param/cosp/icarus-scops-3.7/scops.F /home/wenchang/tigress/CM2.5_tiger2_wy_killtc/exp/../src/atmos_param/cosp/cosp_defs.H /home/wenchang/tigress/CM2.5_tiger2_wy_killtc/exp/../src/atmos_param/cosp/icarus-scops-3.7/congvec.H
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -I$(SRCROOT)atmos_param/cosp -I$(SRCROOT)atmos_param/cosp/icarus-scops-3.7	$(SRCROOT)atmos_param/cosp/icarus-scops-3.7/scops.F
sea_esf_rad.o: $(SRCROOT)atmos_param/sea_esf_rad/sea_esf_rad.F90 mpp.o fms.o time_manager.o rad_utilities.o radiation_diag.o longwave_driver.o shortwave_driver.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)atmos_param/sea_esf_rad/sea_esf_rad.F90
sealw99.o: $(SRCROOT)atmos_param/sea_esf_rad/sealw99.F90 mpp.o fms.o time_manager.o constants.o rad_utilities.o longwave_params.o longwave_clouds.o longwave_fluxes.o longwave_tables.o optical_path.o gas_tf.o lw_gases_stdtf.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)atmos_param/sea_esf_rad/sealw99.F90
shallow_conv.o: $(SRCROOT)atmos_param/shallow_conv/shallow_conv.F90 sat_vapor_pres.o mpp.o fms.o constants.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)atmos_param/shallow_conv/shallow_conv.F90
shortwave.o: $(SRCROOT)atmos_param/fsrad/shortwave.F90 rdparm.o hconst.o fms.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)atmos_param/fsrad/shortwave.F90
shortwave_driver.o: $(SRCROOT)atmos_param/sea_esf_rad/shortwave_driver.F90 mpp.o fms.o rad_utilities.o esfsw_parameters.o lhsw_driver.o esfsw_driver.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)atmos_param/sea_esf_rad/shortwave_driver.F90
sigma_sqd_w_module.o: $(SRCROOT)atmos_param/clubb/CLUBB_core/sigma_sqd_w_module.F90 constants_clubb.o clubb_precision.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)atmos_param/clubb/CLUBB_core/sigma_sqd_w_module.F90
sim_nc_mod.o: $(SRCROOT)atmos_cubed_sphere/tools/sim_nc_mod.F90 mpp.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)atmos_cubed_sphere/tools/sim_nc_mod.F90
simple_pdf.o: $(SRCROOT)atmos_param/strat_cloud/simple_pdf.F90 fms.o betaDistribution.o strat_cloud_utilities.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)atmos_param/strat_cloud/simple_pdf.F90
snow.o: $(SRCROOT)land_lad2/snow/snow.F90 mpp.o fms.o fms_io.o time_manager.o constants.o land_constants.o snow_tile.o land_tile.o land_tile_diag.o land_data.o land_tile_io.o nf_utils.o land_debug.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)land_lad2/snow/snow.F90
snow_tile.o: $(SRCROOT)land_lad2/snow/snow_tile.F90 /home/wenchang/tigress/CM2.5_tiger2_wy_killtc/exp/../src/shared/include/fms_platform.h mpp.o fms.o constants.o land_constants.o land_tile_diag_sel.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -I$(SRCROOT)shared/include	$(SRCROOT)land_lad2/snow/snow_tile.F90
soil.o: $(SRCROOT)land_lad2/soil/soil.F90 /home/wenchang/tigress/CM2.5_tiger2_wy_killtc/exp/../src/land_lad2/canopy_air/../shared/debug.inc mpp.o fms.o fms_io.o time_manager.o diag_manager.o constants.o land_constants.o soil_tile.o land_tile.o land_utils.o land_tile_diag.o land_data.o land_io.o land_tile_io.o nf_utils.o vegn_tile.o land_debug.o uptake.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -I$(SRCROOT)land_lad2/canopy_air/../shared	$(SRCROOT)land_lad2/soil/soil.F90
soil_tile.o: $(SRCROOT)land_lad2/soil/soil_tile.F90 /home/wenchang/tigress/CM2.5_tiger2_wy_killtc/exp/../src/shared/include/fms_platform.h mpp.o fms.o mpp_io.o constants.o land_constants.o land_io.o land_tile_diag_sel.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -I$(SRCROOT)shared/include	$(SRCROOT)land_lad2/soil/soil_tile.F90
sorted_index.o: $(SRCROOT)atmos_cubed_sphere/tools/sorted_index.F90
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)atmos_cubed_sphere/tools/sorted_index.F90
specified_clouds_W.o: $(SRCROOT)atmos_param/sea_esf_rad/specified_clouds_W.F90 time_manager.o cloud_zonal.o cloud_obs.o mpp.o fms.o rad_utilities.o constants.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)atmos_param/sea_esf_rad/specified_clouds_W.F90
sphum.o: $(SRCROOT)land_lad2/shared/sphum.F90 constants.o sat_vapor_pres.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)land_lad2/shared/sphum.F90
sponge_layer_damping.o: $(SRCROOT)atmos_param/clubb/CLUBB_core/sponge_layer_damping.F90 clubb_precision.o grid_class.o constants_clubb.o interpolation.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)atmos_param/clubb/CLUBB_core/sponge_layer_damping.F90
stable_bl_turb.o: $(SRCROOT)atmos_param/stable_bl_turb/stable_bl_turb.F90 mpp.o fms.o diag_manager.o time_manager.o constants.o monin_obukhov.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)atmos_param/stable_bl_turb/stable_bl_turb.F90
standalone_clouds.o: $(SRCROOT)atmos_param/sea_esf_rad/standalone_clouds.F90 mpp.o fms.o rad_utilities.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)atmos_param/sea_esf_rad/standalone_clouds.F90
stat_file_module.o: $(SRCROOT)atmos_param/clubb/CLUBB_core/stat_file_module.F90 clubb_precision.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)atmos_param/clubb/CLUBB_core/stat_file_module.F90
stats_LH_zt.o: $(SRCROOT)atmos_param/clubb/CLUBB_core/stats_LH_zt.F90 constants_clubb.o stats_variables.o stats_type.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)atmos_param/clubb/CLUBB_core/stats_LH_zt.F90
stats_rad_zm.o: $(SRCROOT)atmos_param/clubb/CLUBB_core/stats_rad_zm.F90 constants_clubb.o stats_variables.o stats_type.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)atmos_param/clubb/CLUBB_core/stats_rad_zm.F90
stats_rad_zt.o: $(SRCROOT)atmos_param/clubb/CLUBB_core/stats_rad_zt.F90 constants_clubb.o stats_variables.o stats_type.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)atmos_param/clubb/CLUBB_core/stats_rad_zt.F90
stats_sfc.o: $(SRCROOT)atmos_param/clubb/CLUBB_core/stats_sfc.F90 constants_clubb.o stats_variables.o stats_type.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)atmos_param/clubb/CLUBB_core/stats_sfc.F90
stats_subs.o: $(SRCROOT)atmos_param/clubb/CLUBB_core/stats_subs.F90 stats_variables.o clubb_precision.o output_grads.o output_netcdf.o stats_zm.o stats_zt.o stats_LH_zt.o stats_rad_zt.o stats_rad_zm.o stats_sfc.o error_code.o constants_clubb.o parameters_microphys.o grid_class.o variables_diagnostic_module.o pdf_parameter_module.o T_in_K_module.o parameters_model.o stats_type.o fill_holes.o interpolation.o saturation.o array_index.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)atmos_param/clubb/CLUBB_core/stats_subs.F90
stats_type.o: $(SRCROOT)atmos_param/clubb/CLUBB_core/stats_type.F90 stat_file_module.o clubb_precision.o grid_class.o error_code.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)atmos_param/clubb/CLUBB_core/stats_type.F90
stats_variables.o: $(SRCROOT)atmos_param/clubb/CLUBB_core/stats_variables.F90 stats_type.o clubb_precision.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)atmos_param/clubb/CLUBB_core/stats_variables.F90
stats_zm.o: $(SRCROOT)atmos_param/clubb/CLUBB_core/stats_zm.F90 constants_clubb.o stats_variables.o stats_type.o parameters_model.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)atmos_param/clubb/CLUBB_core/stats_zm.F90
stats_zt.o: $(SRCROOT)atmos_param/clubb/CLUBB_core/stats_zt.F90 constants_clubb.o stats_variables.o stats_type.o parameters_model.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)atmos_param/clubb/CLUBB_core/stats_zt.F90
stock_constants.o: $(SRCROOT)shared/exchange/stock_constants.F90 mpp.o time_manager.o diag_manager.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)shared/exchange/stock_constants.F90
strat_chem_driver.o: $(SRCROOT)atmos_shared/tracer_driver/stratchem/strat_chem_driver.F90 mpp.o fms.o mpp_io.o tracer_manager.o field_manager.o constants.o strat_chem_model.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)atmos_shared/tracer_driver/stratchem/strat_chem_driver.F90
strat_chem_model.o: $(SRCROOT)atmos_shared/tracer_driver/stratchem/strat_chem_model.F90
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)atmos_shared/tracer_driver/stratchem/strat_chem_model.F90
strat_chem_utilities.o: $(SRCROOT)atmos_shared/tracer_driver/tropchem/strat_chem_utilities.F90 mpp_io.o mpp.o fms.o constants.o time_manager.o interpolator.o time_interp.o diag_manager.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)atmos_shared/tracer_driver/tropchem/strat_chem_utilities.F90
strat_cloud.o: $(SRCROOT)atmos_param/strat_cloud/strat_cloud.F90 /home/wenchang/tigress/CM2.5_tiger2_wy_killtc/exp/../src/atmos_param/strat_cloud/strat_nml.h mpp.o fms.o fms_io.o constants.o cloud_rad.o diag_manager.o time_manager.o rad_utilities.o microphysics.o nc_cond.o polysvp.o aerosol_cloud.o strat_cloud_utilities.o strat_netcdf.o strat_cloud_legacy.o check_nan.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -I$(SRCROOT)atmos_param/strat_cloud	$(SRCROOT)atmos_param/strat_cloud/strat_cloud.F90
strat_cloud_legacy.o: $(SRCROOT)atmos_param/strat_cloud/strat_cloud_legacy.F90 sat_vapor_pres.o betaDistribution.o aer_ccn_act.o aer_in_act.o mpp.o fms.o fms_io.o constants.o cloud_rad.o diag_manager.o time_manager.o cloud_generator.o rad_utilities.o strat_cloud_utilities.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)atmos_param/strat_cloud/strat_cloud_legacy.F90
strat_cloud_utilities.o: $(SRCROOT)atmos_param/strat_cloud/strat_cloud_utilities.F90 fms.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)atmos_param/strat_cloud/strat_cloud_utilities.F90
strat_clouds_W.o: $(SRCROOT)atmos_param/sea_esf_rad/strat_clouds_W.F90 constants.o time_manager.o mpp.o fms.o rad_utilities.o esfsw_parameters.o cloud_rad.o random_numbers.o cloud_generator.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)atmos_param/sea_esf_rad/strat_clouds_W.F90
strat_netcdf.o: $(SRCROOT)atmos_param/strat_cloud/strat_netcdf.F90 fms.o diag_manager.o time_manager.o strat_cloud_utilities.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)atmos_param/strat_cloud/strat_netcdf.F90
surface_flux.o: $(SRCROOT)coupler/surface_flux.F90 fms.o monin_obukhov.o sat_vapor_pres.o constants.o mpp.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)coupler/surface_flux.F90
surface_varnce_module.o: $(SRCROOT)atmos_param/clubb/CLUBB_core/surface_varnce_module.F90 parameters_model.o constants_clubb.o numerical_check.o error_code.o array_index.o stats_type.o clubb_precision.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)atmos_param/clubb/CLUBB_core/surface_varnce_module.F90
sw_core.o: $(SRCROOT)atmos_cubed_sphere/model/sw_core.F90 fv_mp_mod.o tp_core.o fv_arrays.o a2b_edge.o test_cases.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)atmos_cubed_sphere/model/sw_core.F90
table_printer.o: $(SRCROOT)land_lad2/shared/table_printer.F90
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)land_lad2/shared/table_printer.F90
test_cases.o: $(SRCROOT)atmos_cubed_sphere/tools/test_cases.F90 constants.o init_hydro.o fv_mp_mod.o fv_grid_utils.o fv_surf_map.o fv_grid_tools.o fv_eta.o mpp.o mpp_domains.o mpp_parameter.o fv_sg.o fv_diagnostics.o fv_arrays.o lin_cloud_microphys.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)atmos_cubed_sphere/tools/test_cases.F90
test_fms_io.o: $(SRCROOT)shared/fms/test_fms_io.F90 /home/wenchang/tigress/CM2.5_tiger2_wy_killtc/exp/../src/shared/include/fms_platform.h mpp.o mpp_domains.o mpp_io.o fms_io.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -I$(SRCROOT)shared/include	$(SRCROOT)shared/fms/test_fms_io.F90
test_horiz_interp.o: $(SRCROOT)shared/horiz_interp/test_horiz_interp.F90 mpp.o mpp_io.o mpp_domains.o horiz_interp.o axis_utils.o fms_io.o fms.o constants.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)shared/horiz_interp/test_horiz_interp.F90
test_mpp.o: $(SRCROOT)shared/mpp/test_mpp.F90 /home/wenchang/tigress/CM2.5_tiger2_wy_killtc/exp/../src/shared/include/fms_platform.h mpp.o mpp_io.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -I$(SRCROOT)shared/include	$(SRCROOT)shared/mpp/test_mpp.F90
test_mpp_domains.o: $(SRCROOT)shared/mpp/test_mpp_domains.F90 /home/wenchang/tigress/CM2.5_tiger2_wy_killtc/exp/../src/shared/include/fms_platform.h mpp.o mpp_domains.o mpp_memutils.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -I$(SRCROOT)shared/include	$(SRCROOT)shared/mpp/test_mpp_domains.F90
test_mpp_io.o: $(SRCROOT)shared/mpp/test_mpp_io.F90 /home/wenchang/tigress/CM2.5_tiger2_wy_killtc/exp/../src/shared/include/fms_platform.h mpp.o mpp_domains.o mpp_io.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -I$(SRCROOT)shared/include	$(SRCROOT)shared/mpp/test_mpp_io.F90
test_mpp_pset.o: $(SRCROOT)shared/mpp/test_mpp_pset.F90 /home/wenchang/tigress/CM2.5_tiger2_wy_killtc/exp/../src/shared/include/fms_platform.h mpp.o mpp_pset.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -I$(SRCROOT)shared/include	$(SRCROOT)shared/mpp/test_mpp_pset.F90
test_xgrid.o: $(SRCROOT)shared/exchange/test_xgrid.F90 mpp.o mpp_domains.o mpp_io.o fms.o fms_io.o constants.o xgrid.o mosaic.o gradient.o ensemble_manager.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)shared/exchange/test_xgrid.F90
threadloc.o: $(SRCROOT)shared/mpp/threadloc.c
	$(CC) $(CPPDEFS) $(CPPFLAGS) $(CFLAGS) $(OTHERFLAGS) $(OTHER_CFLAGS) -c	$(SRCROOT)shared/mpp/threadloc.c
time.o: $(SRCROOT)mom5/ocean_param/gotm-4.0/util/time.F90 /home/wenchang/tigress/CM2.5_tiger2_wy_killtc/exp/../src/mom5/ocean_param/gotm-4.0/include/cppdefs.h /home/wenchang/tigress/CM2.5_tiger2_wy_killtc/exp/../src/mom5/ocean_param/gotm-4.0/include/version.h
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -I$(SRCROOT)mom5/ocean_param/gotm-4.0/include	$(SRCROOT)mom5/ocean_param/gotm-4.0/util/time.F90
time_interp.o: $(SRCROOT)shared/time_interp/time_interp.F90 time_manager.o fms.o mpp.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)shared/time_interp/time_interp.F90
time_interp_external.o: $(SRCROOT)shared/time_interp/time_interp_external.F90 /home/wenchang/tigress/CM2.5_tiger2_wy_killtc/exp/../src/shared/include/fms_platform.h mpp.o mpp_io.o time_manager.o get_cal_time.o mpp_domains.o time_interp.o axis_utils.o fms.o platform.o horiz_interp.o constants.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -I$(SRCROOT)shared/include	$(SRCROOT)shared/time_interp/time_interp_external.F90
time_manager.o: $(SRCROOT)shared/time_manager/time_manager.F90 /home/wenchang/tigress/CM2.5_tiger2_wy_killtc/exp/../src/shared/include/fms_platform.h constants.o fms.o mpp.o fms_io.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -I$(SRCROOT)shared/include	$(SRCROOT)shared/time_manager/time_manager.F90
tke_turb.o: $(SRCROOT)atmos_param/tke_turb/tke_turb.F90 mpp.o fms.o constants.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)atmos_param/tke_turb/tke_turb.F90
tkealgebraic.o: $(SRCROOT)mom5/ocean_param/gotm-4.0/turbulence/tkealgebraic.F90 /home/wenchang/tigress/CM2.5_tiger2_wy_killtc/exp/../src/mom5/ocean_param/gotm-4.0/include/cppdefs.h /home/wenchang/tigress/CM2.5_tiger2_wy_killtc/exp/../src/mom5/ocean_param/gotm-4.0/include/version.h turbulence.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -I$(SRCROOT)mom5/ocean_param/gotm-4.0/include	$(SRCROOT)mom5/ocean_param/gotm-4.0/turbulence/tkealgebraic.F90
tkeeq.o: $(SRCROOT)mom5/ocean_param/gotm-4.0/turbulence/tkeeq.F90 /home/wenchang/tigress/CM2.5_tiger2_wy_killtc/exp/../src/mom5/ocean_param/gotm-4.0/include/cppdefs.h /home/wenchang/tigress/CM2.5_tiger2_wy_killtc/exp/../src/mom5/ocean_param/gotm-4.0/include/version.h turbulence.o util.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -I$(SRCROOT)mom5/ocean_param/gotm-4.0/include	$(SRCROOT)mom5/ocean_param/gotm-4.0/turbulence/tkeeq.F90
topo_drag.o: $(SRCROOT)atmos_param/topo_drag/topo_drag.F90 mpp.o fms.o fms_io.o constants.o horiz_interp.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)atmos_param/topo_drag/topo_drag.F90
topo_rough.o: $(SRCROOT)land_lad2/topo_rough/topo_rough.F90 time_manager.o mpp_domains.o mpp.o fms.o diag_manager.o topography.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)land_lad2/topo_rough/topo_rough.F90
topography.o: $(SRCROOT)shared/topography/topography.F90 gaussian_topog.o horiz_interp.o fms.o fms_io.o constants.o mpp.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)shared/topography/topography.F90
tp_core.o: $(SRCROOT)atmos_cubed_sphere/model/tp_core.F90 fv_mp_mod.o fv_grid_utils.o fv_arrays.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)atmos_cubed_sphere/model/tp_core.F90
tracer_manager.o: $(SRCROOT)shared/tracer_manager/tracer_manager.F90 mpp.o mpp_io.o fms.o field_manager.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)shared/tracer_manager/tracer_manager.F90
transitions.o: $(SRCROOT)land_lad2/transitions/transitions.F90 /home/wenchang/tigress/CM2.5_tiger2_wy_killtc/exp/../src/shared/include/fms_platform.h /home/wenchang/tigress/CM2.5_tiger2_wy_killtc/exp/../src/land_lad2/canopy_air/../shared/debug.inc constants.o mpp.o fms.o mpp_io.o axis_utils.o time_manager.o get_cal_time.o horiz_interp.o time_interp.o diag_manager.o nfu.o vegn_data.o cana_tile.o snow_tile.o vegn_tile.o soil_tile.o land_tile.o land_tile_io.o land_data.o vegn_harvesting.o land_debug.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -I$(SRCROOT)shared/include -I$(SRCROOT)land_lad2/canopy_air/../shared	$(SRCROOT)land_lad2/transitions/transitions.F90
transport_matrix.o: $(SRCROOT)mom5/ocean_tracers/transport_matrix.F90 /home/wenchang/tigress/CM2.5_tiger2_wy_killtc/exp/../src/shared/include/fms_platform.h field_manager.o mpp.o diag_manager.o ocean_tpm_util.o ocean_types.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -I$(SRCROOT)shared/include	$(SRCROOT)mom5/ocean_tracers/transport_matrix.F90
tridiagonal.o: $(SRCROOT)shared/tridiagonal/tridiagonal.F90
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)shared/tridiagonal/tridiagonal.F90
tridiagonal_gotm.o: $(SRCROOT)mom5/ocean_param/gotm-4.0/util/tridiagonal_gotm.F90 /home/wenchang/tigress/CM2.5_tiger2_wy_killtc/exp/../src/mom5/ocean_param/gotm-4.0/include/cppdefs.h /home/wenchang/tigress/CM2.5_tiger2_wy_killtc/exp/../src/mom5/ocean_param/gotm-4.0/include/version.h
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -I$(SRCROOT)mom5/ocean_param/gotm-4.0/include	$(SRCROOT)mom5/ocean_param/gotm-4.0/util/tridiagonal_gotm.F90
tropchem_driver.o: $(SRCROOT)atmos_shared/tracer_driver/tropchem/tropchem_driver.F90 mpp.o fms.o time_manager.o diag_manager.o tracer_manager.o field_manager.o atmos_tracer_utilities.o constants.o interpolator.o time_interp.o mo_chemdr.o mo_chemini.o m_tracname.o moz.mods.o mo_hook.o strat_chem_utilities.o mo_chem_utls.o atmos_sulfate.o esfsw_parameters.o astronomy.o horiz_interp.o fms_io.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)atmos_shared/tracer_driver/tropchem/tropchem_driver.F90
turbulence.o: $(SRCROOT)mom5/ocean_param/gotm-4.0/turbulence/turbulence.F90 /home/wenchang/tigress/CM2.5_tiger2_wy_killtc/exp/../src/mom5/ocean_param/gotm-4.0/include/cppdefs.h /home/wenchang/tigress/CM2.5_tiger2_wy_killtc/exp/../src/mom5/ocean_param/gotm-4.0/include/version.h
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -I$(SRCROOT)mom5/ocean_param/gotm-4.0/include	$(SRCROOT)mom5/ocean_param/gotm-4.0/turbulence/turbulence.F90
uptake.o: $(SRCROOT)land_lad2/soil/uptake.F90 /home/wenchang/tigress/CM2.5_tiger2_wy_killtc/exp/../src/land_lad2/canopy_air/../shared/debug.inc constants.o fms.o soil_tile.o land_debug.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -I$(SRCROOT)land_lad2/canopy_air/../shared	$(SRCROOT)land_lad2/soil/uptake.F90
util.o: $(SRCROOT)mom5/ocean_param/gotm-4.0/util/util.F90 /home/wenchang/tigress/CM2.5_tiger2_wy_killtc/exp/../src/mom5/ocean_param/gotm-4.0/include/cppdefs.h /home/wenchang/tigress/CM2.5_tiger2_wy_killtc/exp/../src/mom5/ocean_param/gotm-4.0/include/version.h
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -I$(SRCROOT)mom5/ocean_param/gotm-4.0/include	$(SRCROOT)mom5/ocean_param/gotm-4.0/util/util.F90
uw_clouds_W.o: $(SRCROOT)atmos_param/sea_esf_rad/uw_clouds_W.F90 time_manager.o mpp.o fms.o constants.o rad_utilities.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)atmos_param/sea_esf_rad/uw_clouds_W.F90
uw_conv.o: $(SRCROOT)atmos_param/shallow_cu/uw_conv.F90 /home/wenchang/tigress/CM2.5_tiger2_wy_killtc/exp/../src/shared/include/fms_platform.h mpp.o constants.o diag_manager.o time_manager.o fms.o field_manager.o tracer_manager.o sat_vapor_pres.o atmos_tracer_utilities.o rad_utilities.o aer_ccn_act.o conv_utilities.o conv_utilities_k.o conv_plumes_k.o conv_closures.o deep_conv.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -I$(SRCROOT)shared/include	$(SRCROOT)atmos_param/shallow_cu/uw_conv.F90
variables_diagnostic_module.o: $(SRCROOT)atmos_param/clubb/CLUBB_core/variables_diagnostic_module.F90 pdf_parameter_module.o clubb_precision.o constants_clubb.o parameters_model.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)atmos_param/clubb/CLUBB_core/variables_diagnostic_module.F90
variables_prognostic_module.o: $(SRCROOT)atmos_param/clubb/CLUBB_core/variables_prognostic_module.F90 pdf_parameter_module.o clubb_precision.o constants_clubb.o parameters_model.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)atmos_param/clubb/CLUBB_core/variables_prognostic_module.F90
variables_radiation_module.o: $(SRCROOT)atmos_param/clubb/CLUBB_core/variables_radiation_module.F90 clubb_precision.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)atmos_param/clubb/CLUBB_core/variables_radiation_module.F90
variances.o: $(SRCROOT)mom5/ocean_param/gotm-4.0/turbulence/variances.F90 /home/wenchang/tigress/CM2.5_tiger2_wy_killtc/exp/../src/mom5/ocean_param/gotm-4.0/include/cppdefs.h /home/wenchang/tigress/CM2.5_tiger2_wy_killtc/exp/../src/mom5/ocean_param/gotm-4.0/include/version.h turbulence.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -I$(SRCROOT)mom5/ocean_param/gotm-4.0/include	$(SRCROOT)mom5/ocean_param/gotm-4.0/turbulence/variances.F90
vegetation.o: $(SRCROOT)land_lad2/vegetation/vegetation.F90 /home/wenchang/tigress/CM2.5_tiger2_wy_killtc/exp/../src/land_lad2/canopy_air/../shared/debug.inc mpp.o fms.o fms_io.o time_manager.o constants.o sphum.o nf_utils.o vegn_tile.o soil_tile.o land_constants.o land_tile.o land_tile_diag.o land_data.o land_io.o land_tile_io.o vegn_data.o vegn_cohort.o canopy_air.o vegn_cohort_io.o land_debug.o vegn_radiation.o vegn_photosynthesis.o vegn_static_override.o vegn_dynamics.o vegn_disturbance.o vegn_harvesting.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -I$(SRCROOT)land_lad2/canopy_air/../shared	$(SRCROOT)land_lad2/vegetation/vegetation.F90
vegn_cohort.o: $(SRCROOT)land_lad2/vegetation/vegn_cohort.F90 constants.o land_constants.o vegn_data.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)land_lad2/vegetation/vegn_cohort.F90
vegn_cohort_io.o: $(SRCROOT)land_lad2/vegetation/vegn_cohort_io.F90 /home/wenchang/tigress/CM2.5_tiger2_wy_killtc/exp/../src/land_lad2/vegetation/vegn_cohort_io.inc fms.o fms_io.o mpp.o nf_utils.o land_io.o land_tile.o land_tile_io.o vegn_cohort.o land_data.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -I$(SRCROOT)land_lad2/vegetation	$(SRCROOT)land_lad2/vegetation/vegn_cohort_io.F90
vegn_data.o: $(SRCROOT)land_lad2/vegetation/vegn_data.F90 mpp.o fms.o land_constants.o land_tile_diag_sel.o table_printer.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)land_lad2/vegetation/vegn_data.F90
vegn_disturbance.o: $(SRCROOT)land_lad2/vegetation/vegn_disturbance.F90 land_constants.o vegn_data.o vegn_tile.o soil_tile.o vegn_cohort.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)land_lad2/vegetation/vegn_disturbance.F90
vegn_dynamics.o: $(SRCROOT)land_lad2/vegetation/vegn_dynamics.F90 /home/wenchang/tigress/CM2.5_tiger2_wy_killtc/exp/../src/land_lad2/canopy_air/../shared/debug.inc fms.o time_manager.o land_constants.o land_tile_diag.o vegn_data.o vegn_tile.o soil_tile.o vegn_cohort.o land_debug.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -I$(SRCROOT)land_lad2/canopy_air/../shared	$(SRCROOT)land_lad2/vegetation/vegn_dynamics.F90
vegn_harvesting.o: $(SRCROOT)land_lad2/vegetation/vegn_harvesting.F90 mpp.o fms.o mpp_io.o vegn_data.o vegn_tile.o vegn_cohort.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)land_lad2/vegetation/vegn_harvesting.F90
vegn_photosynthesis.o: $(SRCROOT)land_lad2/vegetation/vegn_photosynthesis.F90 /home/wenchang/tigress/CM2.5_tiger2_wy_killtc/exp/../src/land_lad2/canopy_air/../shared/debug.inc fms.o constants.o sphum.o land_constants.o land_debug.o vegn_data.o vegn_tile.o vegn_cohort.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -I$(SRCROOT)land_lad2/canopy_air/../shared	$(SRCROOT)land_lad2/vegetation/vegn_photosynthesis.F90
vegn_radiation.o: $(SRCROOT)land_lad2/vegetation/vegn_radiation.F90 fms.o constants.o land_constants.o vegn_data.o vegn_tile.o vegn_cohort.o snow.o land_debug.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)land_lad2/vegetation/vegn_radiation.F90
vegn_static_override.o: $(SRCROOT)land_lad2/vegetation/vegn_static_override.F90 /home/wenchang/tigress/CM2.5_tiger2_wy_killtc/exp/../src/land_lad2/vegetation/read_remap_cohort_data.inc /home/wenchang/tigress/CM2.5_tiger2_wy_killtc/exp/../src/land_lad2/vegetation/read_remap_cohort_data_new.inc constants.o mpp.o mpp_io.o fms_io.o time_manager.o get_cal_time.o fms.o time_interp.o diag_manager.o nf_utils.o land_data.o land_io.o land_numerics.o land_tile_io.o land_tile.o vegn_cohort.o vegn_cohort_io.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -I$(SRCROOT)land_lad2/vegetation	$(SRCROOT)land_lad2/vegetation/vegn_static_override.F90
vegn_tile.o: $(SRCROOT)land_lad2/vegetation/vegn_tile.F90 fms.o constants.o land_constants.o land_io.o land_tile_diag_sel.o vegn_data.o vegn_cohort.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)land_lad2/vegetation/vegn_tile.F90
vert_advection.o: $(SRCROOT)atmos_shared/vert_advection/vert_advection.F90 fms.o mpp.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)atmos_shared/vert_advection/vert_advection.F90
vert_diff.o: $(SRCROOT)atmos_param/vert_diff/vert_diff.F90 constants.o fms.o field_manager.o tracer_manager.o mpp.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)atmos_param/vert_diff/vert_diff.F90
vert_diff_driver.o: $(SRCROOT)atmos_param/vert_diff_driver/vert_diff_driver.F90 vert_diff.o diag_manager.o time_manager.o mpp.o fms.o constants.o field_manager.o tracer_manager.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)atmos_param/vert_diff_driver/vert_diff_driver.F90
vert_turb_driver.o: $(SRCROOT)atmos_param/vert_turb_driver/vert_turb_driver.F90 my25_turb.o tke_turb.o diffusivity.o edt.o strat_cloud.o shallow_conv.o stable_bl_turb.o entrain.o diag_manager.o time_manager.o constants.o mpp.o fms.o field_manager.o tracer_manager.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)atmos_param/vert_turb_driver/vert_turb_driver.F90
wet_deposition_0D.o: $(SRCROOT)atmos_param/donner_deep/wet_deposition_0D.F90
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)atmos_param/donner_deep/wet_deposition_0D.F90
write_ocean_data.o: $(SRCROOT)shared/oda_tools/write_ocean_data.F90 mpp_io.o mpp.o oda_types.o time_manager.o mpp_domains.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)shared/oda_tools/write_ocean_data.F90
xbt_drop_rate_adjust.o: $(SRCROOT)shared/oda_tools/xbt_drop_rate_adjust.f90 oda_types.o
	$(FC) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)shared/oda_tools/xbt_drop_rate_adjust.f90
xgrid.o: $(SRCROOT)shared/exchange/xgrid.F90 /home/wenchang/tigress/CM2.5_tiger2_wy_killtc/exp/../src/shared/include/fms_platform.h fms.o fms_io.o mpp.o mpp_domains.o mpp_io.o constants.o mosaic.o stock_constants.o gradient.o time_manager.o diag_manager.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -I$(SRCROOT)shared/include	$(SRCROOT)shared/exchange/xgrid.F90
zeff.o: $(SRCROOT)atmos_param/cosp/quickbeam/zeff.F90 /home/wenchang/tigress/CM2.5_tiger2_wy_killtc/exp/../src/atmos_param/cosp/cosp_defs.H math_lib.o optics_lib.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -I$(SRCROOT)atmos_param/cosp	$(SRCROOT)atmos_param/cosp/quickbeam/zeff.F90
zetac_clouds_W.o: $(SRCROOT)atmos_param/sea_esf_rad/zetac_clouds_W.F90 mpp.o fms.o constants.o rad_utilities.o microphys_cloud.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)atmos_param/sea_esf_rad/zetac_clouds_W.F90
./adv_center.F90: $(SRCROOT)mom5/ocean_param/gotm-4.0/util/adv_center.F90
	cp $(SRCROOT)mom5/ocean_param/gotm-4.0/util/adv_center.F90 .
./kbalgebraic.F90: $(SRCROOT)mom5/ocean_param/gotm-4.0/turbulence/kbalgebraic.F90
	cp $(SRCROOT)mom5/ocean_param/gotm-4.0/turbulence/kbalgebraic.F90 .
./stats_zt.F90: $(SRCROOT)atmos_param/clubb/CLUBB_core/stats_zt.F90
	cp $(SRCROOT)atmos_param/clubb/CLUBB_core/stats_zt.F90 .
./ocmip2_abiotic.F90: $(SRCROOT)mom5/ocean_bgc/ocmip2_abiotic.F90
	cp $(SRCROOT)mom5/ocean_bgc/ocmip2_abiotic.F90 .
./radar_simulator_types.F90: $(SRCROOT)atmos_param/cosp/quickbeam/radar_simulator_types.F90
	cp $(SRCROOT)atmos_param/cosp/quickbeam/radar_simulator_types.F90 .
./m_tracname.F90: $(SRCROOT)atmos_shared/tracer_driver/tropchem/m_tracname.F90
	cp $(SRCROOT)atmos_shared/tracer_driver/tropchem/m_tracname.F90 .
./co2_data.F90: $(SRCROOT)atmos_param/fsrad/co2_data.F90
	cp $(SRCROOT)atmos_param/fsrad/co2_data.F90 .
./epsbalgebraic.F90: $(SRCROOT)mom5/ocean_param/gotm-4.0/turbulence/epsbalgebraic.F90
	cp $(SRCROOT)mom5/ocean_param/gotm-4.0/turbulence/epsbalgebraic.F90 .
./cldwat2m_micro.F90: $(SRCROOT)atmos_param/strat_cloud/cldwat2m_micro.F90
	cp $(SRCROOT)atmos_param/strat_cloud/cldwat2m_micro.F90 .
./atmos_dust.F90: $(SRCROOT)atmos_shared/tracer_driver/atmos_dust.F90
	cp $(SRCROOT)atmos_shared/tracer_driver/atmos_dust.F90 .
./ocean_tpm_util.F90: $(SRCROOT)mom5/ocean_tracers/ocean_tpm_util.F90
	cp $(SRCROOT)mom5/ocean_tracers/ocean_tpm_util.F90 .
./ocean_shortwave_gfdl.F90: $(SRCROOT)mom5/ocean_param/sources/ocean_shortwave_gfdl.F90
	cp $(SRCROOT)mom5/ocean_param/sources/ocean_shortwave_gfdl.F90 .
./monin_obukhov_interfaces.h: /home/wenchang/tigress/CM2.5_tiger2_wy_killtc/exp/../src/atmos_param/monin_obukhov/monin_obukhov_interfaces.h
	cp /home/wenchang/tigress/CM2.5_tiger2_wy_killtc/exp/../src/atmos_param/monin_obukhov/monin_obukhov_interfaces.h .
./mpp_pset.F90: $(SRCROOT)shared/mpp/mpp_pset.F90
	cp $(SRCROOT)shared/mpp/mpp_pset.F90 .
./drifters_io.F90: $(SRCROOT)shared/drifters/drifters_io.F90
	cp $(SRCROOT)shared/drifters/drifters_io.F90 .
./array_index.F90: $(SRCROOT)atmos_param/clubb/CLUBB_core/array_index.F90
	cp $(SRCROOT)atmos_param/clubb/CLUBB_core/array_index.F90 .
./algebraiclength.F90: $(SRCROOT)mom5/ocean_param/gotm-4.0/turbulence/algebraiclength.F90
	cp $(SRCROOT)mom5/ocean_param/gotm-4.0/turbulence/algebraiclength.F90 .
./time.F90: $(SRCROOT)mom5/ocean_param/gotm-4.0/util/time.F90
	cp $(SRCROOT)mom5/ocean_param/gotm-4.0/util/time.F90 .
./standalone_clouds.F90: $(SRCROOT)atmos_param/sea_esf_rad/standalone_clouds.F90
	cp $(SRCROOT)atmos_param/sea_esf_rad/standalone_clouds.F90 .
./fv_grid_tools.F90: $(SRCROOT)atmos_cubed_sphere/tools/fv_grid_tools.F90
	cp $(SRCROOT)atmos_cubed_sphere/tools/fv_grid_tools.F90 .
./test_xgrid.F90: $(SRCROOT)shared/exchange/test_xgrid.F90
	cp $(SRCROOT)shared/exchange/test_xgrid.F90 .
./ocean_grids.F90: $(SRCROOT)mom5/ocean_core/ocean_grids.F90
	cp $(SRCROOT)mom5/ocean_core/ocean_grids.F90 .
./bm_massflux.F90: $(SRCROOT)atmos_param/betts_miller/bm_massflux.F90
	cp $(SRCROOT)atmos_param/betts_miller/bm_massflux.F90 .
./mpp_comm_mpi.inc: /home/wenchang/tigress/CM2.5_tiger2_wy_killtc/exp/../src/shared/mpp/include/mpp_comm_mpi.inc
	cp /home/wenchang/tigress/CM2.5_tiger2_wy_killtc/exp/../src/shared/mpp/include/mpp_comm_mpi.inc .
./hydrostatic_module.F90: $(SRCROOT)atmos_param/clubb/CLUBB_core/hydrostatic_module.F90
	cp $(SRCROOT)atmos_param/clubb/CLUBB_core/hydrostatic_module.F90 .
./atmos_age_tracer.F90: $(SRCROOT)atmos_shared/tracer_driver/atmos_age_tracer.F90
	cp $(SRCROOT)atmos_shared/tracer_driver/atmos_age_tracer.F90 .
./optics_lib.F90: $(SRCROOT)atmos_param/cosp/quickbeam/optics_lib.F90
	cp $(SRCROOT)atmos_param/cosp/quickbeam/optics_lib.F90 .
./fv_dynamics.F90: $(SRCROOT)atmos_cubed_sphere/model/fv_dynamics.F90
	cp $(SRCROOT)atmos_cubed_sphere/model/fv_dynamics.F90 .
./cloud_interpolator.F90: $(SRCROOT)shared/drifters/cloud_interpolator.F90
	cp $(SRCROOT)shared/drifters/cloud_interpolator.F90 .
./variables_radiation_module.F90: $(SRCROOT)atmos_param/clubb/CLUBB_core/variables_radiation_module.F90
	cp $(SRCROOT)atmos_param/clubb/CLUBB_core/variables_radiation_module.F90 .
./interpolation.F90: $(SRCROOT)atmos_param/clubb/CLUBB_core/interpolation.F90
	cp $(SRCROOT)atmos_param/clubb/CLUBB_core/interpolation.F90 .
./betts_miller.F90: $(SRCROOT)atmos_param/betts_miller/betts_miller.F90
	cp $(SRCROOT)atmos_param/betts_miller/betts_miller.F90 .
./mpp_group_update.h: /home/wenchang/tigress/CM2.5_tiger2_wy_killtc/exp/../src/shared/mpp/include/mpp_group_update.h
	cp /home/wenchang/tigress/CM2.5_tiger2_wy_killtc/exp/../src/shared/mpp/include/mpp_group_update.h .
./fv_timing.F90: $(SRCROOT)atmos_cubed_sphere/tools/fv_timing.F90
	cp $(SRCROOT)atmos_cubed_sphere/tools/fv_timing.F90 .
./atmos_ocean_fluxes.F90: $(SRCROOT)shared/coupler/atmos_ocean_fluxes.F90
	cp $(SRCROOT)shared/coupler/atmos_ocean_fluxes.F90 .
./mpp_data.F90: $(SRCROOT)shared/mpp/mpp_data.F90
	cp $(SRCROOT)shared/mpp/mpp_data.F90 .
./canopy_air.F90: $(SRCROOT)land_lad2/canopy_air/canopy_air.F90
	cp $(SRCROOT)land_lad2/canopy_air/canopy_air.F90 .
./ocean_nphysics_skew.F90: $(SRCROOT)mom5/ocean_param/neutral/ocean_nphysics_skew.F90
	cp $(SRCROOT)mom5/ocean_param/neutral/ocean_nphysics_skew.F90 .
./amip_interp.F90: $(SRCROOT)shared/amip_interp/amip_interp.F90
	cp $(SRCROOT)shared/amip_interp/amip_interp.F90 .
./atmos_model.F90: $(SRCROOT)atmos_coupled/atmos_model.F90
	cp $(SRCROOT)atmos_coupled/atmos_model.F90 .
./sigma_sqd_w_module.F90: $(SRCROOT)atmos_param/clubb/CLUBB_core/sigma_sqd_w_module.F90
	cp $(SRCROOT)atmos_param/clubb/CLUBB_core/sigma_sqd_w_module.F90 .
./Skw_module.F90: $(SRCROOT)atmos_param/clubb/CLUBB_core/Skw_module.F90
	cp $(SRCROOT)atmos_param/clubb/CLUBB_core/Skw_module.F90 .
./gradient.F90: $(SRCROOT)shared/mosaic/gradient.F90
	cp $(SRCROOT)shared/mosaic/gradient.F90 .
./cosp_modis_simulator.F90: $(SRCROOT)atmos_param/cosp/cosp_modis_simulator.F90
	cp $(SRCROOT)atmos_param/cosp/cosp_modis_simulator.F90 .
./ocean_advection_velocity.F90: $(SRCROOT)mom5/ocean_core/ocean_advection_velocity.F90
	cp $(SRCROOT)mom5/ocean_core/ocean_advection_velocity.F90 .
./fv_surf_map.F90: $(SRCROOT)atmos_cubed_sphere/tools/fv_surf_map.F90
	cp $(SRCROOT)atmos_cubed_sphere/tools/fv_surf_map.F90 .
./advance_xm_wpxp_module.F90: $(SRCROOT)atmos_param/clubb/CLUBB_core/advance_xm_wpxp_module.F90
	cp $(SRCROOT)atmos_param/clubb/CLUBB_core/advance_xm_wpxp_module.F90 .
./fms_switches.h: /home/wenchang/tigress/CM2.5_tiger2_wy_killtc/exp/../src/shared/drifters/fms_switches.h
	cp /home/wenchang/tigress/CM2.5_tiger2_wy_killtc/exp/../src/shared/drifters/fms_switches.h .
./tridiagonal_gotm.F90: $(SRCROOT)mom5/ocean_param/gotm-4.0/util/tridiagonal_gotm.F90
	cp $(SRCROOT)mom5/ocean_param/gotm-4.0/util/tridiagonal_gotm.F90 .
./ocean_topog.F90: $(SRCROOT)mom5/ocean_core/ocean_topog.F90
	cp $(SRCROOT)mom5/ocean_core/ocean_topog.F90 .
./ocean_increment_eta.F90: $(SRCROOT)mom5/ocean_param/sources/ocean_increment_eta.F90
	cp $(SRCROOT)mom5/ocean_param/sources/ocean_increment_eta.F90 .
./FMS_coupler_util.F90: $(SRCROOT)ocean_shared/generic_tracers/FMS_coupler_util.F90
	cp $(SRCROOT)ocean_shared/generic_tracers/FMS_coupler_util.F90 .
./drifters_comm.F90: $(SRCROOT)shared/drifters/drifters_comm.F90
	cp $(SRCROOT)shared/drifters/drifters_comm.F90 .
./ocean_blob_static_bottom.F90: $(SRCROOT)mom5/ocean_blobs/ocean_blob_static_bottom.F90
	cp $(SRCROOT)mom5/ocean_blobs/ocean_blob_static_bottom.F90 .
./monin_obukhov_kernel.F90: $(SRCROOT)atmos_param/monin_obukhov/monin_obukhov_kernel.F90
	cp $(SRCROOT)atmos_param/monin_obukhov/monin_obukhov_kernel.F90 .
./lengthscaleeq.F90: $(SRCROOT)mom5/ocean_param/gotm-4.0/turbulence/lengthscaleeq.F90
	cp $(SRCROOT)mom5/ocean_param/gotm-4.0/turbulence/lengthscaleeq.F90 .
./mpp_global_sum_tl.h: /home/wenchang/tigress/CM2.5_tiger2_wy_killtc/exp/../src/shared/mpp/include/mpp_global_sum_tl.h
	cp /home/wenchang/tigress/CM2.5_tiger2_wy_killtc/exp/../src/shared/mpp/include/mpp_global_sum_tl.h .
./debug.inc: /home/wenchang/tigress/CM2.5_tiger2_wy_killtc/exp/../src/land_lad2/shared/debug.inc
	cp /home/wenchang/tigress/CM2.5_tiger2_wy_killtc/exp/../src/land_lad2/shared/debug.inc .
./mpp_chksum_int.h: /home/wenchang/tigress/CM2.5_tiger2_wy_killtc/exp/../src/shared/mpp/include/mpp_chksum_int.h
	cp /home/wenchang/tigress/CM2.5_tiger2_wy_killtc/exp/../src/shared/mpp/include/mpp_chksum_int.h .
./format_input.F90: $(SRCROOT)atmos_param/cosp/quickbeam/format_input.F90
	cp $(SRCROOT)atmos_param/cosp/quickbeam/format_input.F90 .
./drifters_compute_k.h: /home/wenchang/tigress/CM2.5_tiger2_wy_killtc/exp/../src/shared/drifters/drifters_compute_k.h
	cp /home/wenchang/tigress/CM2.5_tiger2_wy_killtc/exp/../src/shared/drifters/drifters_compute_k.h .
./test_mpp_domains.F90: $(SRCROOT)shared/mpp/test_mpp_domains.F90
	cp $(SRCROOT)shared/mpp/test_mpp_domains.F90 .
./alpha_mnb.F90: $(SRCROOT)mom5/ocean_param/gotm-4.0/turbulence/alpha_mnb.F90
	cp $(SRCROOT)mom5/ocean_param/gotm-4.0/turbulence/alpha_mnb.F90 .
./atmos_sea_salt.F90: $(SRCROOT)atmos_shared/tracer_driver/atmos_sea_salt.F90
	cp $(SRCROOT)atmos_shared/tracer_driver/atmos_sea_salt.F90 .
./vegetation.F90: $(SRCROOT)land_lad2/vegetation/vegetation.F90
	cp $(SRCROOT)land_lad2/vegetation/vegetation.F90 .
./mpp_reduce_sma.h: /home/wenchang/tigress/CM2.5_tiger2_wy_killtc/exp/../src/shared/mpp/include/mpp_reduce_sma.h
	cp /home/wenchang/tigress/CM2.5_tiger2_wy_killtc/exp/../src/shared/mpp/include/mpp_reduce_sma.h .
./lidar_simulator.F90: $(SRCROOT)atmos_param/cosp/actsim/lidar_simulator.F90
	cp $(SRCROOT)atmos_param/cosp/actsim/lidar_simulator.F90 .
./read_remap_cohort_data.inc: /home/wenchang/tigress/CM2.5_tiger2_wy_killtc/exp/../src/land_lad2/vegetation/read_remap_cohort_data.inc
	cp /home/wenchang/tigress/CM2.5_tiger2_wy_killtc/exp/../src/land_lad2/vegetation/read_remap_cohort_data.inc .
./ocean_pressure.F90: $(SRCROOT)mom5/ocean_core/ocean_pressure.F90
	cp $(SRCROOT)mom5/ocean_core/ocean_pressure.F90 .
./ice_dyn.F90: $(SRCROOT)ice_sis/ice_dyn.F90
	cp $(SRCROOT)ice_sis/ice_dyn.F90 .
./convert_fluxes.F90: $(SRCROOT)mom5/ocean_param/gotm-4.0/util/convert_fluxes.F90
	cp $(SRCROOT)mom5/ocean_param/gotm-4.0/util/convert_fluxes.F90 .
./mpp_do_update.h: /home/wenchang/tigress/CM2.5_tiger2_wy_killtc/exp/../src/shared/mpp/include/mpp_do_update.h
	cp /home/wenchang/tigress/CM2.5_tiger2_wy_killtc/exp/../src/shared/mpp/include/mpp_do_update.h .
./ocean_nphysics_new.F90: $(SRCROOT)mom5/ocean_param/neutral/ocean_nphysics_new.F90
	cp $(SRCROOT)mom5/ocean_param/neutral/ocean_nphysics_new.F90 .
./ocean_vert_tidal.F90: $(SRCROOT)mom5/ocean_param/vertical/ocean_vert_tidal.F90
	cp $(SRCROOT)mom5/ocean_param/vertical/ocean_vert_tidal.F90 .
./ocean_vert_gotm.F90: $(SRCROOT)mom5/ocean_param/vertical/ocean_vert_gotm.F90
	cp $(SRCROOT)mom5/ocean_param/vertical/ocean_vert_gotm.F90 .
./ice_grid.F90: $(SRCROOT)ice_sis/ice_grid.F90
	cp $(SRCROOT)ice_sis/ice_grid.F90 .
./mt95.f90: $(SRCROOT)atmos_param/clubb/CLUBB_core/mt95.f90
	cp $(SRCROOT)atmos_param/clubb/CLUBB_core/mt95.f90 .
./fft99.F90: $(SRCROOT)shared/fft/fft99.F90
	cp $(SRCROOT)shared/fft/fft99.F90 .
./strat_chem_utilities.F90: $(SRCROOT)atmos_shared/tracer_driver/tropchem/strat_chem_utilities.F90
	cp $(SRCROOT)atmos_shared/tracer_driver/tropchem/strat_chem_utilities.F90 .
./oda_core_ecda.F90: $(SRCROOT)shared/oda_tools/oda_core_ecda.F90
	cp $(SRCROOT)shared/oda_tools/oda_core_ecda.F90 .
./FMS_ocmip2_co2calc.F90: $(SRCROOT)ocean_shared/generic_tracers/FMS_ocmip2_co2calc.F90
	cp $(SRCROOT)ocean_shared/generic_tracers/FMS_ocmip2_co2calc.F90 .
./vegn_radiation.F90: $(SRCROOT)land_lad2/vegetation/vegn_radiation.F90
	cp $(SRCROOT)land_lad2/vegetation/vegn_radiation.F90 .
./ocean_nphysicsA.F90: $(SRCROOT)mom5/ocean_param/neutral/ocean_nphysicsA.F90
	cp $(SRCROOT)mom5/ocean_param/neutral/ocean_nphysicsA.F90 .
./ocean_nphysics_util_new.F90: $(SRCROOT)mom5/ocean_param/neutral/ocean_nphysics_util_new.F90
	cp $(SRCROOT)mom5/ocean_param/neutral/ocean_nphysics_util_new.F90 .
./shortwave_driver.F90: $(SRCROOT)atmos_param/sea_esf_rad/shortwave_driver.F90
	cp $(SRCROOT)atmos_param/sea_esf_rad/shortwave_driver.F90 .
./getput_compressed.inc: /home/wenchang/tigress/CM2.5_tiger2_wy_killtc/exp/../src/land_lad2/shared/nf_utils/getput_compressed.inc
	cp /home/wenchang/tigress/CM2.5_tiger2_wy_killtc/exp/../src/land_lad2/shared/nf_utils/getput_compressed.inc .
./mpp_domains_misc.inc: /home/wenchang/tigress/CM2.5_tiger2_wy_killtc/exp/../src/shared/mpp/include/mpp_domains_misc.inc
	cp /home/wenchang/tigress/CM2.5_tiger2_wy_killtc/exp/../src/shared/mpp/include/mpp_domains_misc.inc .
./flux_exchange.F90: $(SRCROOT)coupler/flux_exchange.F90
	cp $(SRCROOT)coupler/flux_exchange.F90 .
./cosp_types.F90: $(SRCROOT)atmos_param/cosp/cosp_types.F90
	cp $(SRCROOT)atmos_param/cosp/cosp_types.F90 .
./dsd.F90: $(SRCROOT)atmos_param/cosp/quickbeam/dsd.F90
	cp $(SRCROOT)atmos_param/cosp/quickbeam/dsd.F90 .
./ocean_thickness.F90: $(SRCROOT)mom5/ocean_core/ocean_thickness.F90
	cp $(SRCROOT)mom5/ocean_core/ocean_thickness.F90 .
./longwave.F90: $(SRCROOT)atmos_param/fsrad/longwave.F90
	cp $(SRCROOT)atmos_param/fsrad/longwave.F90 .
./fm_util.F90: $(SRCROOT)shared/field_manager/fm_util.F90
	cp $(SRCROOT)shared/field_manager/fm_util.F90 .
./mg_const.F90: $(SRCROOT)atmos_param/strat_cloud/mg_const.F90
	cp $(SRCROOT)atmos_param/strat_cloud/mg_const.F90 .
./parameters_radiation.F90: $(SRCROOT)atmos_param/clubb/CLUBB_core/parameters_radiation.F90
	cp $(SRCROOT)atmos_param/clubb/CLUBB_core/parameters_radiation.F90 .
./ensemble_manager.F90: $(SRCROOT)shared/coupler/ensemble_manager.F90
	cp $(SRCROOT)shared/coupler/ensemble_manager.F90 .
./monin_obukhov.F90: $(SRCROOT)atmos_param/monin_obukhov/monin_obukhov.F90
	cp $(SRCROOT)atmos_param/monin_obukhov/monin_obukhov.F90 .
./gotm_lib_version.F90: $(SRCROOT)mom5/ocean_param/gotm-4.0/turbulence/gotm_lib_version.F90
	cp $(SRCROOT)mom5/ocean_param/gotm-4.0/turbulence/gotm_lib_version.F90 .
./diag_integral.F90: $(SRCROOT)atmos_param/diag_integral/diag_integral.F90
	cp $(SRCROOT)atmos_param/diag_integral/diag_integral.F90 .
./dissipationeq.F90: $(SRCROOT)mom5/ocean_param/gotm-4.0/turbulence/dissipationeq.F90
	cp $(SRCROOT)mom5/ocean_param/gotm-4.0/turbulence/dissipationeq.F90 .
./my25_turb.F90: $(SRCROOT)atmos_param/my25_turb/my25_turb.F90
	cp $(SRCROOT)atmos_param/my25_turb/my25_turb.F90 .
./mcm_lw.F90: $(SRCROOT)atmos_param/fsrad/mcm_lw.F90
	cp $(SRCROOT)atmos_param/fsrad/mcm_lw.F90 .
./recl.inc: /home/wenchang/tigress/CM2.5_tiger2_wy_killtc/exp/../src/atmos_param/clubb/CLUBB_core/recl.inc
	cp /home/wenchang/tigress/CM2.5_tiger2_wy_killtc/exp/../src/atmos_param/clubb/CLUBB_core/recl.inc .
./ocean_submesoscale.F90: $(SRCROOT)mom5/ocean_param/lateral/ocean_submesoscale.F90
	cp $(SRCROOT)mom5/ocean_param/lateral/ocean_submesoscale.F90 .
./CLUBB_driver_SCM.F90: $(SRCROOT)atmos_param/clubb/CLUBB_driver_SCM.F90
	cp $(SRCROOT)atmos_param/clubb/CLUBB_driver_SCM.F90 .
./moist_conv.F90: $(SRCROOT)atmos_param/moist_conv/moist_conv.F90
	cp $(SRCROOT)atmos_param/moist_conv/moist_conv.F90 .
./input_reader.F90: $(SRCROOT)atmos_param/clubb/CLUBB_core/input_reader.F90
	cp $(SRCROOT)atmos_param/clubb/CLUBB_core/input_reader.F90 .
./land_tile_diag.F90: $(SRCROOT)land_lad2/shared/land_tile_diag.F90
	cp $(SRCROOT)land_lad2/shared/land_tile_diag.F90 .
./turbulence.F90: $(SRCROOT)mom5/ocean_param/gotm-4.0/turbulence/turbulence.F90
	cp $(SRCROOT)mom5/ocean_param/gotm-4.0/turbulence/turbulence.F90 .
./test_mpp_io.F90: $(SRCROOT)shared/mpp/test_mpp_io.F90
	cp $(SRCROOT)shared/mpp/test_mpp_io.F90 .
./ode_solvers.F90: $(SRCROOT)mom5/ocean_param/gotm-4.0/util/ode_solvers.F90
	cp $(SRCROOT)mom5/ocean_param/gotm-4.0/util/ode_solvers.F90 .
./diag_grid.F90: $(SRCROOT)shared/diag_manager/diag_grid.F90
	cp $(SRCROOT)shared/diag_manager/diag_grid.F90 .
./lhsw_driver.F90: $(SRCROOT)atmos_param/sea_esf_rad/lhsw_driver.F90
	cp $(SRCROOT)atmos_param/sea_esf_rad/lhsw_driver.F90 .
./advance_helper_module.F90: $(SRCROOT)atmos_param/clubb/CLUBB_core/advance_helper_module.F90
	cp $(SRCROOT)atmos_param/clubb/CLUBB_core/advance_helper_module.F90 .
./transitions.F90: $(SRCROOT)land_lad2/transitions/transitions.F90
	cp $(SRCROOT)land_lad2/transitions/transitions.F90 .
./ocean_tracer_diag.F90: $(SRCROOT)mom5/ocean_diag/ocean_tracer_diag.F90
	cp $(SRCROOT)mom5/ocean_diag/ocean_tracer_diag.F90 .
./ice_type.F90: $(SRCROOT)ice_sis/ice_type.F90
	cp $(SRCROOT)ice_sis/ice_type.F90 .
./cosp_constants.F90: $(SRCROOT)atmos_param/cosp/cosp_constants.F90
	cp $(SRCROOT)atmos_param/cosp/cosp_constants.F90 .
./physics_driver.F90: $(SRCROOT)atmos_param/physics_driver/physics_driver.F90
	cp $(SRCROOT)atmos_param/physics_driver/physics_driver.F90 .
./a2b_edge.F90: $(SRCROOT)atmos_cubed_sphere/model/a2b_edge.F90
	cp $(SRCROOT)atmos_cubed_sphere/model/a2b_edge.F90 .
./clubb_core.F90: $(SRCROOT)atmos_param/clubb/CLUBB_core/clubb_core.F90
	cp $(SRCROOT)atmos_param/clubb/CLUBB_core/clubb_core.F90 .
./vegn_cohort.F90: $(SRCROOT)land_lad2/vegetation/vegn_cohort.F90
	cp $(SRCROOT)land_lad2/vegetation/vegn_cohort.F90 .
./variables_prognostic_module.F90: $(SRCROOT)atmos_param/clubb/CLUBB_core/variables_prognostic_module.F90
	cp $(SRCROOT)atmos_param/clubb/CLUBB_core/variables_prognostic_module.F90 .
./mpp_do_global_field.h: /home/wenchang/tigress/CM2.5_tiger2_wy_killtc/exp/../src/shared/mpp/include/mpp_do_global_field.h
	cp /home/wenchang/tigress/CM2.5_tiger2_wy_killtc/exp/../src/shared/mpp/include/mpp_do_global_field.h .
./moz.subs.F90: $(SRCROOT)atmos_shared/tracer_driver/tropchem/moz.subs.F90
	cp $(SRCROOT)atmos_shared/tracer_driver/tropchem/moz.subs.F90 .
./constants.F90: $(SRCROOT)shared/constants/constants.F90
	cp $(SRCROOT)shared/constants/constants.F90 .
./ocean_nphysics_flux.F90: $(SRCROOT)mom5/ocean_param/neutral/ocean_nphysics_flux.F90
	cp $(SRCROOT)mom5/ocean_param/neutral/ocean_nphysics_flux.F90 .
./ocean_ibgc.F90: $(SRCROOT)mom5/ocean_bgc/ocean_ibgc.F90
	cp $(SRCROOT)mom5/ocean_bgc/ocean_ibgc.F90 .
./variances.F90: $(SRCROOT)mom5/ocean_param/gotm-4.0/turbulence/variances.F90
	cp $(SRCROOT)mom5/ocean_param/gotm-4.0/turbulence/variances.F90 .
./mo_photo.F90: $(SRCROOT)atmos_shared/tracer_driver/tropchem/mo_photo.F90
	cp $(SRCROOT)atmos_shared/tracer_driver/tropchem/mo_photo.F90 .
./microphys_cloud.F90: $(SRCROOT)atmos_param/sea_esf_rad/microphys_cloud.F90
	cp $(SRCROOT)atmos_param/sea_esf_rad/microphys_cloud.F90 .
./land_tile_diag_sel.F90: $(SRCROOT)land_lad2/shared/land_tile_diag_sel.F90
	cp $(SRCROOT)land_lad2/shared/land_tile_diag_sel.F90 .
./ocean_bihcgrid_friction.F90: $(SRCROOT)mom5/ocean_param/lateral/ocean_bihcgrid_friction.F90
	cp $(SRCROOT)mom5/ocean_param/lateral/ocean_bihcgrid_friction.F90 .
./mpp_do_check.h: /home/wenchang/tigress/CM2.5_tiger2_wy_killtc/exp/../src/shared/mpp/include/mpp_do_check.h
	cp /home/wenchang/tigress/CM2.5_tiger2_wy_killtc/exp/../src/shared/mpp/include/mpp_do_check.h .
./read_mosaic.c: $(SRCROOT)shared/mosaic/read_mosaic.c
	cp $(SRCROOT)shared/mosaic/read_mosaic.c .
./ocean_overexchange.F90: $(SRCROOT)mom5/ocean_param/sources/ocean_overexchange.F90
	cp $(SRCROOT)mom5/ocean_param/sources/ocean_overexchange.F90 .
./interp.c: $(SRCROOT)shared/mosaic/interp.c
	cp $(SRCROOT)shared/mosaic/interp.c .
./mpp_error_a_a.h: /home/wenchang/tigress/CM2.5_tiger2_wy_killtc/exp/../src/shared/mpp/include/mpp_error_a_a.h
	cp /home/wenchang/tigress/CM2.5_tiger2_wy_killtc/exp/../src/shared/mpp/include/mpp_error_a_a.h .
./hconst.F90: $(SRCROOT)atmos_param/fsrad/hconst.F90
	cp $(SRCROOT)atmos_param/fsrad/hconst.F90 .
./genericeq.F90: $(SRCROOT)mom5/ocean_param/gotm-4.0/turbulence/genericeq.F90
	cp $(SRCROOT)mom5/ocean_param/gotm-4.0/turbulence/genericeq.F90 .
./cosp_misr_simulator.F90: $(SRCROOT)atmos_param/cosp/cosp_misr_simulator.F90
	cp $(SRCROOT)atmos_param/cosp/cosp_misr_simulator.F90 .
./r_ratio.F90: $(SRCROOT)mom5/ocean_param/gotm-4.0/turbulence/r_ratio.F90
	cp $(SRCROOT)mom5/ocean_param/gotm-4.0/turbulence/r_ratio.F90 .
./cosp_stats.F90: $(SRCROOT)atmos_param/cosp/cosp_stats.F90
	cp $(SRCROOT)atmos_param/cosp/cosp_stats.F90 .
./surface_flux.F90: $(SRCROOT)coupler/surface_flux.F90
	cp $(SRCROOT)coupler/surface_flux.F90 .
./ocean_vert_kpp_test.F90: $(SRCROOT)mom5/ocean_param/vertical/ocean_vert_kpp_test.F90
	cp $(SRCROOT)mom5/ocean_param/vertical/ocean_vert_kpp_test.F90 .
./ocean_sbc.F90: $(SRCROOT)mom5/ocean_core/ocean_sbc.F90
	cp $(SRCROOT)mom5/ocean_core/ocean_sbc.F90 .
./conv_utilities_k.F90: $(SRCROOT)atmos_param/shallow_cu/conv_utilities_k.F90
	cp $(SRCROOT)atmos_param/shallow_cu/conv_utilities_k.F90 .
./generic_miniBLING.F90: $(SRCROOT)ocean_shared/generic_tracers/generic_miniBLING.F90
	cp $(SRCROOT)ocean_shared/generic_tracers/generic_miniBLING.F90 .
./mpp_get_boundary.h: /home/wenchang/tigress/CM2.5_tiger2_wy_killtc/exp/../src/shared/mpp/include/mpp_get_boundary.h
	cp /home/wenchang/tigress/CM2.5_tiger2_wy_killtc/exp/../src/shared/mpp/include/mpp_get_boundary.h .
./mpp_do_update_nest.h: /home/wenchang/tigress/CM2.5_tiger2_wy_killtc/exp/../src/shared/mpp/include/mpp_do_update_nest.h
	cp /home/wenchang/tigress/CM2.5_tiger2_wy_killtc/exp/../src/shared/mpp/include/mpp_do_update_nest.h .
./fv_control.F90: $(SRCROOT)atmos_cubed_sphere/model/fv_control.F90
	cp $(SRCROOT)atmos_cubed_sphere/model/fv_control.F90 .
./drifters.F90: $(SRCROOT)shared/drifters/drifters.F90
	cp $(SRCROOT)shared/drifters/drifters.F90 .
./constants_clubb.F90: $(SRCROOT)atmos_param/clubb/CLUBB_core/constants_clubb.F90
	cp $(SRCROOT)atmos_param/clubb/CLUBB_core/constants_clubb.F90 .
./mcm_sw_driver.F90: $(SRCROOT)atmos_param/fsrad/mcm_sw_driver.F90
	cp $(SRCROOT)atmos_param/fsrad/mcm_sw_driver.F90 .
./grid.F90: $(SRCROOT)shared/mosaic/grid.F90
	cp $(SRCROOT)shared/mosaic/grid.F90 .
./land_tile_diag_buff.F90: $(SRCROOT)land_lad2/shared/land_tile_diag_buff.F90
	cp $(SRCROOT)land_lad2/shared/land_tile_diag_buff.F90 .
./coupler_types.F90: $(SRCROOT)shared/coupler/coupler_types.F90
	cp $(SRCROOT)shared/coupler/coupler_types.F90 .
./moist_processes_utils.F90: $(SRCROOT)atmos_param/moist_processes/moist_processes_utils.F90
	cp $(SRCROOT)atmos_param/moist_processes/moist_processes_utils.F90 .
./snow_tile.F90: $(SRCROOT)land_lad2/snow/snow_tile.F90
	cp $(SRCROOT)land_lad2/snow/snow_tile.F90 .
./mo_chemini.F90: $(SRCROOT)atmos_shared/tracer_driver/tropchem/mo_chemini.F90
	cp $(SRCROOT)atmos_shared/tracer_driver/tropchem/mo_chemini.F90 .
./fv_mp_mod.F90: $(SRCROOT)atmos_cubed_sphere/tools/fv_mp_mod.F90
	cp $(SRCROOT)atmos_cubed_sphere/tools/fv_mp_mod.F90 .
./time_interp_external.F90: $(SRCROOT)shared/time_interp/time_interp_external.F90
	cp $(SRCROOT)shared/time_interp/time_interp_external.F90 .
./mpp_domains_reduce.inc: /home/wenchang/tigress/CM2.5_tiger2_wy_killtc/exp/../src/shared/mpp/include/mpp_domains_reduce.inc
	cp /home/wenchang/tigress/CM2.5_tiger2_wy_killtc/exp/../src/shared/mpp/include/mpp_domains_reduce.inc .
./mpp_error_s_s.h: /home/wenchang/tigress/CM2.5_tiger2_wy_killtc/exp/../src/shared/mpp/include/mpp_error_s_s.h
	cp /home/wenchang/tigress/CM2.5_tiger2_wy_killtc/exp/../src/shared/mpp/include/mpp_error_s_s.h .
./ocean_increment_velocity.F90: $(SRCROOT)mom5/ocean_param/sources/ocean_increment_velocity.F90
	cp $(SRCROOT)mom5/ocean_param/sources/ocean_increment_velocity.F90 .
./vegn_disturbance.F90: $(SRCROOT)land_lad2/vegetation/vegn_disturbance.F90
	cp $(SRCROOT)land_lad2/vegetation/vegn_disturbance.F90 .
./strat_netcdf.F90: $(SRCROOT)atmos_param/strat_cloud/strat_netcdf.F90
	cp $(SRCROOT)atmos_param/strat_cloud/strat_netcdf.F90 .
./mpp_do_get_boundary.h: /home/wenchang/tigress/CM2.5_tiger2_wy_killtc/exp/../src/shared/mpp/include/mpp_do_get_boundary.h
	cp /home/wenchang/tigress/CM2.5_tiger2_wy_killtc/exp/../src/shared/mpp/include/mpp_do_get_boundary.h .
./ocean_pert_co2.F90: $(SRCROOT)mom5/ocean_bgc/ocean_pert_co2.F90
	cp $(SRCROOT)mom5/ocean_bgc/ocean_pert_co2.F90 .
./mono_flux_limiter.F90: $(SRCROOT)atmos_param/clubb/CLUBB_core/mono_flux_limiter.F90
	cp $(SRCROOT)atmos_param/clubb/CLUBB_core/mono_flux_limiter.F90 .
./ocean_po4_pre.F90: $(SRCROOT)mom5/ocean_bgc/ocean_po4_pre.F90
	cp $(SRCROOT)mom5/ocean_bgc/ocean_po4_pre.F90 .
./fv_grid_utils.F90: $(SRCROOT)atmos_cubed_sphere/model/fv_grid_utils.F90
	cp $(SRCROOT)atmos_cubed_sphere/model/fv_grid_utils.F90 .
./mpp_read_compressed.h: /home/wenchang/tigress/CM2.5_tiger2_wy_killtc/exp/../src/shared/mpp/include/mpp_read_compressed.h
	cp /home/wenchang/tigress/CM2.5_tiger2_wy_killtc/exp/../src/shared/mpp/include/mpp_read_compressed.h .
./ocean_bih_tracer.F90: $(SRCROOT)mom5/ocean_param/lateral/ocean_bih_tracer.F90
	cp $(SRCROOT)mom5/ocean_param/lateral/ocean_bih_tracer.F90 .
./fv_update_phys.F90: $(SRCROOT)atmos_cubed_sphere/model/fv_update_phys.F90
	cp $(SRCROOT)atmos_cubed_sphere/model/fv_update_phys.F90 .
./entrain.F90: $(SRCROOT)atmos_param/entrain/entrain.F90
	cp $(SRCROOT)atmos_param/entrain/entrain.F90 .
./test_horiz_interp.F90: $(SRCROOT)shared/horiz_interp/test_horiz_interp.F90
	cp $(SRCROOT)shared/horiz_interp/test_horiz_interp.F90 .
./fms_donner.F90: $(SRCROOT)atmos_param/donner_deep/fms_donner.F90
	cp $(SRCROOT)atmos_param/donner_deep/fms_donner.F90 .
./MG_microp_3D.F90: $(SRCROOT)atmos_param/clubb/MG_microp_3D.F90
	cp $(SRCROOT)atmos_param/clubb/MG_microp_3D.F90 .
./ocean_residency_meta.F90: $(SRCROOT)mom5/ocean_tracers/ocean_residency_meta.F90
	cp $(SRCROOT)mom5/ocean_tracers/ocean_residency_meta.F90 .
./mpp_comm_sma.inc: /home/wenchang/tigress/CM2.5_tiger2_wy_killtc/exp/../src/shared/mpp/include/mpp_comm_sma.inc
	cp /home/wenchang/tigress/CM2.5_tiger2_wy_killtc/exp/../src/shared/mpp/include/mpp_comm_sma.inc .
./fv_physics.F90: $(SRCROOT)atmos_cubed_sphere/driver/coupled/fv_physics.F90
	cp $(SRCROOT)atmos_cubed_sphere/driver/coupled/fv_physics.F90 .
./ocean_vert_mix.F90: $(SRCROOT)mom5/ocean_param/vertical/ocean_vert_mix.F90
	cp $(SRCROOT)mom5/ocean_param/vertical/ocean_vert_mix.F90 .
./mpp_util_sma.inc: /home/wenchang/tigress/CM2.5_tiger2_wy_killtc/exp/../src/shared/mpp/include/mpp_util_sma.inc
	cp /home/wenchang/tigress/CM2.5_tiger2_wy_killtc/exp/../src/shared/mpp/include/mpp_util_sma.inc .
./clubb_precision.F90: $(SRCROOT)atmos_param/clubb/CLUBB_core/clubb_precision.F90
	cp $(SRCROOT)atmos_param/clubb/CLUBB_core/clubb_precision.F90 .
./ocean_wave.F90: $(SRCROOT)mom5/ocean_wave/ocean_wave.F90
	cp $(SRCROOT)mom5/ocean_wave/ocean_wave.F90 .
./gas_tf.F90: $(SRCROOT)atmos_param/sea_esf_rad/gas_tf.F90
	cp $(SRCROOT)atmos_param/sea_esf_rad/gas_tf.F90 .
./table_printer.F90: $(SRCROOT)land_lad2/shared/table_printer.F90
	cp $(SRCROOT)land_lad2/shared/table_printer.F90 .
./fv_io.F90: $(SRCROOT)atmos_cubed_sphere/tools/fv_io.F90
	cp $(SRCROOT)atmos_cubed_sphere/tools/fv_io.F90 .
./prec_scops.F: $(SRCROOT)atmos_param/cosp/llnl/prec_scops.F
	cp $(SRCROOT)atmos_param/cosp/llnl/prec_scops.F .
./platform.F90: $(SRCROOT)shared/platform/platform.F90
	cp $(SRCROOT)shared/platform/platform.F90 .
./moz.mods.F90: $(SRCROOT)atmos_shared/tracer_driver/tropchem/moz.mods.F90
	cp $(SRCROOT)atmos_shared/tracer_driver/tropchem/moz.mods.F90 .
./land_utils.F90: $(SRCROOT)land_lad2/shared/land_utils.F90
	cp $(SRCROOT)land_lad2/shared/land_utils.F90 .
./endian.F90: $(SRCROOT)atmos_param/clubb/CLUBB_core/endian.F90
	cp $(SRCROOT)atmos_param/clubb/CLUBB_core/endian.F90 .
./ocmip2_co2calc.F90: $(SRCROOT)mom5/ocean_bgc/ocmip2_co2calc.F90
	cp $(SRCROOT)mom5/ocean_bgc/ocmip2_co2calc.F90 .
./vert_diff.F90: $(SRCROOT)atmos_param/vert_diff/vert_diff.F90
	cp $(SRCROOT)atmos_param/vert_diff/vert_diff.F90 .
./mpp_data_mpi.inc: /home/wenchang/tigress/CM2.5_tiger2_wy_killtc/exp/../src/shared/mpp/include/mpp_data_mpi.inc
	cp /home/wenchang/tigress/CM2.5_tiger2_wy_killtc/exp/../src/shared/mpp/include/mpp_data_mpi.inc .
./mpp_transmit_nocomm.h: /home/wenchang/tigress/CM2.5_tiger2_wy_killtc/exp/../src/shared/mpp/include/mpp_transmit_nocomm.h
	cp /home/wenchang/tigress/CM2.5_tiger2_wy_killtc/exp/../src/shared/mpp/include/mpp_transmit_nocomm.h .
./sponge_layer_damping.F90: $(SRCROOT)atmos_param/clubb/CLUBB_core/sponge_layer_damping.F90
	cp $(SRCROOT)atmos_param/clubb/CLUBB_core/sponge_layer_damping.F90 .
./fft.F90: $(SRCROOT)shared/fft/fft.F90
	cp $(SRCROOT)shared/fft/fft.F90 .
./pdf_closure_module.F90: $(SRCROOT)atmos_param/clubb/CLUBB_core/pdf_closure_module.F90
	cp $(SRCROOT)atmos_param/clubb/CLUBB_core/pdf_closure_module.F90 .
./mg_drag.F90: $(SRCROOT)atmos_param/mg_drag/mg_drag.F90
	cp $(SRCROOT)atmos_param/mg_drag/mg_drag.F90 .
./kbeq.F90: $(SRCROOT)mom5/ocean_param/gotm-4.0/turbulence/kbeq.F90
	cp $(SRCROOT)mom5/ocean_param/gotm-4.0/turbulence/kbeq.F90 .
./atmos_ch3i.F90: $(SRCROOT)atmos_shared/tracer_driver/atmos_ch3i.F90
	cp $(SRCROOT)atmos_shared/tracer_driver/atmos_ch3i.F90 .
./parse.inc: /home/wenchang/tigress/CM2.5_tiger2_wy_killtc/exp/../src/shared/field_manager/parse.inc
	cp /home/wenchang/tigress/CM2.5_tiger2_wy_killtc/exp/../src/shared/field_manager/parse.inc .
./vert_turb_driver.F90: $(SRCROOT)atmos_param/vert_turb_driver/vert_turb_driver.F90
	cp $(SRCROOT)atmos_param/vert_turb_driver/vert_turb_driver.F90 .
./nonfms_donner.F90: $(SRCROOT)atmos_param/donner_deep/nonfms_donner.F90
	cp $(SRCROOT)atmos_param/donner_deep/nonfms_donner.F90 .
./mpp_io_misc.inc: /home/wenchang/tigress/CM2.5_tiger2_wy_killtc/exp/../src/shared/mpp/include/mpp_io_misc.inc
	cp /home/wenchang/tigress/CM2.5_tiger2_wy_killtc/exp/../src/shared/mpp/include/mpp_io_misc.inc .
./mpp_global_field.h: /home/wenchang/tigress/CM2.5_tiger2_wy_killtc/exp/../src/shared/mpp/include/mpp_global_field.h
	cp /home/wenchang/tigress/CM2.5_tiger2_wy_killtc/exp/../src/shared/mpp/include/mpp_global_field.h .
./gamma_mg.F90: $(SRCROOT)atmos_param/strat_cloud/gamma_mg.F90
	cp $(SRCROOT)atmos_param/strat_cloud/gamma_mg.F90 .
./fv_arrays.F90: $(SRCROOT)atmos_cubed_sphere/model/fv_arrays.F90
	cp $(SRCROOT)atmos_cubed_sphere/model/fv_arrays.F90 .
./stat_file_module.F90: $(SRCROOT)atmos_param/clubb/CLUBB_core/stat_file_module.F90
	cp $(SRCROOT)atmos_param/clubb/CLUBB_core/stat_file_module.F90 .
./oda_types.F90: $(SRCROOT)shared/oda_tools/oda_types.F90
	cp $(SRCROOT)shared/oda_tools/oda_types.F90 .
./generic_COBALT.F90: $(SRCROOT)ocean_shared/generic_tracers/generic_COBALT.F90
	cp $(SRCROOT)ocean_shared/generic_tracers/generic_COBALT.F90 .
./cosp_lidar.F90: $(SRCROOT)atmos_param/cosp/cosp_lidar.F90
	cp $(SRCROOT)atmos_param/cosp/cosp_lidar.F90 .
./atmos_radon.F90: $(SRCROOT)atmos_shared/tracer_driver/atmos_radon.F90
	cp $(SRCROOT)atmos_shared/tracer_driver/atmos_radon.F90 .
./gmres_wrap.F90: $(SRCROOT)atmos_param/clubb/CLUBB_core/gmres_wrap.F90
	cp $(SRCROOT)atmos_param/clubb/CLUBB_core/gmres_wrap.F90 .
./vegn_photosynthesis.F90: $(SRCROOT)land_lad2/vegetation/vegn_photosynthesis.F90
	cp $(SRCROOT)land_lad2/vegetation/vegn_photosynthesis.F90 .
./output_netcdf.F90: $(SRCROOT)atmos_param/clubb/CLUBB_core/output_netcdf.F90
	cp $(SRCROOT)atmos_param/clubb/CLUBB_core/output_netcdf.F90 .
./aerosolrad_package.F90: $(SRCROOT)atmos_param/sea_esf_rad/aerosolrad_package.F90
	cp $(SRCROOT)atmos_param/sea_esf_rad/aerosolrad_package.F90 .
./fms_platform.h: /home/wenchang/tigress/CM2.5_tiger2_wy_killtc/exp/../src/shared/include/fms_platform.h
	cp /home/wenchang/tigress/CM2.5_tiger2_wy_killtc/exp/../src/shared/include/fms_platform.h .
./rh_based_clouds.F90: $(SRCROOT)atmos_param/sea_esf_rad/rh_based_clouds.F90
	cp $(SRCROOT)atmos_param/sea_esf_rad/rh_based_clouds.F90 .
./ocean_obc.F90: $(SRCROOT)mom5/ocean_core/ocean_obc.F90
	cp $(SRCROOT)mom5/ocean_core/ocean_obc.F90 .
./mo_chem_utls.F90: $(SRCROOT)atmos_shared/tracer_driver/tropchem/mo_chem_utls.F90
	cp $(SRCROOT)atmos_shared/tracer_driver/tropchem/mo_chem_utls.F90 .
./horiz_interp_type.F90: $(SRCROOT)shared/horiz_interp/horiz_interp_type.F90
	cp $(SRCROOT)shared/horiz_interp/horiz_interp_type.F90 .
./stable_bl_turb.F90: $(SRCROOT)atmos_param/stable_bl_turb/stable_bl_turb.F90
	cp $(SRCROOT)atmos_param/stable_bl_turb/stable_bl_turb.F90 .
./mpp_sum_nocomm.h: /home/wenchang/tigress/CM2.5_tiger2_wy_killtc/exp/../src/shared/mpp/include/mpp_sum_nocomm.h
	cp /home/wenchang/tigress/CM2.5_tiger2_wy_killtc/exp/../src/shared/mpp/include/mpp_sum_nocomm.h .
./river.F90: $(SRCROOT)land_lad2/river/river.F90
	cp $(SRCROOT)land_lad2/river/river.F90 .
./ocean_velocity_diag.F90: $(SRCROOT)mom5/ocean_diag/ocean_velocity_diag.F90
	cp $(SRCROOT)mom5/ocean_diag/ocean_velocity_diag.F90 .
./ocean_residency_integrand.F90: $(SRCROOT)mom5/ocean_tracers/ocean_residency_integrand.F90
	cp $(SRCROOT)mom5/ocean_tracers/ocean_residency_integrand.F90 .
./ocean_mixdownslope.F90: $(SRCROOT)mom5/ocean_param/lateral/ocean_mixdownslope.F90
	cp $(SRCROOT)mom5/ocean_param/lateral/ocean_mixdownslope.F90 .
./ice_model.F90: $(SRCROOT)ice_sis/ice_model.F90
	cp $(SRCROOT)ice_sis/ice_model.F90 .
./parameter_indices.F90: $(SRCROOT)atmos_param/clubb/CLUBB_core/parameter_indices.F90
	cp $(SRCROOT)atmos_param/clubb/CLUBB_core/parameter_indices.F90 .
./ocean_tracer.F90: $(SRCROOT)mom5/ocean_tracers/ocean_tracer.F90
	cp $(SRCROOT)mom5/ocean_tracers/ocean_tracer.F90 .
./mpp_utilities.F90: $(SRCROOT)shared/mpp/mpp_utilities.F90
	cp $(SRCROOT)shared/mpp/mpp_utilities.F90 .
./mpp_do_updateV_nonblock.h: /home/wenchang/tigress/CM2.5_tiger2_wy_killtc/exp/../src/shared/mpp/include/mpp_do_updateV_nonblock.h
	cp /home/wenchang/tigress/CM2.5_tiger2_wy_killtc/exp/../src/shared/mpp/include/mpp_do_updateV_nonblock.h .
./land_data.F90: $(SRCROOT)land_lad2/land_data.F90
	cp $(SRCROOT)land_lad2/land_data.F90 .
./ocean_sigma_transport.F90: $(SRCROOT)mom5/ocean_param/lateral/ocean_sigma_transport.F90
	cp $(SRCROOT)mom5/ocean_param/lateral/ocean_sigma_transport.F90 .
./cosp_diagnostics.F90: $(SRCROOT)atmos_param/cosp/cosp_diagnostics.F90
	cp $(SRCROOT)atmos_param/cosp/cosp_diagnostics.F90 .
./data_override.F90: $(SRCROOT)shared/data_override/data_override.F90
	cp $(SRCROOT)shared/data_override/data_override.F90 .
./ras.F90: $(SRCROOT)atmos_param/ras/ras.F90
	cp $(SRCROOT)atmos_param/ras/ras.F90 .
./mpp_global_sum.h: /home/wenchang/tigress/CM2.5_tiger2_wy_killtc/exp/../src/shared/mpp/include/mpp_global_sum.h
	cp /home/wenchang/tigress/CM2.5_tiger2_wy_killtc/exp/../src/shared/mpp/include/mpp_global_sum.h .
./time_manager.F90: $(SRCROOT)shared/time_manager/time_manager.F90
	cp $(SRCROOT)shared/time_manager/time_manager.F90 .
./soil.F90: $(SRCROOT)land_lad2/soil/soil.F90
	cp $(SRCROOT)land_lad2/soil/soil.F90 .
./isccp_clouds.F90: $(SRCROOT)atmos_param/sea_esf_rad/isccp_clouds.F90
	cp $(SRCROOT)atmos_param/sea_esf_rad/isccp_clouds.F90 .
./quicksort.F90: $(SRCROOT)shared/drifters/quicksort.F90
	cp $(SRCROOT)shared/drifters/quicksort.F90 .
./microphys_rad.F90: $(SRCROOT)atmos_param/sea_esf_rad/microphys_rad.F90
	cp $(SRCROOT)atmos_param/sea_esf_rad/microphys_rad.F90 .
./stats_LH_zt.F90: $(SRCROOT)atmos_param/clubb/CLUBB_core/stats_LH_zt.F90
	cp $(SRCROOT)atmos_param/clubb/CLUBB_core/stats_LH_zt.F90 .
./ocmip2_he.F90: $(SRCROOT)mom5/ocean_bgc/ocmip2_he.F90
	cp $(SRCROOT)mom5/ocean_bgc/ocmip2_he.F90 .
./diff_center.F90: $(SRCROOT)mom5/ocean_param/gotm-4.0/util/diff_center.F90
	cp $(SRCROOT)mom5/ocean_param/gotm-4.0/util/diff_center.F90 .
./ocean_momentum_source.F90: $(SRCROOT)mom5/ocean_param/sources/ocean_momentum_source.F90
	cp $(SRCROOT)mom5/ocean_param/sources/ocean_momentum_source.F90 .
./nc_cond.F90: $(SRCROOT)atmos_param/strat_cloud/nc_cond.F90
	cp $(SRCROOT)atmos_param/strat_cloud/nc_cond.F90 .
./cppdefs.h: /home/wenchang/tigress/CM2.5_tiger2_wy_killtc/exp/../src/mom5/ocean_param/gotm-4.0/include/cppdefs.h
	cp /home/wenchang/tigress/CM2.5_tiger2_wy_killtc/exp/../src/mom5/ocean_param/gotm-4.0/include/cppdefs.h .
./dyn_core.F90: $(SRCROOT)atmos_cubed_sphere/model/dyn_core.F90
	cp $(SRCROOT)atmos_cubed_sphere/model/dyn_core.F90 .
./atmos_co2.F90: $(SRCROOT)atmos_shared/tracer_driver/atmos_co2.F90
	cp $(SRCROOT)atmos_shared/tracer_driver/atmos_co2.F90 .
./read_mosaic.h: /home/wenchang/tigress/CM2.5_tiger2_wy_killtc/exp/../src/shared/mosaic/read_mosaic.h
	cp /home/wenchang/tigress/CM2.5_tiger2_wy_killtc/exp/../src/shared/mosaic/read_mosaic.h .
./ocean_sponges_eta.F90: $(SRCROOT)mom5/ocean_param/sources/ocean_sponges_eta.F90
	cp $(SRCROOT)mom5/ocean_param/sources/ocean_sponges_eta.F90 .
./ocean_obc_barotrop.F90: $(SRCROOT)mom5/ocean_core/ocean_obc_barotrop.F90
	cp $(SRCROOT)mom5/ocean_core/ocean_obc_barotrop.F90 .
./rh_clouds.F90: $(SRCROOT)atmos_param/rh_clouds/rh_clouds.F90
	cp $(SRCROOT)atmos_param/rh_clouds/rh_clouds.F90 .
./diag_util.F90: $(SRCROOT)shared/diag_manager/diag_util.F90
	cp $(SRCROOT)shared/diag_manager/diag_util.F90 .
./zeff.F90: $(SRCROOT)atmos_param/cosp/quickbeam/zeff.F90
	cp $(SRCROOT)atmos_param/cosp/quickbeam/zeff.F90 .
./column_diagnostics.F90: $(SRCROOT)shared/column_diagnostics/column_diagnostics.F90
	cp $(SRCROOT)shared/column_diagnostics/column_diagnostics.F90 .
./MISR_simulator.F: $(SRCROOT)atmos_param/cosp/MISR_simulator/MISR_simulator.F
	cp $(SRCROOT)atmos_param/cosp/MISR_simulator/MISR_simulator.F .
./mcm_swnew.F90: $(SRCROOT)atmos_param/fsrad/mcm_swnew.F90
	cp $(SRCROOT)atmos_param/fsrad/mcm_swnew.F90 .
./fv_climate_nudge.F90: $(SRCROOT)atmos_cubed_sphere/tools/fv_climate_nudge.F90
	cp $(SRCROOT)atmos_cubed_sphere/tools/fv_climate_nudge.F90 .
./atmos_carbon_aerosol.F90: $(SRCROOT)atmos_shared/tracer_driver/atmos_carbon_aerosol.F90
	cp $(SRCROOT)atmos_shared/tracer_driver/atmos_carbon_aerosol.F90 .
./mpp_util_nocomm.inc: /home/wenchang/tigress/CM2.5_tiger2_wy_killtc/exp/../src/shared/mpp/include/mpp_util_nocomm.inc
	cp /home/wenchang/tigress/CM2.5_tiger2_wy_killtc/exp/../src/shared/mpp/include/mpp_util_nocomm.inc .
./land_tile_io.F90: $(SRCROOT)land_lad2/shared/land_tile_io.F90
	cp $(SRCROOT)land_lad2/shared/land_tile_io.F90 .
./T_in_K_module.F90: $(SRCROOT)atmos_param/clubb/CLUBB_core/T_in_K_module.F90
	cp $(SRCROOT)atmos_param/clubb/CLUBB_core/T_in_K_module.F90 .
./field_manager.F90: $(SRCROOT)shared/field_manager/field_manager.F90
	cp $(SRCROOT)shared/field_manager/field_manager.F90 .
./longwave_params.F90: $(SRCROOT)atmos_param/sea_esf_rad/longwave_params.F90
	cp $(SRCROOT)atmos_param/sea_esf_rad/longwave_params.F90 .
./input_names.F90: $(SRCROOT)atmos_param/clubb/CLUBB_core/input_names.F90
	cp $(SRCROOT)atmos_param/clubb/CLUBB_core/input_names.F90 .
./sim_nc_mod.F90: $(SRCROOT)atmos_cubed_sphere/tools/sim_nc_mod.F90
	cp $(SRCROOT)atmos_cubed_sphere/tools/sim_nc_mod.F90 .
./sea_esf_rad.F90: $(SRCROOT)atmos_param/sea_esf_rad/sea_esf_rad.F90
	cp $(SRCROOT)atmos_param/sea_esf_rad/sea_esf_rad.F90 .
./land_constants.F90: $(SRCROOT)land_lad2/land_constants.F90
	cp $(SRCROOT)land_lad2/land_constants.F90 .
./cmue_d.F90: $(SRCROOT)mom5/ocean_param/gotm-4.0/turbulence/cmue_d.F90
	cp $(SRCROOT)mom5/ocean_param/gotm-4.0/turbulence/cmue_d.F90 .
./mpp.F90: $(SRCROOT)shared/mpp/mpp.F90
	cp $(SRCROOT)shared/mpp/mpp.F90 .
./tridiagonal.F90: $(SRCROOT)shared/tridiagonal/tridiagonal.F90
	cp $(SRCROOT)shared/tridiagonal/tridiagonal.F90 .
./optical_path.F90: $(SRCROOT)atmos_param/sea_esf_rad/optical_path.F90
	cp $(SRCROOT)atmos_param/sea_esf_rad/optical_path.F90 .
./ocean_tempsalt.F90: $(SRCROOT)mom5/ocean_tracers/ocean_tempsalt.F90
	cp $(SRCROOT)mom5/ocean_tracers/ocean_tempsalt.F90 .
./advance_windm_edsclrm_module.F90: $(SRCROOT)atmos_param/clubb/CLUBB_core/advance_windm_edsclrm_module.F90
	cp $(SRCROOT)atmos_param/clubb/CLUBB_core/advance_windm_edsclrm_module.F90 .
./file_functions.F90: $(SRCROOT)atmos_param/clubb/CLUBB_core/file_functions.F90
	cp $(SRCROOT)atmos_param/clubb/CLUBB_core/file_functions.F90 .
./ocean_lap_tracer.F90: $(SRCROOT)mom5/ocean_param/lateral/ocean_lap_tracer.F90
	cp $(SRCROOT)mom5/ocean_param/lateral/ocean_lap_tracer.F90 .
./cmue_rf.F90: $(SRCROOT)mom5/ocean_param/gotm-4.0/turbulence/cmue_rf.F90
	cp $(SRCROOT)mom5/ocean_param/gotm-4.0/turbulence/cmue_rf.F90 .
./ocean_residency.F90: $(SRCROOT)mom5/ocean_tracers/ocean_residency.F90
	cp $(SRCROOT)mom5/ocean_tracers/ocean_residency.F90 .
./variables_diagnostic_module.F90: $(SRCROOT)atmos_param/clubb/CLUBB_core/variables_diagnostic_module.F90
	cp $(SRCROOT)atmos_param/clubb/CLUBB_core/variables_diagnostic_module.F90 .
./mpp_io_write.inc: /home/wenchang/tigress/CM2.5_tiger2_wy_killtc/exp/../src/shared/mpp/include/mpp_io_write.inc
	cp /home/wenchang/tigress/CM2.5_tiger2_wy_killtc/exp/../src/shared/mpp/include/mpp_io_write.inc .
./read_data_4d.inc: /home/wenchang/tigress/CM2.5_tiger2_wy_killtc/exp/../src/shared/fms/read_data_4d.inc
	cp /home/wenchang/tigress/CM2.5_tiger2_wy_killtc/exp/../src/shared/fms/read_data_4d.inc .
./ocean_residency_range.F90: $(SRCROOT)mom5/ocean_tracers/ocean_residency_range.F90
	cp $(SRCROOT)mom5/ocean_tracers/ocean_residency_range.F90 .
./ocean_blob_dynamic_free.F90: $(SRCROOT)mom5/ocean_blobs/ocean_blob_dynamic_free.F90
	cp $(SRCROOT)mom5/ocean_blobs/ocean_blob_dynamic_free.F90 .
./read_climate_nudge_data.F90: $(SRCROOT)atmos_cubed_sphere/tools/read_climate_nudge_data.F90
	cp $(SRCROOT)atmos_cubed_sphere/tools/read_climate_nudge_data.F90 .
./cloud_generator.F90: $(SRCROOT)atmos_param/cloud_generator/cloud_generator.F90
	cp $(SRCROOT)atmos_param/cloud_generator/cloud_generator.F90 .
./sat_vapor_pres_k.F90: $(SRCROOT)shared/sat_vapor_pres/sat_vapor_pres_k.F90
	cp $(SRCROOT)shared/sat_vapor_pres/sat_vapor_pres_k.F90 .
./gradient_c2l.c: $(SRCROOT)shared/mosaic/gradient_c2l.c
	cp $(SRCROOT)shared/mosaic/gradient_c2l.c .
./zetac_clouds_W.F90: $(SRCROOT)atmos_param/sea_esf_rad/zetac_clouds_W.F90
	cp $(SRCROOT)atmos_param/sea_esf_rad/zetac_clouds_W.F90 .
./generic_CFC.F90: $(SRCROOT)ocean_shared/generic_tracers/generic_CFC.F90
	cp $(SRCROOT)ocean_shared/generic_tracers/generic_CFC.F90 .
./ocean_xlandinsert.F90: $(SRCROOT)mom5/ocean_param/sources/ocean_xlandinsert.F90
	cp $(SRCROOT)mom5/ocean_param/sources/ocean_xlandinsert.F90 .
./cmue_b.F90: $(SRCROOT)mom5/ocean_param/gotm-4.0/turbulence/cmue_b.F90
	cp $(SRCROOT)mom5/ocean_param/gotm-4.0/turbulence/cmue_b.F90 .
./stats_subs.F90: $(SRCROOT)atmos_param/clubb/CLUBB_core/stats_subs.F90
	cp $(SRCROOT)atmos_param/clubb/CLUBB_core/stats_subs.F90 .
./ocean_overflow_OFP.F90: $(SRCROOT)mom5/ocean_param/sources/ocean_overflow_OFP.F90
	cp $(SRCROOT)mom5/ocean_param/sources/ocean_overflow_OFP.F90 .
./mpp_parameter.F90: $(SRCROOT)shared/mpp/mpp_parameter.F90
	cp $(SRCROOT)shared/mpp/mpp_parameter.F90 .
./stats_variables.F90: $(SRCROOT)atmos_param/clubb/CLUBB_core/stats_variables.F90
	cp $(SRCROOT)atmos_param/clubb/CLUBB_core/stats_variables.F90 .
./lin_cloud_microphys.F90: $(SRCROOT)atmos_param/lin_cloud_microphys/lin_cloud_microphys.F90
	cp $(SRCROOT)atmos_param/lin_cloud_microphys/lin_cloud_microphys.F90 .
./cloudrad_diagnostics.F90: $(SRCROOT)atmos_param/sea_esf_rad/cloudrad_diagnostics.F90
	cp $(SRCROOT)atmos_param/sea_esf_rad/cloudrad_diagnostics.F90 .
./system_clock.h: /home/wenchang/tigress/CM2.5_tiger2_wy_killtc/exp/../src/shared/mpp/include/system_clock.h
	cp /home/wenchang/tigress/CM2.5_tiger2_wy_killtc/exp/../src/shared/mpp/include/system_clock.h .
./isccp_cloud_types.F: $(SRCROOT)atmos_param/cosp/icarus-scops-3.7/isccp_cloud_types.F
	cp $(SRCROOT)atmos_param/cosp/icarus-scops-3.7/isccp_cloud_types.F .
./MersenneTwister.F90: $(SRCROOT)shared/random_numbers/MersenneTwister.F90
	cp $(SRCROOT)shared/random_numbers/MersenneTwister.F90 .
./nfc.F90: $(SRCROOT)land_lad2/shared/nf_utils/nfc.F90
	cp $(SRCROOT)land_lad2/shared/nf_utils/nfc.F90 .
./mo_usrrxt.F90: $(SRCROOT)atmos_shared/tracer_driver/tropchem/mo_usrrxt.F90
	cp $(SRCROOT)atmos_shared/tracer_driver/tropchem/mo_usrrxt.F90 .
./vert_advection.F90: $(SRCROOT)atmos_shared/vert_advection/vert_advection.F90
	cp $(SRCROOT)atmos_shared/vert_advection/vert_advection.F90 .
./topography.F90: $(SRCROOT)shared/topography/topography.F90
	cp $(SRCROOT)shared/topography/topography.F90 .
./test_mpp.F90: $(SRCROOT)shared/mpp/test_mpp.F90
	cp $(SRCROOT)shared/mpp/test_mpp.F90 .
./mpp_do_updateV.h: /home/wenchang/tigress/CM2.5_tiger2_wy_killtc/exp/../src/shared/mpp/include/mpp_do_updateV.h
	cp /home/wenchang/tigress/CM2.5_tiger2_wy_killtc/exp/../src/shared/mpp/include/mpp_do_updateV.h .
./mpp_error_a_s.h: /home/wenchang/tigress/CM2.5_tiger2_wy_killtc/exp/../src/shared/mpp/include/mpp_error_a_s.h
	cp /home/wenchang/tigress/CM2.5_tiger2_wy_killtc/exp/../src/shared/mpp/include/mpp_error_a_s.h .
./ocean_coriolis.F90: $(SRCROOT)mom5/ocean_core/ocean_coriolis.F90
	cp $(SRCROOT)mom5/ocean_core/ocean_coriolis.F90 .
./affinity.c: $(SRCROOT)shared/mpp/affinity.c
	cp $(SRCROOT)shared/mpp/affinity.c .
./transport_matrix.F90: $(SRCROOT)mom5/ocean_tracers/transport_matrix.F90
	cp $(SRCROOT)mom5/ocean_tracers/transport_matrix.F90 .
./mpp_write_compressed.h: /home/wenchang/tigress/CM2.5_tiger2_wy_killtc/exp/../src/shared/mpp/include/mpp_write_compressed.h
	cp /home/wenchang/tigress/CM2.5_tiger2_wy_killtc/exp/../src/shared/mpp/include/mpp_write_compressed.h .
./fv_mapz.F90: $(SRCROOT)atmos_cubed_sphere/model/fv_mapz.F90
	cp $(SRCROOT)atmos_cubed_sphere/model/fv_mapz.F90 .
./mo_chemdr.F90: $(SRCROOT)atmos_shared/tracer_driver/tropchem/mo_chemdr.F90
	cp $(SRCROOT)atmos_shared/tracer_driver/tropchem/mo_chemdr.F90 .
./conv_closures.F90: $(SRCROOT)atmos_param/shallow_cu/conv_closures.F90
	cp $(SRCROOT)atmos_param/shallow_cu/conv_closures.F90 .
./donner_types.h: /home/wenchang/tigress/CM2.5_tiger2_wy_killtc/exp/../src/atmos_param/donner_deep/donner_types.h
	cp /home/wenchang/tigress/CM2.5_tiger2_wy_killtc/exp/../src/atmos_param/donner_deep/donner_types.h .
./advance_wp2_wp3_module.F90: $(SRCROOT)atmos_param/clubb/CLUBB_core/advance_wp2_wp3_module.F90
	cp $(SRCROOT)atmos_param/clubb/CLUBB_core/advance_wp2_wp3_module.F90 .
./output_grads.F90: $(SRCROOT)atmos_param/clubb/CLUBB_core/output_grads.F90
	cp $(SRCROOT)atmos_param/clubb/CLUBB_core/output_grads.F90 .
./mo_fastjx.F90: $(SRCROOT)atmos_shared/tracer_driver/tropchem/mo_fastjx.F90
	cp $(SRCROOT)atmos_shared/tracer_driver/tropchem/mo_fastjx.F90 .
./load_hydrometeor_classes.F90: $(SRCROOT)atmos_param/cosp/quickbeam/load_hydrometeor_classes.F90
	cp $(SRCROOT)atmos_param/cosp/quickbeam/load_hydrometeor_classes.F90 .
./ocean_workspace.F90: $(SRCROOT)mom5/ocean_core/ocean_workspace.F90
	cp $(SRCROOT)mom5/ocean_core/ocean_workspace.F90 .
./cosp_io.F90: $(SRCROOT)atmos_param/cosp/cosp_io.F90
	cp $(SRCROOT)atmos_param/cosp/cosp_io.F90 .
./cloud_zonal.F90: $(SRCROOT)atmos_param/cloud_zonal/cloud_zonal.F90
	cp $(SRCROOT)atmos_param/cloud_zonal/cloud_zonal.F90 .
./atmos_sulfur_hex.F90: $(SRCROOT)atmos_shared/tracer_driver/atmos_sulfur_hex.F90
	cp $(SRCROOT)atmos_shared/tracer_driver/atmos_sulfur_hex.F90 .
./ocean_shortwave_csiro.F90: $(SRCROOT)mom5/ocean_param/sources/ocean_shortwave_csiro.F90
	cp $(SRCROOT)mom5/ocean_param/sources/ocean_shortwave_csiro.F90 .
./ocmip2_cfc.F90: $(SRCROOT)mom5/ocean_bgc/ocmip2_cfc.F90
	cp $(SRCROOT)mom5/ocean_bgc/ocmip2_cfc.F90 .
./glacier.F90: $(SRCROOT)land_lad2/glacier/glacier.F90
	cp $(SRCROOT)land_lad2/glacier/glacier.F90 .
./mpp_sum_mpi.h: /home/wenchang/tigress/CM2.5_tiger2_wy_killtc/exp/../src/shared/mpp/include/mpp_sum_mpi.h
	cp /home/wenchang/tigress/CM2.5_tiger2_wy_killtc/exp/../src/shared/mpp/include/mpp_sum_mpi.h .
./mo_rodas_slv.F90: $(SRCROOT)atmos_shared/tracer_driver/tropchem/mo_rodas_slv.F90
	cp $(SRCROOT)atmos_shared/tracer_driver/tropchem/mo_rodas_slv.F90 .
./land_model.F90: $(SRCROOT)land_lad2/land_model.F90
	cp $(SRCROOT)land_lad2/land_model.F90 .
./oda_core.F90: $(SRCROOT)shared/oda_tools/oda_core.F90
	cp $(SRCROOT)shared/oda_tools/oda_core.F90 .
./conv_utilities.F90: $(SRCROOT)atmos_param/shallow_cu/conv_utilities.F90
	cp $(SRCROOT)atmos_param/shallow_cu/conv_utilities.F90 .
./get_cal_time.F90: $(SRCROOT)shared/time_manager/get_cal_time.F90
	cp $(SRCROOT)shared/time_manager/get_cal_time.F90 .
./hyper_diffusion_4th_ord.F90: $(SRCROOT)atmos_param/clubb/CLUBB_core/hyper_diffusion_4th_ord.F90
	cp $(SRCROOT)atmos_param/clubb/CLUBB_core/hyper_diffusion_4th_ord.F90 .
./fv_nesting.F90: $(SRCROOT)atmos_cubed_sphere/model/fv_nesting.F90
	cp $(SRCROOT)atmos_cubed_sphere/model/fv_nesting.F90 .
./compute_cpsi3.F90: $(SRCROOT)mom5/ocean_param/gotm-4.0/turbulence/compute_cpsi3.F90
	cp $(SRCROOT)mom5/ocean_param/gotm-4.0/turbulence/compute_cpsi3.F90 .
./ocean_form_drag.F90: $(SRCROOT)mom5/ocean_param/vertical/ocean_form_drag.F90
	cp $(SRCROOT)mom5/ocean_param/vertical/ocean_form_drag.F90 .
./cu_mo_trans.F90: $(SRCROOT)atmos_param/cu_mo_trans/cu_mo_trans.F90
	cp $(SRCROOT)atmos_param/cu_mo_trans/cu_mo_trans.F90 .
./cloud_spec.F90: $(SRCROOT)atmos_param/sea_esf_rad/cloud_spec.F90
	cp $(SRCROOT)atmos_param/sea_esf_rad/cloud_spec.F90 .
./lake_tile.F90: $(SRCROOT)land_lad2/lake/lake_tile.F90
	cp $(SRCROOT)land_lad2/lake/lake_tile.F90 .
./diag_output.F90: $(SRCROOT)shared/diag_manager/diag_output.F90
	cp $(SRCROOT)shared/diag_manager/diag_output.F90 .
./clip_explicit.F90: $(SRCROOT)atmos_param/clubb/CLUBB_core/clip_explicit.F90
	cp $(SRCROOT)atmos_param/clubb/CLUBB_core/clip_explicit.F90 .
./ocean_vert_tidal_test.F90: $(SRCROOT)mom5/ocean_param/vertical/ocean_vert_tidal_test.F90
	cp $(SRCROOT)mom5/ocean_param/vertical/ocean_vert_tidal_test.F90 .
./ocean_shortwave_jerlov.F90: $(SRCROOT)mom5/ocean_param/sources/ocean_shortwave_jerlov.F90
	cp $(SRCROOT)mom5/ocean_param/sources/ocean_shortwave_jerlov.F90 .
./stats_sfc.F90: $(SRCROOT)atmos_param/clubb/CLUBB_core/stats_sfc.F90
	cp $(SRCROOT)atmos_param/clubb/CLUBB_core/stats_sfc.F90 .
./modis_simulator.F90: $(SRCROOT)atmos_param/cosp/MODIS_simulator/modis_simulator.F90
	cp $(SRCROOT)atmos_param/cosp/MODIS_simulator/modis_simulator.F90 .
./saturation.F90: $(SRCROOT)atmos_param/clubb/CLUBB_core/saturation.F90
	cp $(SRCROOT)atmos_param/clubb/CLUBB_core/saturation.F90 .
./external_ic.F90: $(SRCROOT)atmos_cubed_sphere/tools/external_ic.F90
	cp $(SRCROOT)atmos_cubed_sphere/tools/external_ic.F90 .
./donner_nml.h: /home/wenchang/tigress/CM2.5_tiger2_wy_killtc/exp/../src/atmos_param/donner_deep/donner_nml.h
	cp /home/wenchang/tigress/CM2.5_tiger2_wy_killtc/exp/../src/atmos_param/donner_deep/donner_nml.h .
./read_remap_cohort_data_new.inc: /home/wenchang/tigress/CM2.5_tiger2_wy_killtc/exp/../src/land_lad2/vegetation/read_remap_cohort_data_new.inc
	cp /home/wenchang/tigress/CM2.5_tiger2_wy_killtc/exp/../src/land_lad2/vegetation/read_remap_cohort_data_new.inc .
./donner_cloud_model_k.F90: $(SRCROOT)atmos_param/donner_deep/donner_cloud_model_k.F90
	cp $(SRCROOT)atmos_param/donner_deep/donner_cloud_model_k.F90 .
./llnl_stats.F90: $(SRCROOT)atmos_param/cosp/llnl/llnl_stats.F90
	cp $(SRCROOT)atmos_param/cosp/llnl/llnl_stats.F90 .
./cosp_driver.F90: $(SRCROOT)atmos_param/cosp/cosp_driver.F90
	cp $(SRCROOT)atmos_param/cosp/cosp_driver.F90 .
./land_debug.F90: $(SRCROOT)land_lad2/shared/land_debug.F90
	cp $(SRCROOT)land_lad2/shared/land_debug.F90 .
./land_tile.F90: $(SRCROOT)land_lad2/land_tile.F90
	cp $(SRCROOT)land_lad2/land_tile.F90 .
./cmue_a.F90: $(SRCROOT)mom5/ocean_param/gotm-4.0/turbulence/cmue_a.F90
	cp $(SRCROOT)mom5/ocean_param/gotm-4.0/turbulence/cmue_a.F90 .
./fv_sg.F90: $(SRCROOT)atmos_cubed_sphere/model/fv_sg.F90
	cp $(SRCROOT)atmos_cubed_sphere/model/fv_sg.F90 .
./icarus.F90: $(SRCROOT)atmos_param/cosp/icarus-scops-3.7/icarus.F90
	cp $(SRCROOT)atmos_param/cosp/icarus-scops-3.7/icarus.F90 .
./atmos_sulfate.F90: $(SRCROOT)atmos_shared/tracer_driver/atmos_sulfate.F90
	cp $(SRCROOT)atmos_shared/tracer_driver/atmos_sulfate.F90 .
./damping_driver.F90: $(SRCROOT)atmos_param/damping_driver/damping_driver.F90
	cp $(SRCROOT)atmos_param/damping_driver/damping_driver.F90 .
./bm_omp.F90: $(SRCROOT)atmos_param/betts_miller/bm_omp.F90
	cp $(SRCROOT)atmos_param/betts_miller/bm_omp.F90 .
./getput.inc: /home/wenchang/tigress/CM2.5_tiger2_wy_killtc/exp/../src/land_lad2/shared/nf_utils/getput.inc
	cp /home/wenchang/tigress/CM2.5_tiger2_wy_killtc/exp/../src/land_lad2/shared/nf_utils/getput.inc .
./eqstate.F90: $(SRCROOT)mom5/ocean_param/gotm-4.0/util/eqstate.F90
	cp $(SRCROOT)mom5/ocean_param/gotm-4.0/util/eqstate.F90 .
./lmd_ipsl_stats.F90: $(SRCROOT)atmos_param/cosp/actsim/lmd_ipsl_stats.F90
	cp $(SRCROOT)atmos_param/cosp/actsim/lmd_ipsl_stats.F90 .
./moz.mat.F90: $(SRCROOT)atmos_shared/tracer_driver/tropchem/moz.mat.F90
	cp $(SRCROOT)atmos_shared/tracer_driver/tropchem/moz.mat.F90 .
./clip_semi_implicit.F90: $(SRCROOT)atmos_param/clubb/CLUBB_core/clip_semi_implicit.F90
	cp $(SRCROOT)atmos_param/clubb/CLUBB_core/clip_semi_implicit.F90 .
./parameters_microphys.F90: $(SRCROOT)atmos_param/clubb/CLUBB_core/parameters_microphys.F90
	cp $(SRCROOT)atmos_param/clubb/CLUBB_core/parameters_microphys.F90 .
./clouds.F90: $(SRCROOT)atmos_param/clouds/clouds.F90
	cp $(SRCROOT)atmos_param/clouds/clouds.F90 .
./mpp_domains_util.inc: /home/wenchang/tigress/CM2.5_tiger2_wy_killtc/exp/../src/shared/mpp/include/mpp_domains_util.inc
	cp /home/wenchang/tigress/CM2.5_tiger2_wy_killtc/exp/../src/shared/mpp/include/mpp_domains_util.inc .
./ocean_rivermix.F90: $(SRCROOT)mom5/ocean_param/sources/ocean_rivermix.F90
	cp $(SRCROOT)mom5/ocean_param/sources/ocean_rivermix.F90 .
./mpp_update_domains2D.h: /home/wenchang/tigress/CM2.5_tiger2_wy_killtc/exp/../src/shared/mpp/include/mpp_update_domains2D.h
	cp /home/wenchang/tigress/CM2.5_tiger2_wy_killtc/exp/../src/shared/mpp/include/mpp_update_domains2D.h .
./mean_adv.F90: $(SRCROOT)atmos_param/clubb/CLUBB_core/mean_adv.F90
	cp $(SRCROOT)atmos_param/clubb/CLUBB_core/mean_adv.F90 .
./test_mpp_pset.F90: $(SRCROOT)shared/mpp/test_mpp_pset.F90
	cp $(SRCROOT)shared/mpp/test_mpp_pset.F90 .
./fms_io.F90: $(SRCROOT)shared/fms/fms_io.F90
	cp $(SRCROOT)shared/fms/fms_io.F90 .
./strat_chem_driver.F90: $(SRCROOT)atmos_shared/tracer_driver/stratchem/strat_chem_driver.F90
	cp $(SRCROOT)atmos_shared/tracer_driver/stratchem/strat_chem_driver.F90 .
./ocean_convect.F90: $(SRCROOT)mom5/ocean_param/vertical/ocean_convect.F90
	cp $(SRCROOT)mom5/ocean_param/vertical/ocean_convect.F90 .
./donner_deep_clouds_W.F90: $(SRCROOT)atmos_param/sea_esf_rad/donner_deep_clouds_W.F90
	cp $(SRCROOT)atmos_param/sea_esf_rad/donner_deep_clouds_W.F90 .
./congvec.H: /home/wenchang/tigress/CM2.5_tiger2_wy_killtc/exp/../src/atmos_param/cosp/icarus-scops-3.7/congvec.H
	cp /home/wenchang/tigress/CM2.5_tiger2_wy_killtc/exp/../src/atmos_param/cosp/icarus-scops-3.7/congvec.H .
./cg_drag.F90: $(SRCROOT)atmos_param/cg_drag/cg_drag.F90
	cp $(SRCROOT)atmos_param/cg_drag/cg_drag.F90 .
./wet_deposition_0D.F90: $(SRCROOT)atmos_param/donner_deep/wet_deposition_0D.F90
	cp $(SRCROOT)atmos_param/donner_deep/wet_deposition_0D.F90 .
./ice_bergs.F90: $(SRCROOT)ice_sis/ice_bergs.F90
	cp $(SRCROOT)ice_sis/ice_bergs.F90 .
./mo_jpl.F90: $(SRCROOT)atmos_shared/tracer_driver/tropchem/mo_jpl.F90
	cp $(SRCROOT)atmos_shared/tracer_driver/tropchem/mo_jpl.F90 .
./rdparm.F90: $(SRCROOT)atmos_param/fsrad/rdparm.F90
	cp $(SRCROOT)atmos_param/fsrad/rdparm.F90 .
./ocean_barotropic.F90: $(SRCROOT)mom5/ocean_core/ocean_barotropic.F90
	cp $(SRCROOT)mom5/ocean_core/ocean_barotropic.F90 .
./generic_tracer_utils.F90: $(SRCROOT)ocean_shared/generic_tracers/generic_tracer_utils.F90
	cp $(SRCROOT)ocean_shared/generic_tracers/generic_tracer_utils.F90 .
./interp.h: /home/wenchang/tigress/CM2.5_tiger2_wy_killtc/exp/../src/shared/mosaic/interp.h
	cp /home/wenchang/tigress/CM2.5_tiger2_wy_killtc/exp/../src/shared/mosaic/interp.h .
./ocean_bgc_restore.F90: $(SRCROOT)mom5/ocean_bgc/ocean_bgc_restore.F90
	cp $(SRCROOT)mom5/ocean_bgc/ocean_bgc_restore.F90 .
./aerosol_params.F90: $(SRCROOT)atmos_shared/tracer_driver/aer_ccn_act/aerosol_params.F90
	cp $(SRCROOT)atmos_shared/tracer_driver/aer_ccn_act/aerosol_params.F90 .
./ocean_bih_friction.F90: $(SRCROOT)mom5/ocean_param/lateral/ocean_bih_friction.F90
	cp $(SRCROOT)mom5/ocean_param/lateral/ocean_bih_friction.F90 .
./pf_to_mr.F: $(SRCROOT)atmos_param/cosp/llnl/pf_to_mr.F
	cp $(SRCROOT)atmos_param/cosp/llnl/pf_to_mr.F .
./vegn_cohort_io.F90: $(SRCROOT)land_lad2/vegetation/vegn_cohort_io.F90
	cp $(SRCROOT)land_lad2/vegetation/vegn_cohort_io.F90 .
./soil_tile.F90: $(SRCROOT)land_lad2/soil/soil_tile.F90
	cp $(SRCROOT)land_lad2/soil/soil_tile.F90 .
./ocean_wave_type.F90: $(SRCROOT)mom5/ocean_wave/ocean_wave_type.F90
	cp $(SRCROOT)mom5/ocean_wave/ocean_wave_type.F90 .
./load_mie_table.F90: $(SRCROOT)atmos_param/cosp/quickbeam/load_mie_table.F90
	cp $(SRCROOT)atmos_param/cosp/quickbeam/load_mie_table.F90 .
./mask.F90: $(SRCROOT)ice_sis/mask.F90
	cp $(SRCROOT)ice_sis/mask.F90 .
./CLUBB_3D_var.F90: $(SRCROOT)atmos_param/clubb/CLUBB_3D_var.F90
	cp $(SRCROOT)atmos_param/clubb/CLUBB_3D_var.F90 .
./mpp_data_sma.inc: /home/wenchang/tigress/CM2.5_tiger2_wy_killtc/exp/../src/shared/mpp/include/mpp_data_sma.inc
	cp /home/wenchang/tigress/CM2.5_tiger2_wy_killtc/exp/../src/shared/mpp/include/mpp_data_sma.inc .
./create_xgrid.h: /home/wenchang/tigress/CM2.5_tiger2_wy_killtc/exp/../src/shared/mosaic/create_xgrid.h
	cp /home/wenchang/tigress/CM2.5_tiger2_wy_killtc/exp/../src/shared/mosaic/create_xgrid.h .
./lw_gases_stdtf.F90: $(SRCROOT)atmos_param/sea_esf_rad/lw_gases_stdtf.F90
	cp $(SRCROOT)atmos_param/sea_esf_rad/lw_gases_stdtf.F90 .
./ocean_nphysicsC.F90: $(SRCROOT)mom5/ocean_param/neutral/ocean_nphysicsC.F90
	cp $(SRCROOT)mom5/ocean_param/neutral/ocean_nphysicsC.F90 .
./ocean_overflow.F90: $(SRCROOT)mom5/ocean_param/sources/ocean_overflow.F90
	cp $(SRCROOT)mom5/ocean_param/sources/ocean_overflow.F90 .
./compute_rist.F90: $(SRCROOT)mom5/ocean_param/gotm-4.0/turbulence/compute_rist.F90
	cp $(SRCROOT)mom5/ocean_param/gotm-4.0/turbulence/compute_rist.F90 .
./test_cases.F90: $(SRCROOT)atmos_cubed_sphere/tools/test_cases.F90
	cp $(SRCROOT)atmos_cubed_sphere/tools/test_cases.F90 .
./aerosol_cloud.F90: $(SRCROOT)atmos_param/strat_cloud/aerosol_cloud.F90
	cp $(SRCROOT)atmos_param/strat_cloud/aerosol_cloud.F90 .
./parameters_tunable.F90: $(SRCROOT)atmos_param/clubb/CLUBB_core/parameters_tunable.F90
	cp $(SRCROOT)atmos_param/clubb/CLUBB_core/parameters_tunable.F90 .
./ocean_bbc.F90: $(SRCROOT)mom5/ocean_core/ocean_bbc.F90
	cp $(SRCROOT)mom5/ocean_core/ocean_bbc.F90 .
./lscale_cond.F90: $(SRCROOT)atmos_param/lscale_cond/lscale_cond.F90
	cp $(SRCROOT)atmos_param/lscale_cond/lscale_cond.F90 .
./shortwave.F90: $(SRCROOT)atmos_param/fsrad/shortwave.F90
	cp $(SRCROOT)atmos_param/fsrad/shortwave.F90 .
./mpp_do_redistribute.h: /home/wenchang/tigress/CM2.5_tiger2_wy_killtc/exp/../src/shared/mpp/include/mpp_do_redistribute.h
	cp /home/wenchang/tigress/CM2.5_tiger2_wy_killtc/exp/../src/shared/mpp/include/mpp_do_redistribute.h .
./fv_tracer2d.F90: $(SRCROOT)atmos_cubed_sphere/model/fv_tracer2d.F90
	cp $(SRCROOT)atmos_cubed_sphere/model/fv_tracer2d.F90 .
./mo_read_sim_chm.F90: $(SRCROOT)atmos_shared/tracer_driver/tropchem/mo_read_sim_chm.F90
	cp $(SRCROOT)atmos_shared/tracer_driver/tropchem/mo_read_sim_chm.F90 .
./ocean_lap_friction.F90: $(SRCROOT)mom5/ocean_param/lateral/ocean_lap_friction.F90
	cp $(SRCROOT)mom5/ocean_param/lateral/ocean_lap_friction.F90 .
./ocean_density.F90: $(SRCROOT)mom5/ocean_core/ocean_density.F90
	cp $(SRCROOT)mom5/ocean_core/ocean_density.F90 .
./vegn_tile.F90: $(SRCROOT)land_lad2/vegetation/vegn_tile.F90
	cp $(SRCROOT)land_lad2/vegetation/vegn_tile.F90 .
./original_fms_rad.F90: $(SRCROOT)atmos_param/sea_esf_rad/original_fms_rad.F90
	cp $(SRCROOT)atmos_param/sea_esf_rad/original_fms_rad.F90 .
./cmue_sg.F90: $(SRCROOT)mom5/ocean_param/gotm-4.0/turbulence/cmue_sg.F90
	cp $(SRCROOT)mom5/ocean_param/gotm-4.0/turbulence/cmue_sg.F90 .
./ocean_generic_tracer.F90: $(SRCROOT)mom5/ocean_bgc/ocean_generic_tracer.F90
	cp $(SRCROOT)mom5/ocean_bgc/ocean_generic_tracer.F90 .
./cumulus_closure_k.F90: $(SRCROOT)atmos_param/donner_deep/cumulus_closure_k.F90
	cp $(SRCROOT)atmos_param/donner_deep/cumulus_closure_k.F90 .
./ocean_sponges_velocity.F90: $(SRCROOT)mom5/ocean_param/sources/ocean_sponges_velocity.F90
	cp $(SRCROOT)mom5/ocean_param/sources/ocean_sponges_velocity.F90 .
./ocean_passive.F90: $(SRCROOT)mom5/ocean_tracers/ocean_passive.F90
	cp $(SRCROOT)mom5/ocean_tracers/ocean_passive.F90 .
./longwave_fluxes.F90: $(SRCROOT)atmos_param/sea_esf_rad/longwave_fluxes.F90
	cp $(SRCROOT)atmos_param/sea_esf_rad/longwave_fluxes.F90 .
./sphum.F90: $(SRCROOT)land_lad2/shared/sphum.F90
	cp $(SRCROOT)land_lad2/shared/sphum.F90 .
./ocean_vert_const.F90: $(SRCROOT)mom5/ocean_param/vertical/ocean_vert_const.F90
	cp $(SRCROOT)mom5/ocean_param/vertical/ocean_vert_const.F90 .
./axis_utils.F90: $(SRCROOT)shared/axis_utils/axis_utils.F90
	cp $(SRCROOT)shared/axis_utils/axis_utils.F90 .
./mpp_sum.inc: /home/wenchang/tigress/CM2.5_tiger2_wy_killtc/exp/../src/shared/mpp/include/mpp_sum.inc
	cp /home/wenchang/tigress/CM2.5_tiger2_wy_killtc/exp/../src/shared/mpp/include/mpp_sum.inc .
./conv_plumes_k.F90: $(SRCROOT)atmos_param/shallow_cu/conv_plumes_k.F90
	cp $(SRCROOT)atmos_param/shallow_cu/conv_plumes_k.F90 .
./ice_spec.F90: $(SRCROOT)ice_sis/ice_spec.F90
	cp $(SRCROOT)ice_sis/ice_spec.F90 .
./donner_rad_k.F90: $(SRCROOT)atmos_param/donner_deep/donner_rad_k.F90
	cp $(SRCROOT)atmos_param/donner_deep/donner_rad_k.F90 .
./mpp_util_mpi.inc: /home/wenchang/tigress/CM2.5_tiger2_wy_killtc/exp/../src/shared/mpp/include/mpp_util_mpi.inc
	cp /home/wenchang/tigress/CM2.5_tiger2_wy_killtc/exp/../src/shared/mpp/include/mpp_util_mpi.inc .
./drifters_push.h: /home/wenchang/tigress/CM2.5_tiger2_wy_killtc/exp/../src/shared/drifters/drifters_push.h
	cp /home/wenchang/tigress/CM2.5_tiger2_wy_killtc/exp/../src/shared/drifters/drifters_push.h .
./mo_setsox.F90: $(SRCROOT)atmos_shared/tracer_driver/tropchem/mo_setsox.F90
	cp $(SRCROOT)atmos_shared/tracer_driver/tropchem/mo_setsox.F90 .
./lagrange.F90: $(SRCROOT)mom5/ocean_param/gotm-4.0/util/lagrange.F90
	cp $(SRCROOT)mom5/ocean_param/gotm-4.0/util/lagrange.F90 .
./aer_ccn_act_k.F90: $(SRCROOT)atmos_shared/tracer_driver/aer_ccn_act/aer_ccn_act_k.F90
	cp $(SRCROOT)atmos_shared/tracer_driver/aer_ccn_act/aer_ccn_act_k.F90 .
./mpp_do_checkV.h: /home/wenchang/tigress/CM2.5_tiger2_wy_killtc/exp/../src/shared/mpp/include/mpp_do_checkV.h
	cp /home/wenchang/tigress/CM2.5_tiger2_wy_killtc/exp/../src/shared/mpp/include/mpp_do_checkV.h .
./cmue_c.F90: $(SRCROOT)mom5/ocean_param/gotm-4.0/turbulence/cmue_c.F90
	cp $(SRCROOT)mom5/ocean_param/gotm-4.0/turbulence/cmue_c.F90 .
./diffusion.F90: $(SRCROOT)atmos_param/clubb/CLUBB_core/diffusion.F90
	cp $(SRCROOT)atmos_param/clubb/CLUBB_core/diffusion.F90 .
./radiation_driver.F90: $(SRCROOT)atmos_param/radiation_driver/radiation_driver.F90
	cp $(SRCROOT)atmos_param/radiation_driver/radiation_driver.F90 .
./pdf_parameter_module.F90: $(SRCROOT)atmos_param/clubb/CLUBB_core/pdf_parameter_module.F90
	cp $(SRCROOT)atmos_param/clubb/CLUBB_core/pdf_parameter_module.F90 .
./horiz_interp_conserve.F90: $(SRCROOT)shared/horiz_interp/horiz_interp_conserve.F90
	cp $(SRCROOT)shared/horiz_interp/horiz_interp_conserve.F90 .
./mpp_sum_sma.h: /home/wenchang/tigress/CM2.5_tiger2_wy_killtc/exp/../src/shared/mpp/include/mpp_sum_sma.h
	cp /home/wenchang/tigress/CM2.5_tiger2_wy_killtc/exp/../src/shared/mpp/include/mpp_sum_sma.h .
./strat_cloud_legacy.F90: $(SRCROOT)atmos_param/strat_cloud/strat_cloud_legacy.F90
	cp $(SRCROOT)atmos_param/strat_cloud/strat_cloud_legacy.F90 .
./init_hydro.F90: $(SRCROOT)atmos_cubed_sphere/tools/init_hydro.F90
	cp $(SRCROOT)atmos_cubed_sphere/tools/init_hydro.F90 .
./tropchem_driver.F90: $(SRCROOT)atmos_shared/tracer_driver/tropchem/tropchem_driver.F90
	cp $(SRCROOT)atmos_shared/tracer_driver/tropchem/tropchem_driver.F90 .
./longwave_tables.F90: $(SRCROOT)atmos_param/sea_esf_rad/longwave_tables.F90
	cp $(SRCROOT)atmos_param/sea_esf_rad/longwave_tables.F90 .
./atmosphere.F90: $(SRCROOT)atmos_cubed_sphere/driver/coupled/atmosphere.F90
	cp $(SRCROOT)atmos_cubed_sphere/driver/coupled/atmosphere.F90 .
./mpp_write.h: /home/wenchang/tigress/CM2.5_tiger2_wy_killtc/exp/../src/shared/mpp/include/mpp_write.h
	cp /home/wenchang/tigress/CM2.5_tiger2_wy_killtc/exp/../src/shared/mpp/include/mpp_write.h .
./grey_radiation.F90: $(SRCROOT)atmos_param/grey_radiation/grey_radiation.F90
	cp $(SRCROOT)atmos_param/grey_radiation/grey_radiation.F90 .
./diag_manager.F90: $(SRCROOT)shared/diag_manager/diag_manager.F90
	cp $(SRCROOT)shared/diag_manager/diag_manager.F90 .
./ocean_parameters.F90: $(SRCROOT)mom5/ocean_core/ocean_parameters.F90
	cp $(SRCROOT)mom5/ocean_core/ocean_parameters.F90 .
./ocean_adv_vel_diag.F90: $(SRCROOT)mom5/ocean_diag/ocean_adv_vel_diag.F90
	cp $(SRCROOT)mom5/ocean_diag/ocean_adv_vel_diag.F90 .
./ocean_domains.F90: $(SRCROOT)mom5/ocean_core/ocean_domains.F90
	cp $(SRCROOT)mom5/ocean_core/ocean_domains.F90 .
./mpp_transmit_sma.h: /home/wenchang/tigress/CM2.5_tiger2_wy_killtc/exp/../src/shared/mpp/include/mpp_transmit_sma.h
	cp /home/wenchang/tigress/CM2.5_tiger2_wy_killtc/exp/../src/shared/mpp/include/mpp_transmit_sma.h .
./ocean_util.F90: $(SRCROOT)mom5/ocean_core/ocean_util.F90
	cp $(SRCROOT)mom5/ocean_core/ocean_util.F90 .
./fms.F90: $(SRCROOT)shared/fms/fms.F90
	cp $(SRCROOT)shared/fms/fms.F90 .
./snow.F90: $(SRCROOT)land_lad2/snow/snow.F90
	cp $(SRCROOT)land_lad2/snow/snow.F90 .
./ocean_sponges_tracer.F90: $(SRCROOT)mom5/ocean_param/sources/ocean_sponges_tracer.F90
	cp $(SRCROOT)mom5/ocean_param/sources/ocean_sponges_tracer.F90 .
./alt_cloud.F90: $(SRCROOT)atmos_param/clubb/alt_cloud.F90
	cp $(SRCROOT)atmos_param/clubb/alt_cloud.F90 .
./fv_restart.F90: $(SRCROOT)atmos_cubed_sphere/tools/fv_restart.F90
	cp $(SRCROOT)atmos_cubed_sphere/tools/fv_restart.F90 .
./astronomy.F90: $(SRCROOT)shared/astronomy/astronomy.F90
	cp $(SRCROOT)shared/astronomy/astronomy.F90 .
./river_physics.F90: $(SRCROOT)land_lad2/river/river_physics.F90
	cp $(SRCROOT)land_lad2/river/river_physics.F90 .
./ocean_nphysics_util.F90: $(SRCROOT)mom5/ocean_param/neutral/ocean_nphysics_util.F90
	cp $(SRCROOT)mom5/ocean_param/neutral/ocean_nphysics_util.F90 .
./xbt_drop_rate_adjust.f90: $(SRCROOT)shared/oda_tools/xbt_drop_rate_adjust.f90
	cp $(SRCROOT)shared/oda_tools/xbt_drop_rate_adjust.f90 .
./ocean_vert_pp.F90: $(SRCROOT)mom5/ocean_param/vertical/ocean_vert_pp.F90
	cp $(SRCROOT)mom5/ocean_param/vertical/ocean_vert_pp.F90 .
./ispralength.F90: $(SRCROOT)mom5/ocean_param/gotm-4.0/turbulence/ispralength.F90
	cp $(SRCROOT)mom5/ocean_param/gotm-4.0/turbulence/ispralength.F90 .
./mo_exp_slv.F90: $(SRCROOT)atmos_shared/tracer_driver/tropchem/mo_exp_slv.F90
	cp $(SRCROOT)atmos_shared/tracer_driver/tropchem/mo_exp_slv.F90 .
./mpp_global_reduce.h: /home/wenchang/tigress/CM2.5_tiger2_wy_killtc/exp/../src/shared/mpp/include/mpp_global_reduce.h
	cp /home/wenchang/tigress/CM2.5_tiger2_wy_killtc/exp/../src/shared/mpp/include/mpp_global_reduce.h .
./error_code.F90: $(SRCROOT)atmos_param/clubb/CLUBB_core/error_code.F90
	cp $(SRCROOT)atmos_param/clubb/CLUBB_core/error_code.F90 .
./aer_in_act.F90: $(SRCROOT)atmos_shared/tracer_driver/aer_ccn_act/aer_in_act.F90
	cp $(SRCROOT)atmos_shared/tracer_driver/aer_ccn_act/aer_in_act.F90 .
./fsrad.F90: $(SRCROOT)atmos_param/fsrad/fsrad.F90
	cp $(SRCROOT)atmos_param/fsrad/fsrad.F90 .
./scops.F: $(SRCROOT)atmos_param/cosp/icarus-scops-3.7/scops.F
	cp $(SRCROOT)atmos_param/cosp/icarus-scops-3.7/scops.F .
./rad_diag.F90: $(SRCROOT)atmos_param/fsrad/rad_diag.F90
	cp $(SRCROOT)atmos_param/fsrad/rad_diag.F90 .
./memuse.c: $(SRCROOT)shared/memutils/memuse.c
	cp $(SRCROOT)shared/memutils/memuse.c .
./cosp_radar.F90: $(SRCROOT)atmos_param/cosp/llnl/cosp_radar.F90
	cp $(SRCROOT)atmos_param/cosp/llnl/cosp_radar.F90 .
./debug.inc: /home/wenchang/tigress/CM2.5_tiger2_wy_killtc/exp/../src/land_lad2/canopy_air/../shared/debug.inc
	cp /home/wenchang/tigress/CM2.5_tiger2_wy_killtc/exp/../src/land_lad2/canopy_air/../shared/debug.inc .
./ocean_lapcst_friction.F90: $(SRCROOT)mom5/ocean_param/lateral/ocean_lapcst_friction.F90
	cp $(SRCROOT)mom5/ocean_param/lateral/ocean_lapcst_friction.F90 .
./mpp_comm_nocomm.inc: /home/wenchang/tigress/CM2.5_tiger2_wy_killtc/exp/../src/shared/mpp/include/mpp_comm_nocomm.inc
	cp /home/wenchang/tigress/CM2.5_tiger2_wy_killtc/exp/../src/shared/mpp/include/mpp_comm_nocomm.inc .
./ocean_frazil.F90: $(SRCROOT)mom5/ocean_tracers/ocean_frazil.F90
	cp $(SRCROOT)mom5/ocean_tracers/ocean_frazil.F90 .
./mcm_swtbls.F90: $(SRCROOT)atmos_param/fsrad/mcm_swtbls.F90
	cp $(SRCROOT)atmos_param/fsrad/mcm_swtbls.F90 .
./mpp_chksum.h: /home/wenchang/tigress/CM2.5_tiger2_wy_killtc/exp/../src/shared/mpp/include/mpp_chksum.h
	cp /home/wenchang/tigress/CM2.5_tiger2_wy_killtc/exp/../src/shared/mpp/include/mpp_chksum.h .
./strat_cloud_utilities.F90: $(SRCROOT)atmos_param/strat_cloud/strat_cloud_utilities.F90
	cp $(SRCROOT)atmos_param/strat_cloud/strat_cloud_utilities.F90 .
./grid_class.F90: $(SRCROOT)atmos_param/clubb/CLUBB_core/grid_class.F90
	cp $(SRCROOT)atmos_param/clubb/CLUBB_core/grid_class.F90 .
./donner_deep_k.F90: $(SRCROOT)atmos_param/donner_deep/donner_deep_k.F90
	cp $(SRCROOT)atmos_param/donner_deep/donner_deep_k.F90 .
./generic_tracer.F90: $(SRCROOT)ocean_shared/generic_tracers/generic_tracer.F90
	cp $(SRCROOT)ocean_shared/generic_tracers/generic_tracer.F90 .
./microphysics.F90: $(SRCROOT)atmos_param/strat_cloud/microphysics.F90
	cp $(SRCROOT)atmos_param/strat_cloud/microphysics.F90 .
./strat_chem_model.F90: $(SRCROOT)atmos_shared/tracer_driver/stratchem/strat_chem_model.F90
	cp $(SRCROOT)atmos_shared/tracer_driver/stratchem/strat_chem_model.F90 .
./mo_hook.F90: $(SRCROOT)atmos_shared/tracer_driver/tropchem/mo_hook.F90
	cp $(SRCROOT)atmos_shared/tracer_driver/tropchem/mo_hook.F90 .
./tke_turb.F90: $(SRCROOT)atmos_param/tke_turb/tke_turb.F90
	cp $(SRCROOT)atmos_param/tke_turb/tke_turb.F90 .
./moistproc_kernels.F90: $(SRCROOT)atmos_param/moist_processes/moistproc_kernels.F90
	cp $(SRCROOT)atmos_param/moist_processes/moistproc_kernels.F90 .
./diag_data.F90: $(SRCROOT)shared/diag_manager/diag_data.F90
	cp $(SRCROOT)shared/diag_manager/diag_data.F90 .
./kpp.F90: $(SRCROOT)mom5/ocean_param/gotm-4.0/turbulence/kpp.F90
	cp $(SRCROOT)mom5/ocean_param/gotm-4.0/turbulence/kpp.F90 .
./mpp_memutils.F90: $(SRCROOT)shared/mpp/mpp_memutils.F90
	cp $(SRCROOT)shared/mpp/mpp_memutils.F90 .
./ocean_nphysicsB.F90: $(SRCROOT)mom5/ocean_param/neutral/ocean_nphysicsB.F90
	cp $(SRCROOT)mom5/ocean_param/neutral/ocean_nphysicsB.F90 .
./diag_clouds_W.F90: $(SRCROOT)atmos_param/sea_esf_rad/diag_clouds_W.F90
	cp $(SRCROOT)atmos_param/sea_esf_rad/diag_clouds_W.F90 .
./donner_lscloud_k.F90: $(SRCROOT)atmos_param/donner_deep/donner_lscloud_k.F90
	cp $(SRCROOT)atmos_param/donner_deep/donner_lscloud_k.F90 .
./test_fms_io.F90: $(SRCROOT)shared/fms/test_fms_io.F90
	cp $(SRCROOT)shared/fms/test_fms_io.F90 .
./horiz_interp_spherical.F90: $(SRCROOT)shared/horiz_interp/horiz_interp_spherical.F90
	cp $(SRCROOT)shared/horiz_interp/horiz_interp_spherical.F90 .
./simple_pdf.F90: $(SRCROOT)atmos_param/strat_cloud/simple_pdf.F90
	cp $(SRCROOT)atmos_param/strat_cloud/simple_pdf.F90 .
./glac_tile.F90: $(SRCROOT)land_lad2/glacier/glac_tile.F90
	cp $(SRCROOT)land_lad2/glacier/glac_tile.F90 .
./ocean_vert_kpp_mom4p1.F90: $(SRCROOT)mom5/ocean_param/vertical/ocean_vert_kpp_mom4p1.F90
	cp $(SRCROOT)mom5/ocean_param/vertical/ocean_vert_kpp_mom4p1.F90 .
./atmos_soa.F90: $(SRCROOT)atmos_shared/tracer_driver/atmos_soa.F90
	cp $(SRCROOT)atmos_shared/tracer_driver/atmos_soa.F90 .
./ocean_tpm.F90: $(SRCROOT)mom5/ocean_tracers/ocean_tpm.F90
	cp $(SRCROOT)mom5/ocean_tracers/ocean_tpm.F90 .
./cosp_isccp_simulator.F90: $(SRCROOT)atmos_param/cosp/cosp_isccp_simulator.F90
	cp $(SRCROOT)atmos_param/cosp/cosp_isccp_simulator.F90 .
./tracer_manager.F90: $(SRCROOT)shared/tracer_manager/tracer_manager.F90
	cp $(SRCROOT)shared/tracer_manager/tracer_manager.F90 .
./array_lib.F90: $(SRCROOT)atmos_param/cosp/quickbeam/array_lib.F90
	cp $(SRCROOT)atmos_param/cosp/quickbeam/array_lib.F90 .
./gmres_cache.F90: $(SRCROOT)atmos_param/clubb/CLUBB_core/gmres_cache.F90
	cp $(SRCROOT)atmos_param/clubb/CLUBB_core/gmres_cache.F90 .
./radar_simulator.F90: $(SRCROOT)atmos_param/cosp/quickbeam/radar_simulator.F90
	cp $(SRCROOT)atmos_param/cosp/quickbeam/radar_simulator.F90 .
./cosp.F90: $(SRCROOT)atmos_param/cosp/cosp.F90
	cp $(SRCROOT)atmos_param/cosp/cosp.F90 .
./uptake.F90: $(SRCROOT)land_lad2/soil/uptake.F90
	cp $(SRCROOT)land_lad2/soil/uptake.F90 .
./ocean_riverspread.F90: $(SRCROOT)mom5/ocean_param/sources/ocean_riverspread.F90
	cp $(SRCROOT)mom5/ocean_param/sources/ocean_riverspread.F90 .
./horiz_interp_bilinear.F90: $(SRCROOT)shared/horiz_interp/horiz_interp_bilinear.F90
	cp $(SRCROOT)shared/horiz_interp/horiz_interp_bilinear.F90 .
./vegn_data.F90: $(SRCROOT)land_lad2/vegetation/vegn_data.F90
	cp $(SRCROOT)land_lad2/vegetation/vegn_data.F90 .
./topo_rough.F90: $(SRCROOT)land_lad2/topo_rough/topo_rough.F90
	cp $(SRCROOT)land_lad2/topo_rough/topo_rough.F90 .
./mosaic_util.c: $(SRCROOT)shared/mosaic/mosaic_util.c
	cp $(SRCROOT)shared/mosaic/mosaic_util.c .
./xgrid.F90: $(SRCROOT)shared/exchange/xgrid.F90
	cp $(SRCROOT)shared/exchange/xgrid.F90 .
./donner_cape_k.F90: $(SRCROOT)atmos_param/donner_deep/donner_cape_k.F90
	cp $(SRCROOT)atmos_param/donner_deep/donner_cape_k.F90 .
./advance_sclrm_Nd_module.F90: $(SRCROOT)atmos_param/clubb/advance_sclrm_Nd_module.F90
	cp $(SRCROOT)atmos_param/clubb/advance_sclrm_Nd_module.F90 .
./fv_nudge.F90: $(SRCROOT)atmos_cubed_sphere/tools/fv_nudge.F90
	cp $(SRCROOT)atmos_cubed_sphere/tools/fv_nudge.F90 .
./mpp_update_domains2D_nonblock.h: /home/wenchang/tigress/CM2.5_tiger2_wy_killtc/exp/../src/shared/mpp/include/mpp_update_domains2D_nonblock.h
	cp /home/wenchang/tigress/CM2.5_tiger2_wy_killtc/exp/../src/shared/mpp/include/mpp_update_domains2D_nonblock.h .
./ocmip2_biotic.F90: $(SRCROOT)mom5/ocean_bgc/ocmip2_biotic.F90
	cp $(SRCROOT)mom5/ocean_bgc/ocmip2_biotic.F90 .
./vegn_dynamics.F90: $(SRCROOT)land_lad2/vegetation/vegn_dynamics.F90
	cp $(SRCROOT)land_lad2/vegetation/vegn_dynamics.F90 .
./constant.h: /home/wenchang/tigress/CM2.5_tiger2_wy_killtc/exp/../src/shared/mosaic/constant.h
	cp /home/wenchang/tigress/CM2.5_tiger2_wy_killtc/exp/../src/shared/mosaic/constant.h .
./land_io.F90: $(SRCROOT)land_lad2/shared/land_io.F90
	cp $(SRCROOT)land_lad2/shared/land_io.F90 .
./cloud_rad.F90: $(SRCROOT)atmos_param/cloud_rad/cloud_rad.F90
	cp $(SRCROOT)atmos_param/cloud_rad/cloud_rad.F90 .
./drifters_input.F90: $(SRCROOT)shared/drifters/drifters_input.F90
	cp $(SRCROOT)shared/drifters/drifters_input.F90 .
./ocean_types.F90: $(SRCROOT)mom5/ocean_core/ocean_types.F90
	cp $(SRCROOT)mom5/ocean_core/ocean_types.F90 .
./stats_zm.F90: $(SRCROOT)atmos_param/clubb/CLUBB_core/stats_zm.F90
	cp $(SRCROOT)atmos_param/clubb/CLUBB_core/stats_zm.F90 .
./ocean_rough.F90: $(SRCROOT)ice_param/ocean_rough.F90
	cp $(SRCROOT)ice_param/ocean_rough.F90 .
./strat_nml.h: /home/wenchang/tigress/CM2.5_tiger2_wy_killtc/exp/../src/atmos_param/strat_cloud/strat_nml.h
	cp /home/wenchang/tigress/CM2.5_tiger2_wy_killtc/exp/../src/atmos_param/strat_cloud/strat_nml.h .
./ice_thm.F90: $(SRCROOT)ice_sis/ice_thm.F90
	cp $(SRCROOT)ice_sis/ice_thm.F90 .
./mpp_reduce_nocomm.h: /home/wenchang/tigress/CM2.5_tiger2_wy_killtc/exp/../src/shared/mpp/include/mpp_reduce_nocomm.h
	cp /home/wenchang/tigress/CM2.5_tiger2_wy_killtc/exp/../src/shared/mpp/include/mpp_reduce_nocomm.h .
./mpp_io_util.inc: /home/wenchang/tigress/CM2.5_tiger2_wy_killtc/exp/../src/shared/mpp/include/mpp_io_util.inc
	cp /home/wenchang/tigress/CM2.5_tiger2_wy_killtc/exp/../src/shared/mpp/include/mpp_io_util.inc .
./coupler_main.F90: $(SRCROOT)coupler/coupler_main.F90
	cp $(SRCROOT)coupler/coupler_main.F90 .
./cloudrad_package.F90: $(SRCROOT)atmos_param/sea_esf_rad/cloudrad_package.F90
	cp $(SRCROOT)atmos_param/sea_esf_rad/cloudrad_package.F90 .
./topo_drag.F90: $(SRCROOT)atmos_param/topo_drag/topo_drag.F90
	cp $(SRCROOT)atmos_param/topo_drag/topo_drag.F90 .
./mpp_error_s_a.h: /home/wenchang/tigress/CM2.5_tiger2_wy_killtc/exp/../src/shared/mpp/include/mpp_error_s_a.h
	cp /home/wenchang/tigress/CM2.5_tiger2_wy_killtc/exp/../src/shared/mpp/include/mpp_error_s_a.h .
./atmos_convection_tracer.F90: $(SRCROOT)atmos_shared/tracer_driver/atmos_convection_tracer.F90
	cp $(SRCROOT)atmos_shared/tracer_driver/atmos_convection_tracer.F90 .
./write_ocean_data.F90: $(SRCROOT)shared/oda_tools/write_ocean_data.F90
	cp $(SRCROOT)shared/oda_tools/write_ocean_data.F90 .
./co2int.F90: $(SRCROOT)atmos_param/fsrad/co2int.F90
	cp $(SRCROOT)atmos_param/fsrad/co2int.F90 .
./uw_clouds_W.F90: $(SRCROOT)atmos_param/sea_esf_rad/uw_clouds_W.F90
	cp $(SRCROOT)atmos_param/sea_esf_rad/uw_clouds_W.F90 .
./strat_cloud.F90: $(SRCROOT)atmos_param/strat_cloud/strat_cloud.F90
	cp $(SRCROOT)atmos_param/strat_cloud/strat_cloud.F90 .
./cosp_simulator.F90: $(SRCROOT)atmos_param/cosp/cosp_simulator.F90
	cp $(SRCROOT)atmos_param/cosp/cosp_simulator.F90 .
./boundary.F90: $(SRCROOT)atmos_cubed_sphere/model/boundary.F90
	cp $(SRCROOT)atmos_cubed_sphere/model/boundary.F90 .
./time_interp.F90: $(SRCROOT)shared/time_interp/time_interp.F90
	cp $(SRCROOT)shared/time_interp/time_interp.F90 .
./mixing_length.F90: $(SRCROOT)atmos_param/clubb/CLUBB_core/mixing_length.F90
	cp $(SRCROOT)atmos_param/clubb/CLUBB_core/mixing_length.F90 .
./horiz_interp.F90: $(SRCROOT)shared/horiz_interp/horiz_interp.F90
	cp $(SRCROOT)shared/horiz_interp/horiz_interp.F90 .
./mrgrnk.F90: $(SRCROOT)atmos_param/cosp/quickbeam/mrgrnk.F90
	cp $(SRCROOT)atmos_param/cosp/quickbeam/mrgrnk.F90 .
./specified_clouds_W.F90: $(SRCROOT)atmos_param/sea_esf_rad/specified_clouds_W.F90
	cp $(SRCROOT)atmos_param/sea_esf_rad/specified_clouds_W.F90 .
./river_type.F90: $(SRCROOT)land_lad2/river/river_type.F90
	cp $(SRCROOT)land_lad2/river/river_type.F90 .
./mpp_util.inc: /home/wenchang/tigress/CM2.5_tiger2_wy_killtc/exp/../src/shared/mpp/include/mpp_util.inc
	cp /home/wenchang/tigress/CM2.5_tiger2_wy_killtc/exp/../src/shared/mpp/include/mpp_util.inc .
./stats_rad_zt.F90: $(SRCROOT)atmos_param/clubb/CLUBB_core/stats_rad_zt.F90
	cp $(SRCROOT)atmos_param/clubb/CLUBB_core/stats_rad_zt.F90 .
./atmos_nudge.F90: $(SRCROOT)atmos_shared/atmos_nudge/atmos_nudge.F90
	cp $(SRCROOT)atmos_shared/atmos_nudge/atmos_nudge.F90 .
./ice_nucl.F90: $(SRCROOT)atmos_shared/tracer_driver/aer_ccn_act/ice_nucl.F90
	cp $(SRCROOT)atmos_shared/tracer_driver/aer_ccn_act/ice_nucl.F90 .
./ocean_nphysics_diff.F90: $(SRCROOT)mom5/ocean_param/neutral/ocean_nphysics_diff.F90
	cp $(SRCROOT)mom5/ocean_param/neutral/ocean_nphysics_diff.F90 .
./ocean_blob_static_free.F90: $(SRCROOT)mom5/ocean_blobs/ocean_blob_static_free.F90
	cp $(SRCROOT)mom5/ocean_blobs/ocean_blob_static_free.F90 .
./ocean_nphysics.F90: $(SRCROOT)mom5/ocean_param/neutral/ocean_nphysics.F90
	cp $(SRCROOT)mom5/ocean_param/neutral/ocean_nphysics.F90 .
./create_xgrid.c: $(SRCROOT)shared/mosaic/create_xgrid.c
	cp $(SRCROOT)shared/mosaic/create_xgrid.c .
./mpp_transmit_mpi.h: /home/wenchang/tigress/CM2.5_tiger2_wy_killtc/exp/../src/shared/mpp/include/mpp_transmit_mpi.h
	cp /home/wenchang/tigress/CM2.5_tiger2_wy_killtc/exp/../src/shared/mpp/include/mpp_transmit_mpi.h .
./land_numerics.F90: $(SRCROOT)land_lad2/shared/land_numerics.F90
	cp $(SRCROOT)land_lad2/shared/land_numerics.F90 .
./rotstayn_klein_mp.F90: $(SRCROOT)atmos_param/strat_cloud/rotstayn_klein_mp.F90
	cp $(SRCROOT)atmos_param/strat_cloud/rotstayn_klein_mp.F90 .
./deep_conv.F90: $(SRCROOT)atmos_param/shallow_cu/deep_conv.F90
	cp $(SRCROOT)atmos_param/shallow_cu/deep_conv.F90 .
./tp_core.F90: $(SRCROOT)atmos_cubed_sphere/model/tp_core.F90
	cp $(SRCROOT)atmos_cubed_sphere/model/tp_core.F90 .
./model_flags.F90: $(SRCROOT)atmos_param/clubb/CLUBB_core/model_flags.F90
	cp $(SRCROOT)atmos_param/clubb/CLUBB_core/model_flags.F90 .
./tkeeq.F90: $(SRCROOT)mom5/ocean_param/gotm-4.0/turbulence/tkeeq.F90
	cp $(SRCROOT)mom5/ocean_param/gotm-4.0/turbulence/tkeeq.F90 .
./interpolator.F90: $(SRCROOT)atmos_shared/interpolator/interpolator.F90
	cp $(SRCROOT)atmos_shared/interpolator/interpolator.F90 .
./mpp_io_read.inc: /home/wenchang/tigress/CM2.5_tiger2_wy_killtc/exp/../src/shared/mpp/include/mpp_io_read.inc
	cp /home/wenchang/tigress/CM2.5_tiger2_wy_killtc/exp/../src/shared/mpp/include/mpp_io_read.inc .
./lake.F90: $(SRCROOT)land_lad2/lake/lake.F90
	cp $(SRCROOT)land_lad2/lake/lake.F90 .
./ozone.F90: $(SRCROOT)atmos_param/sea_esf_rad/ozone.F90
	cp $(SRCROOT)atmos_param/sea_esf_rad/ozone.F90 .
./math_lib.F90: $(SRCROOT)atmos_param/cosp/quickbeam/math_lib.F90
	cp $(SRCROOT)atmos_param/cosp/quickbeam/math_lib.F90 .
./ocean_bihcst_friction.F90: $(SRCROOT)mom5/ocean_param/lateral/ocean_bihcst_friction.F90
	cp $(SRCROOT)mom5/ocean_param/lateral/ocean_bihcst_friction.F90 .
./atmos_tracer_utilities.F90: $(SRCROOT)atmos_shared/tracer_driver/atmos_tracer_utilities.F90
	cp $(SRCROOT)atmos_shared/tracer_driver/atmos_tracer_utilities.F90 .
./production.F90: $(SRCROOT)mom5/ocean_param/gotm-4.0/turbulence/production.F90
	cp $(SRCROOT)mom5/ocean_param/gotm-4.0/turbulence/production.F90 .
./mo_imp_slv.F90: $(SRCROOT)atmos_shared/tracer_driver/tropchem/mo_imp_slv.F90
	cp $(SRCROOT)atmos_shared/tracer_driver/tropchem/mo_imp_slv.F90 .
./donner_types.F90: $(SRCROOT)atmos_param/donner_deep/donner_types.F90
	cp $(SRCROOT)atmos_param/donner_deep/donner_types.F90 .
./donner_lite_k.F90: $(SRCROOT)atmos_param/donner_deep/donner_lite_k.F90
	cp $(SRCROOT)atmos_param/donner_deep/donner_lite_k.F90 .
./longwave_clouds.F90: $(SRCROOT)atmos_param/sea_esf_rad/longwave_clouds.F90
	cp $(SRCROOT)atmos_param/sea_esf_rad/longwave_clouds.F90 .
./radiation_diag.F90: $(SRCROOT)atmos_param/sea_esf_rad/radiation_diag.F90
	cp $(SRCROOT)atmos_param/sea_esf_rad/radiation_diag.F90 .
./ocean_albedo.F90: $(SRCROOT)ice_param/ocean_albedo.F90
	cp $(SRCROOT)ice_param/ocean_albedo.F90 .
./generic_ERGOM.F90: $(SRCROOT)ocean_shared/generic_tracers/generic_ERGOM.F90
	cp $(SRCROOT)ocean_shared/generic_tracers/generic_ERGOM.F90 .
./mosaic.F90: $(SRCROOT)shared/mosaic/mosaic.F90
	cp $(SRCROOT)shared/mosaic/mosaic.F90 .
./internal_wave.F90: $(SRCROOT)mom5/ocean_param/gotm-4.0/turbulence/internal_wave.F90
	cp $(SRCROOT)mom5/ocean_param/gotm-4.0/turbulence/internal_wave.F90 .
./stats_type.F90: $(SRCROOT)atmos_param/clubb/CLUBB_core/stats_type.F90
	cp $(SRCROOT)atmos_param/clubb/CLUBB_core/stats_type.F90 .
./lapack_wrap.F90: $(SRCROOT)atmos_param/clubb/CLUBB_core/lapack_wrap.F90
	cp $(SRCROOT)atmos_param/clubb/CLUBB_core/lapack_wrap.F90 .
./gaussian_topog.F90: $(SRCROOT)shared/topography/gaussian_topog.F90
	cp $(SRCROOT)shared/topography/gaussian_topog.F90 .
./version.h: /home/wenchang/tigress/CM2.5_tiger2_wy_killtc/exp/../src/mom5/ocean_param/gotm-4.0/include/version.h
	cp /home/wenchang/tigress/CM2.5_tiger2_wy_killtc/exp/../src/mom5/ocean_param/gotm-4.0/include/version.h .
./atmos_tracer_driver.F90: $(SRCROOT)atmos_shared/tracer_driver/atmos_tracer_driver.F90
	cp $(SRCROOT)atmos_shared/tracer_driver/atmos_tracer_driver.F90 .
./diff_face.F90: $(SRCROOT)mom5/ocean_param/gotm-4.0/util/diff_face.F90
	cp $(SRCROOT)mom5/ocean_param/gotm-4.0/util/diff_face.F90 .
./atmos_lib.F90: $(SRCROOT)atmos_param/cosp/quickbeam/atmos_lib.F90
	cp $(SRCROOT)atmos_param/cosp/quickbeam/atmos_lib.F90 .
./vegn_static_override.F90: $(SRCROOT)land_lad2/vegetation/vegn_static_override.F90
	cp $(SRCROOT)land_lad2/vegetation/vegn_static_override.F90 .
./mpp_define_nest_domains.inc: /home/wenchang/tigress/CM2.5_tiger2_wy_killtc/exp/../src/shared/mpp/include/mpp_define_nest_domains.inc
	cp /home/wenchang/tigress/CM2.5_tiger2_wy_killtc/exp/../src/shared/mpp/include/mpp_define_nest_domains.inc .
./polysvp.F90: $(SRCROOT)atmos_param/strat_cloud/polysvp.F90
	cp $(SRCROOT)atmos_param/strat_cloud/polysvp.F90 .
./donner_utilities_k.F90: $(SRCROOT)atmos_param/donner_deep/donner_utilities_k.F90
	cp $(SRCROOT)atmos_param/donner_deep/donner_utilities_k.F90 .
./esfsw_driver.F90: $(SRCROOT)atmos_param/sea_esf_rad/esfsw_driver.F90
	cp $(SRCROOT)atmos_param/sea_esf_rad/esfsw_driver.F90 .
./mpp_domains_comm.inc: /home/wenchang/tigress/CM2.5_tiger2_wy_killtc/exp/../src/shared/mpp/include/mpp_domains_comm.inc
	cp /home/wenchang/tigress/CM2.5_tiger2_wy_killtc/exp/../src/shared/mpp/include/mpp_domains_comm.inc .
./ocean_xlandmix.F90: $(SRCROOT)mom5/ocean_param/sources/ocean_xlandmix.F90
	cp $(SRCROOT)mom5/ocean_param/sources/ocean_xlandmix.F90 .
./cloud_obs.F90: $(SRCROOT)atmos_param/cloud_obs/cloud_obs.F90
	cp $(SRCROOT)atmos_param/cloud_obs/cloud_obs.F90 .
./cosp_rttov_simulator.F90: $(SRCROOT)atmos_param/cosp/cosp_rttov_simulator.F90
	cp $(SRCROOT)atmos_param/cosp/cosp_rttov_simulator.F90 .
./mpp_gather.h: /home/wenchang/tigress/CM2.5_tiger2_wy_killtc/exp/../src/shared/mpp/include/mpp_gather.h
	cp /home/wenchang/tigress/CM2.5_tiger2_wy_killtc/exp/../src/shared/mpp/include/mpp_gather.h .
./cosp_utils.F90: $(SRCROOT)atmos_param/cosp/cosp_utils.F90
	cp $(SRCROOT)atmos_param/cosp/cosp_utils.F90 .
./mpp_transmit.inc: /home/wenchang/tigress/CM2.5_tiger2_wy_killtc/exp/../src/shared/mpp/include/mpp_transmit.inc
	cp /home/wenchang/tigress/CM2.5_tiger2_wy_killtc/exp/../src/shared/mpp/include/mpp_transmit.inc .
./ocean_age_tracer.F90: $(SRCROOT)mom5/ocean_tracers/ocean_age_tracer.F90
	cp $(SRCROOT)mom5/ocean_tracers/ocean_age_tracer.F90 .
./vegn_harvesting.F90: $(SRCROOT)land_lad2/vegetation/vegn_harvesting.F90
	cp $(SRCROOT)land_lad2/vegetation/vegn_harvesting.F90 .
./detr_ice_num.F90: $(SRCROOT)atmos_param/moist_processes/detr_ice_num.F90
	cp $(SRCROOT)atmos_param/moist_processes/detr_ice_num.F90 .
./ocean_model.F90: $(SRCROOT)mom5/ocean_core/ocean_model.F90
	cp $(SRCROOT)mom5/ocean_core/ocean_model.F90 .
./util.F90: $(SRCROOT)mom5/ocean_param/gotm-4.0/util/util.F90
	cp $(SRCROOT)mom5/ocean_param/gotm-4.0/util/util.F90 .
./cana_tile.F90: $(SRCROOT)land_lad2/canopy_air/cana_tile.F90
	cp $(SRCROOT)land_lad2/canopy_air/cana_tile.F90 .
./ocean_lapgen_friction.F90: $(SRCROOT)mom5/ocean_param/lateral/ocean_lapgen_friction.F90
	cp $(SRCROOT)mom5/ocean_param/lateral/ocean_lapgen_friction.F90 .
./random_numbers.F90: $(SRCROOT)shared/random_numbers/random_numbers.F90
	cp $(SRCROOT)shared/random_numbers/random_numbers.F90 .
./ocean_residency_ml.F90: $(SRCROOT)mom5/ocean_tracers/ocean_residency_ml.F90
	cp $(SRCROOT)mom5/ocean_tracers/ocean_residency_ml.F90 .
./vegn_cohort_io.inc: /home/wenchang/tigress/CM2.5_tiger2_wy_killtc/exp/../src/land_lad2/vegetation/vegn_cohort_io.inc
	cp /home/wenchang/tigress/CM2.5_tiger2_wy_killtc/exp/../src/land_lad2/vegetation/vegn_cohort_io.inc .
./sw_core.F90: $(SRCROOT)atmos_cubed_sphere/model/sw_core.F90
	cp $(SRCROOT)atmos_cubed_sphere/model/sw_core.F90 .
./ocean_vert_chen.F90: $(SRCROOT)mom5/ocean_param/vertical/ocean_vert_chen.F90
	cp $(SRCROOT)mom5/ocean_param/vertical/ocean_vert_chen.F90 .
./read_data_3d.inc: /home/wenchang/tigress/CM2.5_tiger2_wy_killtc/exp/../src/shared/fms/read_data_3d.inc
	cp /home/wenchang/tigress/CM2.5_tiger2_wy_killtc/exp/../src/shared/fms/read_data_3d.inc .
./ocean_bihgen_friction.F90: $(SRCROOT)mom5/ocean_param/lateral/ocean_bihgen_friction.F90
	cp $(SRCROOT)mom5/ocean_param/lateral/ocean_bihgen_friction.F90 .
./sealw99.F90: $(SRCROOT)atmos_param/sea_esf_rad/sealw99.F90
	cp $(SRCROOT)atmos_param/sea_esf_rad/sealw99.F90 .
./mosaic_util.h: /home/wenchang/tigress/CM2.5_tiger2_wy_killtc/exp/../src/shared/mosaic/mosaic_util.h
	cp /home/wenchang/tigress/CM2.5_tiger2_wy_killtc/exp/../src/shared/mosaic/mosaic_util.h .
./mo_fphoto.F90: $(SRCROOT)atmos_shared/tracer_driver/tropchem/mo_fphoto.F90
	cp $(SRCROOT)atmos_shared/tracer_driver/tropchem/mo_fphoto.F90 .
./mgrp_prscr_clds.F90: $(SRCROOT)atmos_param/sea_esf_rad/mgrp_prscr_clds.F90
	cp $(SRCROOT)atmos_param/sea_esf_rad/mgrp_prscr_clds.F90 .
./radiative_gases.F90: $(SRCROOT)atmos_param/sea_esf_rad/radiative_gases.F90
	cp $(SRCROOT)atmos_param/sea_esf_rad/radiative_gases.F90 .
./tkealgebraic.F90: $(SRCROOT)mom5/ocean_param/gotm-4.0/turbulence/tkealgebraic.F90
	cp $(SRCROOT)mom5/ocean_param/gotm-4.0/turbulence/tkealgebraic.F90 .
./oda_driver.F90: $(SRCROOT)mom5/ocean_core/oda_driver.F90
	cp $(SRCROOT)mom5/ocean_core/oda_driver.F90 .
./drifters_set_field.h: /home/wenchang/tigress/CM2.5_tiger2_wy_killtc/exp/../src/shared/drifters/drifters_set_field.h
	cp /home/wenchang/tigress/CM2.5_tiger2_wy_killtc/exp/../src/shared/drifters/drifters_set_field.h .
./ocean_blob_dynamic_bottom.F90: $(SRCROOT)mom5/ocean_blobs/ocean_blob_dynamic_bottom.F90
	cp $(SRCROOT)mom5/ocean_blobs/ocean_blob_dynamic_bottom.F90 .
./mpp_read_2Ddecomp.h: /home/wenchang/tigress/CM2.5_tiger2_wy_killtc/exp/../src/shared/mpp/include/mpp_read_2Ddecomp.h
	cp /home/wenchang/tigress/CM2.5_tiger2_wy_killtc/exp/../src/shared/mpp/include/mpp_read_2Ddecomp.h .
./ocean_shortwave.F90: $(SRCROOT)mom5/ocean_param/sources/ocean_shortwave.F90
	cp $(SRCROOT)mom5/ocean_param/sources/ocean_shortwave.F90 .
./bulkphys_rad.F90: $(SRCROOT)atmos_param/sea_esf_rad/bulkphys_rad.F90
	cp $(SRCROOT)atmos_param/sea_esf_rad/bulkphys_rad.F90 .
./fv_ada_nudge.F90: $(SRCROOT)atmos_cubed_sphere/tools/fv_ada_nudge.F90
	cp $(SRCROOT)atmos_cubed_sphere/tools/fv_ada_nudge.F90 .
./ocean_vert_kpp_mom4p0.F90: $(SRCROOT)mom5/ocean_param/vertical/ocean_vert_kpp_mom4p0.F90
	cp $(SRCROOT)mom5/ocean_param/vertical/ocean_vert_kpp_mom4p0.F90 .
./edt.F90: $(SRCROOT)atmos_param/edt/edt.F90
	cp $(SRCROOT)atmos_param/edt/edt.F90 .
./strat_clouds_W.F90: $(SRCROOT)atmos_param/sea_esf_rad/strat_clouds_W.F90
	cp $(SRCROOT)atmos_param/sea_esf_rad/strat_clouds_W.F90 .
./cosp_defs.H: /home/wenchang/tigress/CM2.5_tiger2_wy_killtc/exp/../src/atmos_param/cosp/cosp_defs.H
	cp /home/wenchang/tigress/CM2.5_tiger2_wy_killtc/exp/../src/atmos_param/cosp/cosp_defs.H .
./ocean_blob_util.F90: $(SRCROOT)mom5/ocean_blobs/ocean_blob_util.F90
	cp $(SRCROOT)mom5/ocean_blobs/ocean_blob_util.F90 .
./vert_diff_driver.F90: $(SRCROOT)atmos_param/vert_diff_driver/vert_diff_driver.F90
	cp $(SRCROOT)atmos_param/vert_diff_driver/vert_diff_driver.F90 .
./mpp_comm.inc: /home/wenchang/tigress/CM2.5_tiger2_wy_killtc/exp/../src/shared/mpp/include/mpp_comm.inc
	cp /home/wenchang/tigress/CM2.5_tiger2_wy_killtc/exp/../src/shared/mpp/include/mpp_comm.inc .
./rot.F90: $(SRCROOT)ice_sis/rot.F90
	cp $(SRCROOT)ice_sis/rot.F90 .
./betaDistribution.F90: $(SRCROOT)atmos_param/cloud_generator/betaDistribution.F90
	cp $(SRCROOT)atmos_param/cloud_generator/betaDistribution.F90 .
./donner_deep_miz.F90: $(SRCROOT)atmos_param/donner_deep/donner_deep_miz.F90
	cp $(SRCROOT)atmos_param/donner_deep/donner_deep_miz.F90 .
./diag_axis.F90: $(SRCROOT)shared/diag_manager/diag_axis.F90
	cp $(SRCROOT)shared/diag_manager/diag_axis.F90 .
./surface_varnce_module.F90: $(SRCROOT)atmos_param/clubb/CLUBB_core/surface_varnce_module.F90
	cp $(SRCROOT)atmos_param/clubb/CLUBB_core/surface_varnce_module.F90 .
./ocean_velocity_advect.F90: $(SRCROOT)mom5/ocean_core/ocean_velocity_advect.F90
	cp $(SRCROOT)mom5/ocean_core/ocean_velocity_advect.F90 .
./ocean_lapcgrid_friction.F90: $(SRCROOT)mom5/ocean_param/lateral/ocean_lapcgrid_friction.F90
	cp $(SRCROOT)mom5/ocean_param/lateral/ocean_lapcgrid_friction.F90 .
./micro_mg.F90: $(SRCROOT)atmos_param/strat_cloud/micro_mg.F90
	cp $(SRCROOT)atmos_param/strat_cloud/micro_mg.F90 .
./nfu.F90: $(SRCROOT)land_lad2/shared/nf_utils/nfu.F90
	cp $(SRCROOT)land_lad2/shared/nf_utils/nfu.F90 .
./write_data.inc: /home/wenchang/tigress/CM2.5_tiger2_wy_killtc/exp/../src/shared/fms/write_data.inc
	cp /home/wenchang/tigress/CM2.5_tiger2_wy_killtc/exp/../src/shared/fms/write_data.inc .
./longwave_driver.F90: $(SRCROOT)atmos_param/sea_esf_rad/longwave_driver.F90
	cp $(SRCROOT)atmos_param/sea_esf_rad/longwave_driver.F90 .
./csr_matrix_class_3array.F90: $(SRCROOT)atmos_param/clubb/CLUBB_core/csr_matrix_class_3array.F90
	cp $(SRCROOT)atmos_param/clubb/CLUBB_core/csr_matrix_class_3array.F90 .
./q2over2eq.F90: $(SRCROOT)mom5/ocean_param/gotm-4.0/turbulence/q2over2eq.F90
	cp $(SRCROOT)mom5/ocean_param/gotm-4.0/turbulence/q2over2eq.F90 .
./numerical_check.F90: $(SRCROOT)atmos_param/clubb/CLUBB_core/numerical_check.F90
	cp $(SRCROOT)atmos_param/clubb/CLUBB_core/numerical_check.F90 .
./ocean_tracer_util.F90: $(SRCROOT)mom5/ocean_diag/ocean_tracer_util.F90
	cp $(SRCROOT)mom5/ocean_diag/ocean_tracer_util.F90 .
./ocean_velocity.F90: $(SRCROOT)mom5/ocean_core/ocean_velocity.F90
	cp $(SRCROOT)mom5/ocean_core/ocean_velocity.F90 .
./fv_fill.F90: $(SRCROOT)atmos_cubed_sphere/model/fv_fill.F90
	cp $(SRCROOT)atmos_cubed_sphere/model/fv_fill.F90 .
./ocean_tracer_advect.F90: $(SRCROOT)mom5/ocean_tracers/ocean_tracer_advect.F90
	cp $(SRCROOT)mom5/ocean_tracers/ocean_tracer_advect.F90 .
./horiz_interp_bicubic.F90: $(SRCROOT)shared/horiz_interp/horiz_interp_bicubic.F90
	cp $(SRCROOT)shared/horiz_interp/horiz_interp_bicubic.F90 .
./donner_meso_k.F90: $(SRCROOT)atmos_param/donner_deep/donner_meso_k.F90
	cp $(SRCROOT)atmos_param/donner_deep/donner_meso_k.F90 .
./esfsw_parameters.F90: $(SRCROOT)atmos_param/sea_esf_rad/esfsw_parameters.F90
	cp $(SRCROOT)atmos_param/sea_esf_rad/esfsw_parameters.F90 .
./ocean_operators.F90: $(SRCROOT)mom5/ocean_core/ocean_operators.F90
	cp $(SRCROOT)mom5/ocean_core/ocean_operators.F90 .
./aerosol.F90: $(SRCROOT)atmos_param/sea_esf_rad/aerosol.F90
	cp $(SRCROOT)atmos_param/sea_esf_rad/aerosol.F90 .
./matrix_operations.F90: $(SRCROOT)atmos_param/clubb/CLUBB_core/matrix_operations.F90
	cp $(SRCROOT)atmos_param/clubb/CLUBB_core/matrix_operations.F90 .
./sat_vapor_pres.F90: $(SRCROOT)shared/sat_vapor_pres/sat_vapor_pres.F90
	cp $(SRCROOT)shared/sat_vapor_pres/sat_vapor_pres.F90 .
./ocean_nphysics_tensor.F90: $(SRCROOT)mom5/ocean_param/neutral/ocean_nphysics_tensor.F90
	cp $(SRCROOT)mom5/ocean_param/neutral/ocean_nphysics_tensor.F90 .
./mpp_io.F90: $(SRCROOT)shared/mpp/mpp_io.F90
	cp $(SRCROOT)shared/mpp/mpp_io.F90 .
./aer_ccn_act.F90: $(SRCROOT)atmos_shared/tracer_driver/aer_ccn_act/aer_ccn_act.F90
	cp $(SRCROOT)atmos_shared/tracer_driver/aer_ccn_act/aer_ccn_act.F90 .
./ocean_vert_util.F90: $(SRCROOT)mom5/ocean_param/vertical/ocean_vert_util.F90
	cp $(SRCROOT)mom5/ocean_param/vertical/ocean_vert_util.F90 .
./fs_profile.F90: $(SRCROOT)atmos_param/fsrad/fs_profile.F90
	cp $(SRCROOT)atmos_param/fsrad/fs_profile.F90 .
./drifters_core.F90: $(SRCROOT)shared/drifters/drifters_core.F90
	cp $(SRCROOT)shared/drifters/drifters_core.F90 .
./morrison_gettelman_microp.F90: $(SRCROOT)atmos_param/strat_cloud/morrison_gettelman_microp.F90
	cp $(SRCROOT)atmos_param/strat_cloud/morrison_gettelman_microp.F90 .
./nf_utils.F90: $(SRCROOT)land_lad2/shared/nf_utils/nf_utils.F90
	cp $(SRCROOT)land_lad2/shared/nf_utils/nf_utils.F90 .
./donner_deep.F90: $(SRCROOT)atmos_param/donner_deep/donner_deep.F90
	cp $(SRCROOT)atmos_param/donner_deep/donner_deep.F90 .
./diag_table.F90: $(SRCROOT)shared/diag_manager/diag_table.F90
	cp $(SRCROOT)shared/diag_manager/diag_table.F90 .
./ocean_blob.F90: $(SRCROOT)mom5/ocean_blobs/ocean_blob.F90
	cp $(SRCROOT)mom5/ocean_blobs/ocean_blob.F90 .
./rad_utilities.F90: $(SRCROOT)atmos_param/sea_esf_rad/rad_utilities.F90
	cp $(SRCROOT)atmos_param/sea_esf_rad/rad_utilities.F90 .
./mpp_do_update_nonblock.h: /home/wenchang/tigress/CM2.5_tiger2_wy_killtc/exp/../src/shared/mpp/include/mpp_do_update_nonblock.h
	cp /home/wenchang/tigress/CM2.5_tiger2_wy_killtc/exp/../src/shared/mpp/include/mpp_do_update_nonblock.h .
./generic_BLING.F90: $(SRCROOT)ocean_shared/generic_tracers/generic_BLING.F90
	cp $(SRCROOT)ocean_shared/generic_tracers/generic_BLING.F90 .
./mpp_data_nocomm.inc: /home/wenchang/tigress/CM2.5_tiger2_wy_killtc/exp/../src/shared/mpp/include/mpp_data_nocomm.inc
	cp /home/wenchang/tigress/CM2.5_tiger2_wy_killtc/exp/../src/shared/mpp/include/mpp_data_nocomm.inc .
./mpp_update_nest_domains.h: /home/wenchang/tigress/CM2.5_tiger2_wy_killtc/exp/../src/shared/mpp/include/mpp_update_nest_domains.h
	cp /home/wenchang/tigress/CM2.5_tiger2_wy_killtc/exp/../src/shared/mpp/include/mpp_update_nest_domains.h .
./mpp_domains.F90: $(SRCROOT)shared/mpp/mpp_domains.F90
	cp $(SRCROOT)shared/mpp/mpp_domains.F90 .
./ocean_blob_diag.F90: $(SRCROOT)mom5/ocean_blobs/ocean_blob_diag.F90
	cp $(SRCROOT)mom5/ocean_blobs/ocean_blob_diag.F90 .
./mpp_chksum_scalar.h: /home/wenchang/tigress/CM2.5_tiger2_wy_killtc/exp/../src/shared/mpp/include/mpp_chksum_scalar.h
	cp /home/wenchang/tigress/CM2.5_tiger2_wy_killtc/exp/../src/shared/mpp/include/mpp_chksum_scalar.h .
./dry_adj.F90: $(SRCROOT)atmos_param/dry_adj/dry_adj.F90
	cp $(SRCROOT)atmos_param/dry_adj/dry_adj.F90 .
./generic_TOPAZ.F90: $(SRCROOT)ocean_shared/generic_tracers/generic_TOPAZ.F90
	cp $(SRCROOT)ocean_shared/generic_tracers/generic_TOPAZ.F90 .
./parameters_model.F90: $(SRCROOT)atmos_param/clubb/CLUBB_core/parameters_model.F90
	cp $(SRCROOT)atmos_param/clubb/CLUBB_core/parameters_model.F90 .
./diag_cloud.F90: $(SRCROOT)atmos_param/diag_cloud/diag_cloud.F90
	cp $(SRCROOT)atmos_param/diag_cloud/diag_cloud.F90 .
./stats_rad_zm.F90: $(SRCROOT)atmos_param/clubb/CLUBB_core/stats_rad_zm.F90
	cp $(SRCROOT)atmos_param/clubb/CLUBB_core/stats_rad_zm.F90 .
./advance_xp2_xpyp_module.F90: $(SRCROOT)atmos_param/clubb/CLUBB_core/advance_xp2_xpyp_module.F90
	cp $(SRCROOT)atmos_param/clubb/CLUBB_core/advance_xp2_xpyp_module.F90 .
./read_data_2d.inc: /home/wenchang/tigress/CM2.5_tiger2_wy_killtc/exp/../src/shared/fms/read_data_2d.inc
	cp /home/wenchang/tigress/CM2.5_tiger2_wy_killtc/exp/../src/shared/fms/read_data_2d.inc .
./mpp_domains_define.inc: /home/wenchang/tigress/CM2.5_tiger2_wy_killtc/exp/../src/shared/mpp/include/mpp_domains_define.inc
	cp /home/wenchang/tigress/CM2.5_tiger2_wy_killtc/exp/../src/shared/mpp/include/mpp_domains_define.inc .
./ocean_memory.h: /home/wenchang/tigress/CM2.5_tiger2_wy_killtc/exp/../src/mom5/ocean_core/ocean_memory.h
	cp /home/wenchang/tigress/CM2.5_tiger2_wy_killtc/exp/../src/mom5/ocean_core/ocean_memory.h .
./ocean_diagnostics.F90: $(SRCROOT)mom5/ocean_diag/ocean_diagnostics.F90
	cp $(SRCROOT)mom5/ocean_diag/ocean_diagnostics.F90 .
./potentialml.F90: $(SRCROOT)mom5/ocean_param/gotm-4.0/turbulence/potentialml.F90
	cp $(SRCROOT)mom5/ocean_param/gotm-4.0/turbulence/potentialml.F90 .
./fv_diagnostics.F90: $(SRCROOT)atmos_cubed_sphere/tools/fv_diagnostics.F90
	cp $(SRCROOT)atmos_cubed_sphere/tools/fv_diagnostics.F90 .
./memutils.F90: $(SRCROOT)shared/memutils/memutils.F90
	cp $(SRCROOT)shared/memutils/memutils.F90 .
./mo_setinv.F90: $(SRCROOT)atmos_shared/tracer_driver/tropchem/mo_setinv.F90
	cp $(SRCROOT)atmos_shared/tracer_driver/tropchem/mo_setinv.F90 .
./fk_craig.F90: $(SRCROOT)mom5/ocean_param/gotm-4.0/turbulence/fk_craig.F90
	cp $(SRCROOT)mom5/ocean_param/gotm-4.0/turbulence/fk_craig.F90 .
./mpp_scatter.h: /home/wenchang/tigress/CM2.5_tiger2_wy_killtc/exp/../src/shared/mpp/include/mpp_scatter.h
	cp /home/wenchang/tigress/CM2.5_tiger2_wy_killtc/exp/../src/shared/mpp/include/mpp_scatter.h .
./gridinterpol.F90: $(SRCROOT)mom5/ocean_param/gotm-4.0/util/gridinterpol.F90
	cp $(SRCROOT)mom5/ocean_param/gotm-4.0/util/gridinterpol.F90 .
./mpp_reduce_mpi.h: /home/wenchang/tigress/CM2.5_tiger2_wy_killtc/exp/../src/shared/mpp/include/mpp_reduce_mpi.h
	cp /home/wenchang/tigress/CM2.5_tiger2_wy_killtc/exp/../src/shared/mpp/include/mpp_reduce_mpi.h .
./diffusivity.F90: $(SRCROOT)atmos_param/diffusivity/diffusivity.F90
	cp $(SRCROOT)atmos_param/diffusivity/diffusivity.F90 .
./ocean_drifters.F90: $(SRCROOT)mom5/ocean_diag/ocean_drifters.F90
	cp $(SRCROOT)mom5/ocean_diag/ocean_drifters.F90 .
./fill_holes.F90: $(SRCROOT)atmos_param/clubb/CLUBB_core/fill_holes.F90
	cp $(SRCROOT)atmos_param/clubb/CLUBB_core/fill_holes.F90 .
./check_nan.F90: $(SRCROOT)atmos_param/strat_cloud/check_nan.F90
	cp $(SRCROOT)atmos_param/strat_cloud/check_nan.F90 .
./sorted_index.F90: $(SRCROOT)atmos_cubed_sphere/tools/sorted_index.F90
	cp $(SRCROOT)atmos_cubed_sphere/tools/sorted_index.F90 .
./gases.F90: $(SRCROOT)atmos_param/cosp/quickbeam/gases.F90
	cp $(SRCROOT)atmos_param/cosp/quickbeam/gases.F90 .
./threadloc.c: $(SRCROOT)shared/mpp/threadloc.c
	cp $(SRCROOT)shared/mpp/threadloc.c .
./shallow_conv.F90: $(SRCROOT)atmos_param/shallow_conv/shallow_conv.F90
	cp $(SRCROOT)atmos_param/shallow_conv/shallow_conv.F90 .
./gradient_c2l.h: /home/wenchang/tigress/CM2.5_tiger2_wy_killtc/exp/../src/shared/mosaic/gradient_c2l.h
	cp /home/wenchang/tigress/CM2.5_tiger2_wy_killtc/exp/../src/shared/mosaic/gradient_c2l.h .
./nh_core.F90: $(SRCROOT)atmos_cubed_sphere/model/null/nh_core.F90
	cp $(SRCROOT)atmos_cubed_sphere/model/null/nh_core.F90 .
./uw_conv.F90: $(SRCROOT)atmos_param/shallow_cu/uw_conv.F90
	cp $(SRCROOT)atmos_param/shallow_cu/uw_conv.F90 .
./corr_matrix_module.F90: $(SRCROOT)atmos_param/clubb/CLUBB_core/corr_matrix_module.F90
	cp $(SRCROOT)atmos_param/clubb/CLUBB_core/corr_matrix_module.F90 .
./calendar.F90: $(SRCROOT)atmos_param/clubb/CLUBB_core/calendar.F90
	cp $(SRCROOT)atmos_param/clubb/CLUBB_core/calendar.F90 .
./moist_processes.F90: $(SRCROOT)atmos_param/moist_processes/moist_processes.F90
	cp $(SRCROOT)atmos_param/moist_processes/moist_processes.F90 .
./fv_eta.F90: $(SRCROOT)atmos_cubed_sphere/tools/fv_eta.F90
	cp $(SRCROOT)atmos_cubed_sphere/tools/fv_eta.F90 .
./rad_output_file.F90: $(SRCROOT)atmos_param/sea_esf_rad/rad_output_file.F90
	cp $(SRCROOT)atmos_param/sea_esf_rad/rad_output_file.F90 .
./mpp_write_2Ddecomp.h: /home/wenchang/tigress/CM2.5_tiger2_wy_killtc/exp/../src/shared/mpp/include/mpp_write_2Ddecomp.h
	cp /home/wenchang/tigress/CM2.5_tiger2_wy_killtc/exp/../src/shared/mpp/include/mpp_write_2Ddecomp.h .
./nsclock.c: $(SRCROOT)shared/mpp/nsclock.c
	cp $(SRCROOT)shared/mpp/nsclock.c .
./pos_definite_module.F90: $(SRCROOT)atmos_param/clubb/CLUBB_core/pos_definite_module.F90
	cp $(SRCROOT)atmos_param/clubb/CLUBB_core/pos_definite_module.F90 .
./extrapolation.F90: $(SRCROOT)atmos_param/clubb/CLUBB_core/extrapolation.F90
	cp $(SRCROOT)atmos_param/clubb/CLUBB_core/extrapolation.F90 .
./anl_erf.F90: $(SRCROOT)atmos_param/clubb/CLUBB_core/anl_erf.F90
	cp $(SRCROOT)atmos_param/clubb/CLUBB_core/anl_erf.F90 .
./ice_albedo.F90: $(SRCROOT)ice_param/ice_albedo.F90
	cp $(SRCROOT)ice_param/ice_albedo.F90 .
./stock_constants.F90: $(SRCROOT)shared/exchange/stock_constants.F90
	cp $(SRCROOT)shared/exchange/stock_constants.F90 .
./external_sst.F90: $(SRCROOT)atmos_cubed_sphere/tools/external_sst.F90
	cp $(SRCROOT)atmos_cubed_sphere/tools/external_sst.F90 .
./mpp_io_connect.inc: /home/wenchang/tigress/CM2.5_tiger2_wy_killtc/exp/../src/shared/mpp/include/mpp_io_connect.inc
	cp /home/wenchang/tigress/CM2.5_tiger2_wy_killtc/exp/../src/shared/mpp/include/mpp_io_connect.inc .
./ocean_increment_tracer.F90: $(SRCROOT)mom5/ocean_param/sources/ocean_increment_tracer.F90
	cp $(SRCROOT)mom5/ocean_param/sources/ocean_increment_tracer.F90 .
./diag_cloud_rad.F90: $(SRCROOT)atmos_param/diag_cloud_rad/diag_cloud_rad.F90
	cp $(SRCROOT)atmos_param/diag_cloud_rad/diag_cloud_rad.F90 .
./cmue_ma.F90: $(SRCROOT)mom5/ocean_param/gotm-4.0/turbulence/cmue_ma.F90
	cp $(SRCROOT)mom5/ocean_param/gotm-4.0/turbulence/cmue_ma.F90 .
SRC = $(SRCROOT)atmos_param/sea_esf_rad/rad_output_file.F90 $(SRCROOT)mom5/ocean_param/neutral/ocean_nphysics.F90 $(SRCROOT)shared/time_interp/time_interp_external.F90 $(SRCROOT)atmos_shared/tracer_driver/atmos_tracer_driver.F90 $(SRCROOT)mom5/ocean_blobs/ocean_blob_diag.F90 $(SRCROOT)mom5/ocean_diag/ocean_velocity_diag.F90 $(SRCROOT)shared/exchange/xgrid.F90 $(SRCROOT)atmos_param/clubb/CLUBB_core/diffusion.F90 $(SRCROOT)mom5/ocean_param/gotm-4.0/turbulence/kbalgebraic.F90 $(SRCROOT)atmos_param/sea_esf_rad/cloud_spec.F90 $(SRCROOT)atmos_param/cosp/cosp_rttov_simulator.F90 $(SRCROOT)mom5/ocean_param/sources/ocean_xlandmix.F90 $(SRCROOT)mom5/ocean_tracers/ocean_residency.F90 $(SRCROOT)atmos_param/clubb/CLUBB_core/parameters_tunable.F90 $(SRCROOT)atmos_param/clubb/CLUBB_core/output_netcdf.F90 $(SRCROOT)atmos_param/moist_processes/moist_processes_utils.F90 $(SRCROOT)mom5/ocean_tracers/ocean_age_tracer.F90 $(SRCROOT)mom5/ocean_param/sources/ocean_overflow_OFP.F90 $(SRCROOT)mom5/ocean_param/gotm-4.0/turbulence/r_ratio.F90 $(SRCROOT)atmos_param/diag_cloud_rad/diag_cloud_rad.F90 $(SRCROOT)shared/mpp/mpp_parameter.F90 $(SRCROOT)atmos_param/cosp/cosp_modis_simulator.F90 $(SRCROOT)shared/mpp/mpp_domains.F90 $(SRCROOT)mom5/ocean_param/gotm-4.0/turbulence/cmue_c.F90 $(SRCROOT)atmos_shared/tracer_driver/aer_ccn_act/ice_nucl.F90 $(SRCROOT)mom5/ocean_param/lateral/ocean_bih_tracer.F90 $(SRCROOT)atmos_param/clubb/CLUBB_core/input_reader.F90 $(SRCROOT)mom5/ocean_param/vertical/ocean_vert_const.F90 $(SRCROOT)mom5/ocean_param/gotm-4.0/turbulence/alpha_mnb.F90 $(SRCROOT)mom5/ocean_core/ocean_workspace.F90 $(SRCROOT)atmos_param/cosp/quickbeam/atmos_lib.F90 $(SRCROOT)atmos_cubed_sphere/tools/sorted_index.F90 $(SRCROOT)mom5/ocean_param/gotm-4.0/turbulence/potentialml.F90 $(SRCROOT)atmos_cubed_sphere/tools/sim_nc_mod.F90 $(SRCROOT)atmos_param/sea_esf_rad/standalone_clouds.F90 $(SRCROOT)land_lad2/vegetation/vegetation.F90 $(SRCROOT)land_lad2/shared/land_tile_diag_buff.F90 $(SRCROOT)atmos_param/moist_processes/moistproc_kernels.F90 $(SRCROOT)coupler/coupler_main.F90 $(SRCROOT)atmos_param/sea_esf_rad/sea_esf_rad.F90 $(SRCROOT)atmos_cubed_sphere/tools/fv_mp_mod.F90 $(SRCROOT)shared/diag_manager/diag_data.F90 $(SRCROOT)atmos_shared/tracer_driver/tropchem/mo_fphoto.F90 $(SRCROOT)atmos_param/sea_esf_rad/lhsw_driver.F90 $(SRCROOT)shared/time_manager/get_cal_time.F90 $(SRCROOT)atmos_param/cosp/quickbeam/format_input.F90 $(SRCROOT)atmos_param/clubb/CLUBB_core/mixing_length.F90 $(SRCROOT)ocean_shared/generic_tracers/generic_CFC.F90 $(SRCROOT)atmos_param/cosp/llnl/llnl_stats.F90 $(SRCROOT)shared/fms/test_fms_io.F90 $(SRCROOT)atmos_param/clubb/CLUBB_core/endian.F90 $(SRCROOT)mom5/ocean_param/neutral/ocean_nphysics_skew.F90 $(SRCROOT)atmos_param/vert_turb_driver/vert_turb_driver.F90 $(SRCROOT)atmos_param/clubb/CLUBB_core/matrix_operations.F90 $(SRCROOT)atmos_param/clubb/CLUBB_core/parameter_indices.F90 $(SRCROOT)atmos_shared/tracer_driver/atmos_tracer_utilities.F90 $(SRCROOT)atmos_param/physics_driver/physics_driver.F90 $(SRCROOT)mom5/ocean_param/gotm-4.0/util/gridinterpol.F90 $(SRCROOT)atmos_param/cosp/quickbeam/gases.F90 $(SRCROOT)atmos_param/strat_cloud/check_nan.F90 $(SRCROOT)atmos_param/clubb/CLUBB_core/pdf_parameter_module.F90 $(SRCROOT)mom5/ocean_param/gotm-4.0/util/adv_center.F90 $(SRCROOT)land_lad2/vegetation/vegn_harvesting.F90 $(SRCROOT)atmos_param/sea_esf_rad/diag_clouds_W.F90 $(SRCROOT)mom5/ocean_param/lateral/ocean_lapcgrid_friction.F90 $(SRCROOT)ice_sis/mask.F90 $(SRCROOT)atmos_cubed_sphere/model/tp_core.F90 $(SRCROOT)atmos_param/cosp/quickbeam/dsd.F90 $(SRCROOT)atmos_shared/tracer_driver/atmos_age_tracer.F90 $(SRCROOT)atmos_param/shallow_cu/conv_utilities.F90 $(SRCROOT)atmos_param/donner_deep/donner_lite_k.F90 $(SRCROOT)atmos_shared/tracer_driver/tropchem/mo_setsox.F90 $(SRCROOT)atmos_param/donner_deep/donner_types.F90 $(SRCROOT)atmos_shared/tracer_driver/atmos_sea_salt.F90 $(SRCROOT)atmos_param/strat_cloud/polysvp.F90 $(SRCROOT)atmos_param/mg_drag/mg_drag.F90 $(SRCROOT)mom5/ocean_core/ocean_velocity.F90 $(SRCROOT)mom5/ocean_core/ocean_topog.F90 $(SRCROOT)mom5/ocean_param/gotm-4.0/turbulence/lengthscaleeq.F90 $(SRCROOT)atmos_param/fsrad/hconst.F90 $(SRCROOT)atmos_param/fsrad/mcm_swtbls.F90 $(SRCROOT)shared/astronomy/astronomy.F90 $(SRCROOT)mom5/ocean_param/lateral/ocean_lap_friction.F90 $(SRCROOT)shared/oda_tools/write_ocean_data.F90 $(SRCROOT)mom5/ocean_param/gotm-4.0/turbulence/cmue_b.F90 $(SRCROOT)land_lad2/vegetation/vegn_cohort.F90 $(SRCROOT)atmos_param/vert_diff/vert_diff.F90 $(SRCROOT)atmos_param/cloud_generator/betaDistribution.F90 $(SRCROOT)atmos_param/cloud_generator/cloud_generator.F90 $(SRCROOT)atmos_param/clubb/CLUBB_core/variables_prognostic_module.F90 $(SRCROOT)atmos_cubed_sphere/model/fv_grid_utils.F90 $(SRCROOT)atmos_param/fsrad/fs_profile.F90 $(SRCROOT)mom5/ocean_param/vertical/ocean_vert_kpp_mom4p0.F90 $(SRCROOT)mom5/ocean_blobs/ocean_blob_dynamic_bottom.F90 $(SRCROOT)shared/mpp/mpp_pset.F90 $(SRCROOT)shared/oda_tools/oda_core_ecda.F90 $(SRCROOT)atmos_cubed_sphere/tools/fv_nudge.F90 $(SRCROOT)atmos_param/donner_deep/donner_deep_miz.F90 $(SRCROOT)atmos_param/clubb/CLUBB_core/stats_rad_zm.F90 $(SRCROOT)mom5/ocean_param/gotm-4.0/turbulence/dissipationeq.F90 $(SRCROOT)atmos_param/strat_cloud/gamma_mg.F90 $(SRCROOT)atmos_param/sea_esf_rad/bulkphys_rad.F90 $(SRCROOT)atmos_param/clubb/CLUBB_core/parameters_radiation.F90 $(SRCROOT)atmos_param/clubb/CLUBB_core/stats_zm.F90 $(SRCROOT)atmos_param/sea_esf_rad/rad_utilities.F90 $(SRCROOT)atmos_param/cosp/cosp.F90 $(SRCROOT)land_lad2/shared/land_tile_io.F90 $(SRCROOT)mom5/ocean_param/lateral/ocean_mixdownslope.F90 $(SRCROOT)atmos_shared/tracer_driver/atmos_ch3i.F90 $(SRCROOT)atmos_param/shallow_cu/deep_conv.F90 $(SRCROOT)atmos_cubed_sphere/model/fv_arrays.F90 $(SRCROOT)land_lad2/shared/sphum.F90 $(SRCROOT)mom5/ocean_blobs/ocean_blob_static_free.F90 $(SRCROOT)mom5/ocean_param/lateral/ocean_bihcst_friction.F90 $(SRCROOT)atmos_param/sea_esf_rad/original_fms_rad.F90 $(SRCROOT)mom5/ocean_diag/ocean_tracer_util.F90 $(SRCROOT)atmos_param/sea_esf_rad/mgrp_prscr_clds.F90 $(SRCROOT)land_lad2/land_model.F90 $(SRCROOT)shared/horiz_interp/horiz_interp.F90 $(SRCROOT)atmos_shared/tracer_driver/tropchem/moz.mods.F90 $(SRCROOT)atmos_param/shallow_cu/conv_plumes_k.F90 $(SRCROOT)atmos_param/cosp/cosp_utils.F90 $(SRCROOT)atmos_param/cloud_zonal/cloud_zonal.F90 $(SRCROOT)atmos_shared/tracer_driver/tropchem/moz.mat.F90 $(SRCROOT)mom5/ocean_tracers/ocean_tracer_advect.F90 $(SRCROOT)mom5/ocean_bgc/ocmip2_co2calc.F90 $(SRCROOT)atmos_param/sea_esf_rad/microphys_cloud.F90 $(SRCROOT)shared/sat_vapor_pres/sat_vapor_pres.F90 $(SRCROOT)atmos_param/strat_cloud/microphysics.F90 $(SRCROOT)mom5/ocean_param/gotm-4.0/turbulence/turbulence.F90 $(SRCROOT)atmos_param/strat_cloud/cldwat2m_micro.F90 $(SRCROOT)atmos_param/clubb/CLUBB_core/stats_zt.F90 $(SRCROOT)atmos_param/strat_cloud/simple_pdf.F90 $(SRCROOT)mom5/ocean_param/gotm-4.0/turbulence/variances.F90 $(SRCROOT)mom5/ocean_param/gotm-4.0/turbulence/fk_craig.F90 $(SRCROOT)land_lad2/shared/table_printer.F90 $(SRCROOT)ocean_shared/generic_tracers/generic_miniBLING.F90 $(SRCROOT)atmos_param/clubb/CLUBB_core/mono_flux_limiter.F90 $(SRCROOT)atmos_param/cosp/actsim/lidar_simulator.F90 $(SRCROOT)atmos_param/cosp/cosp_constants.F90 $(SRCROOT)atmos_param/cosp/icarus-scops-3.7/scops.F $(SRCROOT)atmos_param/cosp/cosp_types.F90 $(SRCROOT)atmos_cubed_sphere/model/sw_core.F90 $(SRCROOT)atmos_cubed_sphere/driver/coupled/atmosphere.F90 $(SRCROOT)atmos_shared/tracer_driver/tropchem/mo_usrrxt.F90 $(SRCROOT)atmos_param/cosp/quickbeam/radar_simulator.F90 $(SRCROOT)atmos_shared/tracer_driver/tropchem/strat_chem_utilities.F90 $(SRCROOT)atmos_param/donner_deep/wet_deposition_0D.F90 $(SRCROOT)shared/coupler/coupler_types.F90 $(SRCROOT)shared/mpp/mpp_memutils.F90 $(SRCROOT)atmos_cubed_sphere/tools/external_ic.F90 $(SRCROOT)atmos_param/cosp/cosp_diagnostics.F90 $(SRCROOT)shared/random_numbers/random_numbers.F90 $(SRCROOT)shared/memutils/memuse.c $(SRCROOT)mom5/ocean_param/gotm-4.0/turbulence/cmue_rf.F90 $(SRCROOT)atmos_param/clubb/CLUBB_core/gmres_wrap.F90 $(SRCROOT)mom5/ocean_param/gotm-4.0/util/util.F90 $(SRCROOT)atmos_param/cosp/cosp_stats.F90 $(SRCROOT)atmos_shared/tracer_driver/tropchem/mo_setinv.F90 $(SRCROOT)shared/mosaic/create_xgrid.c $(SRCROOT)land_lad2/land_constants.F90 $(SRCROOT)shared/axis_utils/axis_utils.F90 $(SRCROOT)land_lad2/river/river_physics.F90 $(SRCROOT)shared/coupler/atmos_ocean_fluxes.F90 $(SRCROOT)atmos_param/strat_cloud/nc_cond.F90 $(SRCROOT)atmos_param/strat_cloud/morrison_gettelman_microp.F90 $(SRCROOT)shared/amip_interp/amip_interp.F90 $(SRCROOT)land_lad2/vegetation/vegn_disturbance.F90 $(SRCROOT)mom5/ocean_core/ocean_velocity_advect.F90 $(SRCROOT)atmos_param/sea_esf_rad/longwave_driver.F90 $(SRCROOT)shared/diag_manager/diag_grid.F90 $(SRCROOT)mom5/ocean_diag/ocean_tracer_diag.F90 $(SRCROOT)mom5/ocean_param/neutral/ocean_nphysics_diff.F90 $(SRCROOT)atmos_cubed_sphere/model/fv_sg.F90 $(SRCROOT)atmos_param/clubb/CLUBB_core/parameters_microphys.F90 $(SRCROOT)atmos_param/fsrad/rdparm.F90 $(SRCROOT)atmos_param/clubb/CLUBB_core/array_index.F90 $(SRCROOT)atmos_cubed_sphere/model/fv_nesting.F90 $(SRCROOT)atmos_param/clubb/CLUBB_core/model_flags.F90 $(SRCROOT)shared/tridiagonal/tridiagonal.F90 $(SRCROOT)atmos_cubed_sphere/tools/fv_ada_nudge.F90 $(SRCROOT)land_lad2/river/river.F90 $(SRCROOT)atmos_param/clubb/CLUBB_core/input_names.F90 $(SRCROOT)atmos_cubed_sphere/model/fv_control.F90 $(SRCROOT)land_lad2/land_tile.F90 $(SRCROOT)atmos_cubed_sphere/tools/fv_timing.F90 $(SRCROOT)shared/horiz_interp/test_horiz_interp.F90 $(SRCROOT)atmos_cubed_sphere/driver/coupled/fv_physics.F90 $(SRCROOT)mom5/ocean_param/lateral/ocean_bihcgrid_friction.F90 $(SRCROOT)mom5/ocean_param/gotm-4.0/util/ode_solvers.F90 $(SRCROOT)atmos_param/clubb/CLUBB_core/extrapolation.F90 $(SRCROOT)mom5/ocean_param/sources/ocean_momentum_source.F90 $(SRCROOT)shared/oda_tools/oda_types.F90 $(SRCROOT)land_lad2/topo_rough/topo_rough.F90 $(SRCROOT)mom5/ocean_bgc/ocean_pert_co2.F90 $(SRCROOT)mom5/ocean_core/ocean_pressure.F90 $(SRCROOT)mom5/ocean_param/sources/ocean_overflow.F90 $(SRCROOT)atmos_param/clubb/CLUBB_core/clubb_core.F90 $(SRCROOT)shared/mpp/test_mpp_pset.F90 $(SRCROOT)atmos_param/clubb/CLUBB_core/stats_subs.F90 $(SRCROOT)atmos_param/clubb/CLUBB_3D_var.F90 $(SRCROOT)mom5/ocean_param/gotm-4.0/turbulence/cmue_d.F90 $(SRCROOT)atmos_shared/atmos_nudge/atmos_nudge.F90 $(SRCROOT)atmos_param/moist_processes/moist_processes.F90 $(SRCROOT)mom5/ocean_param/gotm-4.0/turbulence/tkealgebraic.F90 $(SRCROOT)mom5/ocean_tracers/ocean_tpm_util.F90 $(SRCROOT)atmos_param/sea_esf_rad/shortwave_driver.F90 $(SRCROOT)ice_param/ocean_rough.F90 $(SRCROOT)atmos_shared/tracer_driver/tropchem/mo_photo.F90 $(SRCROOT)shared/mpp/mpp_data.F90 $(SRCROOT)atmos_param/fsrad/rad_diag.F90 $(SRCROOT)mom5/ocean_param/gotm-4.0/turbulence/compute_rist.F90 $(SRCROOT)atmos_param/cosp/icarus-scops-3.7/isccp_cloud_types.F $(SRCROOT)mom5/ocean_param/neutral/ocean_nphysics_new.F90 $(SRCROOT)atmos_param/fsrad/mcm_swnew.F90 $(SRCROOT)mom5/ocean_bgc/ocean_bgc_restore.F90 $(SRCROOT)atmos_param/strat_cloud/strat_cloud.F90 $(SRCROOT)atmos_shared/tracer_driver/tropchem/mo_chem_utls.F90 $(SRCROOT)mom5/ocean_core/ocean_domains.F90 $(SRCROOT)atmos_param/clubb/CLUBB_core/stats_type.F90 $(SRCROOT)mom5/ocean_bgc/ocmip2_abiotic.F90 $(SRCROOT)atmos_param/clubb/CLUBB_core/advance_wp2_wp3_module.F90 $(SRCROOT)shared/mpp/threadloc.c $(SRCROOT)atmos_param/donner_deep/donner_cape_k.F90 $(SRCROOT)atmos_shared/tracer_driver/atmos_convection_tracer.F90 $(SRCROOT)mom5/ocean_core/ocean_util.F90 $(SRCROOT)mom5/ocean_param/vertical/ocean_vert_pp.F90 $(SRCROOT)atmos_param/cloud_rad/cloud_rad.F90 $(SRCROOT)ice_param/ice_albedo.F90 $(SRCROOT)atmos_param/fsrad/longwave.F90 $(SRCROOT)mom5/ocean_tracers/ocean_tracer.F90 $(SRCROOT)land_lad2/shared/nf_utils/nfc.F90 $(SRCROOT)ice_sis/ice_thm.F90 $(SRCROOT)atmos_cubed_sphere/model/fv_update_phys.F90 $(SRCROOT)mom5/ocean_tracers/ocean_residency_ml.F90 $(SRCROOT)mom5/ocean_core/ocean_obc_barotrop.F90 $(SRCROOT)ice_sis/ice_spec.F90 $(SRCROOT)atmos_param/sea_esf_rad/optical_path.F90 $(SRCROOT)atmos_param/grey_radiation/grey_radiation.F90 $(SRCROOT)atmos_param/clubb/CLUBB_driver_SCM.F90 $(SRCROOT)mom5/ocean_param/vertical/ocean_convect.F90 $(SRCROOT)atmos_shared/tracer_driver/tropchem/m_tracname.F90 $(SRCROOT)atmos_param/cg_drag/cg_drag.F90 $(SRCROOT)atmos_param/cosp/MISR_simulator/MISR_simulator.F $(SRCROOT)shared/diag_manager/diag_table.F90 $(SRCROOT)shared/time_manager/time_manager.F90 $(SRCROOT)mom5/ocean_param/gotm-4.0/turbulence/cmue_sg.F90 $(SRCROOT)atmos_param/sea_esf_rad/sealw99.F90 $(SRCROOT)shared/exchange/stock_constants.F90 $(SRCROOT)atmos_shared/vert_advection/vert_advection.F90 $(SRCROOT)atmos_param/strat_cloud/strat_cloud_utilities.F90 $(SRCROOT)mom5/ocean_param/gotm-4.0/turbulence/cmue_a.F90 $(SRCROOT)atmos_param/clouds/clouds.F90 $(SRCROOT)mom5/ocean_param/gotm-4.0/util/convert_fluxes.F90 $(SRCROOT)atmos_param/donner_deep/donner_rad_k.F90 $(SRCROOT)shared/horiz_interp/horiz_interp_spherical.F90 $(SRCROOT)shared/topography/topography.F90 $(SRCROOT)atmos_param/clubb/CLUBB_core/pos_definite_module.F90 $(SRCROOT)atmos_param/monin_obukhov/monin_obukhov_kernel.F90 $(SRCROOT)atmos_shared/tracer_driver/stratchem/strat_chem_model.F90 $(SRCROOT)atmos_param/cosp/llnl/pf_to_mr.F $(SRCROOT)atmos_param/strat_cloud/strat_netcdf.F90 $(SRCROOT)atmos_cubed_sphere/tools/fv_diagnostics.F90 $(SRCROOT)land_lad2/shared/nf_utils/nf_utils.F90 $(SRCROOT)atmos_cubed_sphere/model/null/nh_core.F90 $(SRCROOT)mom5/ocean_param/sources/ocean_rivermix.F90 $(SRCROOT)ice_sis/ice_type.F90 $(SRCROOT)ocean_shared/generic_tracers/generic_tracer.F90 $(SRCROOT)atmos_param/betts_miller/betts_miller.F90 $(SRCROOT)atmos_param/donner_deep/donner_deep_k.F90 $(SRCROOT)mom5/ocean_core/ocean_bbc.F90 $(SRCROOT)atmos_shared/tracer_driver/tropchem/tropchem_driver.F90 $(SRCROOT)atmos_shared/tracer_driver/tropchem/mo_hook.F90 $(SRCROOT)atmos_param/clubb/CLUBB_core/mt95.f90 $(SRCROOT)atmos_param/donner_deep/fms_donner.F90 $(SRCROOT)atmos_param/rh_clouds/rh_clouds.F90 $(SRCROOT)land_lad2/vegetation/vegn_static_override.F90 $(SRCROOT)atmos_shared/interpolator/interpolator.F90 $(SRCROOT)atmos_param/lscale_cond/lscale_cond.F90 $(SRCROOT)atmos_param/cloud_obs/cloud_obs.F90 $(SRCROOT)atmos_param/clubb/CLUBB_core/advance_xp2_xpyp_module.F90 $(SRCROOT)atmos_param/donner_deep/nonfms_donner.F90 $(SRCROOT)land_lad2/vegetation/vegn_photosynthesis.F90 $(SRCROOT)atmos_param/clubb/CLUBB_core/clip_explicit.F90 $(SRCROOT)atmos_shared/tracer_driver/atmos_sulfate.F90 $(SRCROOT)shared/mpp/test_mpp_domains.F90 $(SRCROOT)land_lad2/river/river_type.F90 $(SRCROOT)atmos_param/clubb/CLUBB_core/sigma_sqd_w_module.F90 $(SRCROOT)atmos_param/donner_deep/cumulus_closure_k.F90 $(SRCROOT)atmos_shared/tracer_driver/atmos_soa.F90 $(SRCROOT)atmos_param/clubb/CLUBB_core/Skw_module.F90 $(SRCROOT)atmos_param/sea_esf_rad/specified_clouds_W.F90 $(SRCROOT)mom5/ocean_param/gotm-4.0/turbulence/ispralength.F90 $(SRCROOT)shared/oda_tools/xbt_drop_rate_adjust.f90 $(SRCROOT)shared/fms/fms_io.F90 $(SRCROOT)atmos_param/clubb/CLUBB_core/interpolation.F90 $(SRCROOT)atmos_param/clubb/CLUBB_core/output_grads.F90 $(SRCROOT)mom5/ocean_core/ocean_model.F90 $(SRCROOT)ice_sis/ice_grid.F90 $(SRCROOT)atmos_param/clubb/CLUBB_core/hyper_diffusion_4th_ord.F90 $(SRCROOT)shared/diag_manager/diag_util.F90 $(SRCROOT)mom5/ocean_param/gotm-4.0/turbulence/kpp.F90 $(SRCROOT)ocean_shared/generic_tracers/generic_COBALT.F90 $(SRCROOT)atmos_param/clubb/CLUBB_core/calendar.F90 $(SRCROOT)mom5/ocean_param/sources/ocean_shortwave_gfdl.F90 $(SRCROOT)shared/time_interp/time_interp.F90 $(SRCROOT)land_lad2/snow/snow.F90 $(SRCROOT)shared/memutils/memutils.F90 $(SRCROOT)mom5/ocean_param/gotm-4.0/util/time.F90 $(SRCROOT)shared/topography/gaussian_topog.F90 $(SRCROOT)atmos_shared/tracer_driver/aer_ccn_act/aerosol_params.F90 $(SRCROOT)shared/diag_manager/diag_axis.F90 $(SRCROOT)atmos_shared/tracer_driver/stratchem/strat_chem_driver.F90 $(SRCROOT)mom5/ocean_param/gotm-4.0/util/lagrange.F90 $(SRCROOT)mom5/ocean_param/gotm-4.0/util/eqstate.F90 $(SRCROOT)atmos_param/clubb/CLUBB_core/variables_diagnostic_module.F90 $(SRCROOT)mom5/ocean_param/vertical/ocean_form_drag.F90 $(SRCROOT)mom5/ocean_param/neutral/ocean_nphysicsC.F90 $(SRCROOT)mom5/ocean_param/vertical/ocean_vert_chen.F90 $(SRCROOT)atmos_param/fsrad/shortwave.F90 $(SRCROOT)shared/mosaic/interp.c $(SRCROOT)mom5/ocean_tracers/ocean_tempsalt.F90 $(SRCROOT)atmos_cubed_sphere/tools/fv_io.F90 $(SRCROOT)shared/horiz_interp/horiz_interp_conserve.F90 $(SRCROOT)mom5/ocean_param/gotm-4.0/turbulence/kbeq.F90 $(SRCROOT)ocean_shared/generic_tracers/generic_ERGOM.F90 $(SRCROOT)atmos_cubed_sphere/tools/external_sst.F90 $(SRCROOT)shared/platform/platform.F90 $(SRCROOT)atmos_param/clubb/CLUBB_core/surface_varnce_module.F90 $(SRCROOT)land_lad2/land_data.F90 $(SRCROOT)mom5/ocean_param/sources/ocean_xlandinsert.F90 $(SRCROOT)mom5/ocean_core/ocean_thickness.F90 $(SRCROOT)mom5/ocean_core/oda_driver.F90 $(SRCROOT)mom5/ocean_param/vertical/ocean_vert_kpp_test.F90 $(SRCROOT)atmos_cubed_sphere/tools/test_cases.F90 $(SRCROOT)shared/drifters/drifters.F90 $(SRCROOT)mom5/ocean_core/ocean_density.F90 $(SRCROOT)atmos_param/strat_cloud/micro_mg.F90 $(SRCROOT)atmos_param/diag_cloud/diag_cloud.F90 $(SRCROOT)atmos_shared/tracer_driver/tropchem/mo_exp_slv.F90 $(SRCROOT)land_lad2/vegetation/vegn_radiation.F90 $(SRCROOT)atmos_param/clubb/CLUBB_core/anl_erf.F90 $(SRCROOT)shared/sat_vapor_pres/sat_vapor_pres_k.F90 $(SRCROOT)atmos_cubed_sphere/model/fv_tracer2d.F90 $(SRCROOT)atmos_param/monin_obukhov/monin_obukhov.F90 $(SRCROOT)shared/mosaic/gradient.F90 $(SRCROOT)land_lad2/lake/lake_tile.F90 $(SRCROOT)mom5/ocean_diag/ocean_drifters.F90 $(SRCROOT)atmos_cubed_sphere/model/dyn_core.F90 $(SRCROOT)shared/fft/fft.F90 $(SRCROOT)atmos_param/moist_conv/moist_conv.F90 $(SRCROOT)ocean_shared/generic_tracers/generic_tracer_utils.F90 $(SRCROOT)atmos_param/vert_diff_driver/vert_diff_driver.F90 $(SRCROOT)mom5/ocean_param/neutral/ocean_nphysics_tensor.F90 $(SRCROOT)shared/constants/constants.F90 $(SRCROOT)atmos_param/betts_miller/bm_massflux.F90 $(SRCROOT)mom5/ocean_param/lateral/ocean_lapgen_friction.F90 $(SRCROOT)atmos_param/clubb/advance_sclrm_Nd_module.F90 $(SRCROOT)shared/column_diagnostics/column_diagnostics.F90 $(SRCROOT)atmos_param/clubb/CLUBB_core/clip_semi_implicit.F90 $(SRCROOT)mom5/ocean_param/gotm-4.0/turbulence/internal_wave.F90 $(SRCROOT)mom5/ocean_param/sources/ocean_sponges_tracer.F90 $(SRCROOT)shared/mpp/nsclock.c $(SRCROOT)land_lad2/shared/land_io.F90 $(SRCROOT)mom5/ocean_param/sources/ocean_increment_eta.F90 $(SRCROOT)mom5/ocean_param/gotm-4.0/util/diff_center.F90 $(SRCROOT)atmos_param/lin_cloud_microphys/lin_cloud_microphys.F90 $(SRCROOT)atmos_param/cosp/cosp_lidar.F90 $(SRCROOT)atmos_param/shallow_cu/uw_conv.F90 $(SRCROOT)mom5/ocean_wave/ocean_wave_type.F90 $(SRCROOT)atmos_param/sea_esf_rad/zetac_clouds_W.F90 $(SRCROOT)atmos_param/sea_esf_rad/radiation_diag.F90 $(SRCROOT)land_lad2/shared/land_numerics.F90 $(SRCROOT)atmos_cubed_sphere/tools/read_climate_nudge_data.F90 $(SRCROOT)atmos_shared/tracer_driver/tropchem/moz.subs.F90 $(SRCROOT)atmos_param/clubb/CLUBB_core/error_code.F90 $(SRCROOT)land_lad2/soil/soil.F90 $(SRCROOT)shared/tracer_manager/tracer_manager.F90 $(SRCROOT)mom5/ocean_wave/ocean_wave.F90 $(SRCROOT)atmos_param/clubb/MG_microp_3D.F90 $(SRCROOT)atmos_param/sea_esf_rad/isccp_clouds.F90 $(SRCROOT)atmos_cubed_sphere/model/a2b_edge.F90 $(SRCROOT)mom5/ocean_param/sources/ocean_shortwave.F90 $(SRCROOT)atmos_param/cosp/actsim/lmd_ipsl_stats.F90 $(SRCROOT)atmos_shared/tracer_driver/tropchem/mo_rodas_slv.F90 $(SRCROOT)atmos_param/sea_esf_rad/aerosol.F90 $(SRCROOT)mom5/ocean_bgc/ocmip2_he.F90 $(SRCROOT)ocean_shared/generic_tracers/generic_BLING.F90 $(SRCROOT)mom5/ocean_param/sources/ocean_overexchange.F90 $(SRCROOT)atmos_param/clubb/CLUBB_core/csr_matrix_class_3array.F90 $(SRCROOT)shared/mpp/mpp_utilities.F90 $(SRCROOT)land_lad2/shared/nf_utils/nfu.F90 $(SRCROOT)atmos_param/clubb/CLUBB_core/advance_helper_module.F90 $(SRCROOT)mom5/ocean_param/sources/ocean_sponges_eta.F90 $(SRCROOT)atmos_param/clubb/CLUBB_core/fill_holes.F90 $(SRCROOT)atmos_param/clubb/CLUBB_core/numerical_check.F90 $(SRCROOT)mom5/ocean_param/gotm-4.0/turbulence/compute_cpsi3.F90 $(SRCROOT)atmos_param/fsrad/mcm_sw_driver.F90 $(SRCROOT)shared/diag_manager/diag_output.F90 $(SRCROOT)atmos_param/clubb/CLUBB_core/hydrostatic_module.F90 $(SRCROOT)atmos_param/cosp/MODIS_simulator/modis_simulator.F90 $(SRCROOT)atmos_param/cosp/quickbeam/array_lib.F90 $(SRCROOT)atmos_param/sea_esf_rad/uw_clouds_W.F90 $(SRCROOT)atmos_param/clubb/CLUBB_core/mean_adv.F90 $(SRCROOT)atmos_cubed_sphere/model/fv_dynamics.F90 $(SRCROOT)atmos_cubed_sphere/model/fv_mapz.F90 $(SRCROOT)mom5/ocean_core/ocean_sbc.F90 $(SRCROOT)atmos_param/cosp/cosp_simulator.F90 $(SRCROOT)shared/exchange/test_xgrid.F90 $(SRCROOT)atmos_coupled/atmos_model.F90 $(SRCROOT)atmos_param/stable_bl_turb/stable_bl_turb.F90 $(SRCROOT)land_lad2/glacier/glacier.F90 $(SRCROOT)atmos_param/donner_deep/donner_deep.F90 $(SRCROOT)atmos_cubed_sphere/tools/fv_restart.F90 $(SRCROOT)atmos_param/sea_esf_rad/strat_clouds_W.F90 $(SRCROOT)atmos_param/fsrad/co2int.F90 $(SRCROOT)atmos_shared/tracer_driver/aer_ccn_act/aer_ccn_act_k.F90 $(SRCROOT)mom5/ocean_param/lateral/ocean_bih_friction.F90 $(SRCROOT)atmos_param/clubb/CLUBB_core/grid_class.F90 $(SRCROOT)atmos_param/sea_esf_rad/esfsw_parameters.F90 $(SRCROOT)atmos_param/strat_cloud/aerosol_cloud.F90 $(SRCROOT)mom5/ocean_tracers/ocean_frazil.F90 $(SRCROOT)mom5/ocean_param/lateral/ocean_bihgen_friction.F90 $(SRCROOT)mom5/ocean_param/gotm-4.0/turbulence/gotm_lib_version.F90 $(SRCROOT)mom5/ocean_diag/ocean_diagnostics.F90 $(SRCROOT)atmos_shared/tracer_driver/tropchem/mo_fastjx.F90 $(SRCROOT)shared/drifters/drifters_io.F90 $(SRCROOT)shared/mosaic/grid.F90 $(SRCROOT)mom5/ocean_param/gotm-4.0/turbulence/cmue_ma.F90 $(SRCROOT)mom5/ocean_param/lateral/ocean_submesoscale.F90 $(SRCROOT)atmos_param/shallow_cu/conv_closures.F90 $(SRCROOT)atmos_param/strat_cloud/strat_cloud_legacy.F90 $(SRCROOT)land_lad2/vegetation/vegn_tile.F90 $(SRCROOT)atmos_shared/tracer_driver/aer_ccn_act/aer_in_act.F90 $(SRCROOT)atmos_cubed_sphere/tools/fv_grid_tools.F90 $(SRCROOT)atmos_cubed_sphere/model/fv_fill.F90 $(SRCROOT)land_lad2/vegetation/vegn_data.F90 $(SRCROOT)shared/field_manager/field_manager.F90 $(SRCROOT)land_lad2/canopy_air/canopy_air.F90 $(SRCROOT)mom5/ocean_param/gotm-4.0/turbulence/q2over2eq.F90 $(SRCROOT)mom5/ocean_tracers/transport_matrix.F90 $(SRCROOT)atmos_param/clubb/CLUBB_core/stats_LH_zt.F90 $(SRCROOT)mom5/ocean_param/neutral/ocean_nphysics_util.F90 $(SRCROOT)atmos_param/clubb/CLUBB_core/gmres_cache.F90 $(SRCROOT)atmos_shared/tracer_driver/atmos_carbon_aerosol.F90 $(SRCROOT)atmos_param/ras/ras.F90 $(SRCROOT)atmos_param/cosp/quickbeam/load_hydrometeor_classes.F90 $(SRCROOT)atmos_param/fsrad/co2_data.F90 $(SRCROOT)mom5/ocean_core/ocean_advection_velocity.F90 $(SRCROOT)shared/horiz_interp/horiz_interp_type.F90 $(SRCROOT)atmos_cubed_sphere/tools/fv_surf_map.F90 $(SRCROOT)shared/mpp/test_mpp.F90 $(SRCROOT)mom5/ocean_tracers/ocean_residency_range.F90 $(SRCROOT)atmos_param/cosp/llnl/prec_scops.F $(SRCROOT)mom5/ocean_param/neutral/ocean_nphysics_util_new.F90 $(SRCROOT)atmos_param/cosp/quickbeam/zeff.F90 $(SRCROOT)atmos_param/moist_processes/detr_ice_num.F90 $(SRCROOT)atmos_param/clubb/CLUBB_core/pdf_closure_module.F90 $(SRCROOT)atmos_shared/tracer_driver/tropchem/mo_chemdr.F90 $(SRCROOT)shared/mpp/affinity.c $(SRCROOT)shared/mosaic/gradient_c2l.c $(SRCROOT)mom5/ocean_blobs/ocean_blob_dynamic_free.F90 $(SRCROOT)mom5/ocean_param/lateral/ocean_sigma_transport.F90 $(SRCROOT)land_lad2/canopy_air/cana_tile.F90 $(SRCROOT)land_lad2/lake/lake.F90 $(SRCROOT)atmos_param/edt/edt.F90 $(SRCROOT)atmos_param/shallow_cu/conv_utilities_k.F90 $(SRCROOT)mom5/ocean_param/vertical/ocean_vert_gotm.F90 $(SRCROOT)land_lad2/soil/uptake.F90 $(SRCROOT)mom5/ocean_tracers/ocean_residency_integrand.F90 $(SRCROOT)mom5/ocean_param/sources/ocean_increment_tracer.F90 $(SRCROOT)land_lad2/glacier/glac_tile.F90 $(SRCROOT)mom5/ocean_bgc/ocmip2_biotic.F90 $(SRCROOT)shared/drifters/cloud_interpolator.F90 $(SRCROOT)atmos_param/clubb/CLUBB_core/saturation.F90 $(SRCROOT)atmos_param/clubb/alt_cloud.F90 $(SRCROOT)ocean_shared/generic_tracers/FMS_coupler_util.F90 $(SRCROOT)ice_param/ocean_albedo.F90 $(SRCROOT)mom5/ocean_param/gotm-4.0/turbulence/tkeeq.F90 $(SRCROOT)atmos_param/clubb/CLUBB_core/variables_radiation_module.F90 $(SRCROOT)atmos_shared/tracer_driver/atmos_dust.F90 $(SRCROOT)mom5/ocean_param/gotm-4.0/turbulence/genericeq.F90 $(SRCROOT)shared/mosaic/mosaic.F90 $(SRCROOT)atmos_shared/tracer_driver/tropchem/mo_jpl.F90 $(SRCROOT)shared/horiz_interp/horiz_interp_bicubic.F90 $(SRCROOT)mom5/ocean_param/lateral/ocean_lap_tracer.F90 $(SRCROOT)atmos_param/cosp/quickbeam/mrgrnk.F90 $(SRCROOT)mom5/ocean_param/vertical/ocean_vert_kpp_mom4p1.F90 $(SRCROOT)atmos_param/donner_deep/donner_utilities_k.F90 $(SRCROOT)shared/oda_tools/oda_core.F90 $(SRCROOT)shared/diag_manager/diag_manager.F90 $(SRCROOT)land_lad2/shared/land_tile_diag.F90 $(SRCROOT)atmos_param/strat_cloud/mg_const.F90 $(SRCROOT)atmos_shared/tracer_driver/tropchem/mo_read_sim_chm.F90 $(SRCROOT)mom5/ocean_param/lateral/ocean_lapcst_friction.F90 $(SRCROOT)atmos_param/clubb/CLUBB_core/stat_file_module.F90 $(SRCROOT)ice_sis/ice_bergs.F90 $(SRCROOT)ice_sis/ice_model.F90 $(SRCROOT)shared/drifters/quicksort.F90 $(SRCROOT)atmos_param/clubb/CLUBB_core/sponge_layer_damping.F90 $(SRCROOT)mom5/ocean_param/neutral/ocean_nphysics_flux.F90 $(SRCROOT)coupler/surface_flux.F90 $(SRCROOT)atmos_param/donner_deep/donner_meso_k.F90 $(SRCROOT)atmos_param/cosp/cosp_isccp_simulator.F90 $(SRCROOT)atmos_param/cosp/quickbeam/math_lib.F90 $(SRCROOT)mom5/ocean_param/sources/ocean_sponges_velocity.F90 $(SRCROOT)atmos_param/diffusivity/diffusivity.F90 $(SRCROOT)mom5/ocean_tracers/ocean_residency_meta.F90 $(SRCROOT)mom5/ocean_core/ocean_operators.F90 $(SRCROOT)ocean_shared/generic_tracers/FMS_ocmip2_co2calc.F90 $(SRCROOT)shared/mpp/mpp_io.F90 $(SRCROOT)mom5/ocean_core/ocean_coriolis.F90 $(SRCROOT)land_lad2/vegetation/vegn_cohort_io.F90 $(SRCROOT)atmos_param/clubb/CLUBB_core/lapack_wrap.F90 $(SRCROOT)atmos_cubed_sphere/tools/init_hydro.F90 $(SRCROOT)shared/field_manager/fm_util.F90 $(SRCROOT)atmos_param/sea_esf_rad/ozone.F90 $(SRCROOT)atmos_param/betts_miller/bm_omp.F90 $(SRCROOT)ice_sis/ice_dyn.F90 $(SRCROOT)atmos_param/clubb/CLUBB_core/file_functions.F90 $(SRCROOT)mom5/ocean_diag/ocean_adv_vel_diag.F90 $(SRCROOT)shared/mosaic/read_mosaic.c $(SRCROOT)atmos_param/clubb/CLUBB_core/advance_xm_wpxp_module.F90 $(SRCROOT)mom5/ocean_param/neutral/ocean_nphysicsB.F90 $(SRCROOT)atmos_param/my25_turb/my25_turb.F90 $(SRCROOT)ocean_shared/generic_tracers/generic_TOPAZ.F90 $(SRCROOT)atmos_param/radiation_driver/radiation_driver.F90 $(SRCROOT)ice_sis/rot.F90 $(SRCROOT)mom5/ocean_param/vertical/ocean_vert_util.F90 $(SRCROOT)atmos_param/shallow_conv/shallow_conv.F90 $(SRCROOT)land_lad2/transitions/transitions.F90 $(SRCROOT)atmos_param/donner_deep/donner_cloud_model_k.F90 $(SRCROOT)atmos_param/fsrad/mcm_lw.F90 $(SRCROOT)mom5/ocean_param/gotm-4.0/turbulence/production.F90 $(SRCROOT)atmos_param/sea_esf_rad/donner_deep_clouds_W.F90 $(SRCROOT)shared/fft/fft99.F90 $(SRCROOT)mom5/ocean_param/vertical/ocean_vert_tidal_test.F90 $(SRCROOT)atmos_param/donner_deep/donner_lscloud_k.F90 $(SRCROOT)shared/coupler/ensemble_manager.F90 $(SRCROOT)land_lad2/shared/land_utils.F90 $(SRCROOT)atmos_param/diag_integral/diag_integral.F90 $(SRCROOT)atmos_param/cosp/cosp_driver.F90 $(SRCROOT)land_lad2/soil/soil_tile.F90 $(SRCROOT)atmos_param/cu_mo_trans/cu_mo_trans.F90 $(SRCROOT)atmos_cubed_sphere/tools/fv_climate_nudge.F90 $(SRCROOT)atmos_param/clubb/CLUBB_core/corr_matrix_module.F90 $(SRCROOT)atmos_param/tke_turb/tke_turb.F90 $(SRCROOT)mom5/ocean_tracers/ocean_passive.F90 $(SRCROOT)atmos_shared/tracer_driver/tropchem/mo_imp_slv.F90 $(SRCROOT)shared/data_override/data_override.F90 $(SRCROOT)mom5/ocean_core/ocean_obc.F90 $(SRCROOT)land_lad2/shared/land_debug.F90 $(SRCROOT)mom5/ocean_param/sources/ocean_increment_velocity.F90 $(SRCROOT)mom5/ocean_param/vertical/ocean_vert_mix.F90 $(SRCROOT)shared/mosaic/mosaic_util.c $(SRCROOT)mom5/ocean_core/ocean_grids.F90 $(SRCROOT)atmos_param/strat_cloud/rotstayn_klein_mp.F90 $(SRCROOT)atmos_param/sea_esf_rad/longwave_fluxes.F90 $(SRCROOT)atmos_param/clubb/CLUBB_core/parameters_model.F90 $(SRCROOT)mom5/ocean_param/gotm-4.0/util/diff_face.F90 $(SRCROOT)mom5/ocean_param/sources/ocean_shortwave_csiro.F90 $(SRCROOT)land_lad2/vegetation/vegn_dynamics.F90 $(SRCROOT)mom5/ocean_bgc/ocean_generic_tracer.F90 $(SRCROOT)mom5/ocean_bgc/ocmip2_cfc.F90 $(SRCROOT)atmos_shared/tracer_driver/tropchem/mo_chemini.F90 $(SRCROOT)mom5/ocean_core/ocean_parameters.F90 $(SRCROOT)atmos_param/cosp/quickbeam/radar_simulator_types.F90 $(SRCROOT)mom5/ocean_param/sources/ocean_shortwave_jerlov.F90 $(SRCROOT)atmos_param/sea_esf_rad/aerosolrad_package.F90 $(SRCROOT)atmos_param/sea_esf_rad/cloudrad_package.F90 $(SRCROOT)shared/fms/fms.F90 $(SRCROOT)mom5/ocean_param/gotm-4.0/util/tridiagonal_gotm.F90 $(SRCROOT)shared/mpp/test_mpp_io.F90 $(SRCROOT)mom5/ocean_blobs/ocean_blob_util.F90 $(SRCROOT)atmos_param/cosp/cosp_io.F90 $(SRCROOT)mom5/ocean_param/neutral/ocean_nphysicsA.F90 $(SRCROOT)atmos_param/clubb/CLUBB_core/constants_clubb.F90 $(SRCROOT)shared/drifters/drifters_comm.F90 $(SRCROOT)atmos_param/cosp/icarus-scops-3.7/icarus.F90 $(SRCROOT)mom5/ocean_param/sources/ocean_riverspread.F90 $(SRCROOT)atmos_param/sea_esf_rad/microphys_rad.F90 $(SRCROOT)atmos_param/clubb/CLUBB_core/stats_variables.F90 $(SRCROOT)atmos_param/sea_esf_rad/rh_based_clouds.F90 $(SRCROOT)atmos_param/sea_esf_rad/cloudrad_diagnostics.F90 $(SRCROOT)atmos_param/clubb/CLUBB_core/stats_sfc.F90 $(SRCROOT)atmos_cubed_sphere/model/boundary.F90 $(SRCROOT)atmos_param/sea_esf_rad/longwave_tables.F90 $(SRCROOT)atmos_param/damping_driver/damping_driver.F90 $(SRCROOT)atmos_param/clubb/CLUBB_core/advance_windm_edsclrm_module.F90 $(SRCROOT)mom5/ocean_core/ocean_barotropic.F90 $(SRCROOT)atmos_param/sea_esf_rad/esfsw_driver.F90 $(SRCROOT)atmos_param/sea_esf_rad/gas_tf.F90 $(SRCROOT)atmos_param/dry_adj/dry_adj.F90 $(SRCROOT)atmos_param/topo_drag/topo_drag.F90 $(SRCROOT)mom5/ocean_bgc/ocean_po4_pre.F90 $(SRCROOT)mom5/ocean_tracers/ocean_tpm.F90 $(SRCROOT)atmos_param/cosp/cosp_misr_simulator.F90 $(SRCROOT)atmos_shared/tracer_driver/atmos_sulfur_hex.F90 $(SRCROOT)mom5/ocean_blobs/ocean_blob.F90 $(SRCROOT)mom5/ocean_param/gotm-4.0/turbulence/algebraiclength.F90 $(SRCROOT)shared/drifters/drifters_input.F90 $(SRCROOT)atmos_param/entrain/entrain.F90 $(SRCROOT)atmos_param/clubb/CLUBB_core/T_in_K_module.F90 $(SRCROOT)mom5/ocean_blobs/ocean_blob_static_bottom.F90 $(SRCROOT)shared/horiz_interp/horiz_interp_bilinear.F90 $(SRCROOT)atmos_param/sea_esf_rad/radiative_gases.F90 $(SRCROOT)atmos_param/sea_esf_rad/longwave_clouds.F90 $(SRCROOT)atmos_param/sea_esf_rad/longwave_params.F90 $(SRCROOT)mom5/ocean_param/gotm-4.0/turbulence/epsbalgebraic.F90 $(SRCROOT)atmos_param/clubb/CLUBB_core/clubb_precision.F90 $(SRCROOT)coupler/flux_exchange.F90 $(SRCROOT)mom5/ocean_param/vertical/ocean_vert_tidal.F90 $(SRCROOT)atmos_cubed_sphere/tools/fv_eta.F90 $(SRCROOT)land_lad2/snow/snow_tile.F90 $(SRCROOT)atmos_shared/tracer_driver/atmos_radon.F90 $(SRCROOT)atmos_shared/tracer_driver/aer_ccn_act/aer_ccn_act.F90 $(SRCROOT)mom5/ocean_core/ocean_types.F90 $(SRCROOT)land_lad2/shared/land_tile_diag_sel.F90 $(SRCROOT)atmos_param/cosp/llnl/cosp_radar.F90 $(SRCROOT)atmos_param/fsrad/fsrad.F90 $(SRCROOT)atmos_param/clubb/CLUBB_core/stats_rad_zt.F90 $(SRCROOT)mom5/ocean_bgc/ocean_ibgc.F90 $(SRCROOT)atmos_shared/tracer_driver/atmos_co2.F90 $(SRCROOT)atmos_param/sea_esf_rad/lw_gases_stdtf.F90 $(SRCROOT)shared/mpp/mpp.F90 $(SRCROOT)shared/random_numbers/MersenneTwister.F90 $(SRCROOT)shared/drifters/drifters_core.F90 $(SRCROOT)atmos_param/cosp/quickbeam/optics_lib.F90 $(SRCROOT)atmos_param/cosp/quickbeam/load_mie_table.F90 /home/wenchang/tigress/CM2.5_tiger2_wy_killtc/exp/../src/shared/mpp/include/mpp_domains_reduce.inc /home/wenchang/tigress/CM2.5_tiger2_wy_killtc/exp/../src/shared/mpp/include/mpp_error_s_s.h /home/wenchang/tigress/CM2.5_tiger2_wy_killtc/exp/../src/land_lad2/shared/nf_utils/getput_compressed.inc /home/wenchang/tigress/CM2.5_tiger2_wy_killtc/exp/../src/shared/mpp/include/mpp_domains_misc.inc /home/wenchang/tigress/CM2.5_tiger2_wy_killtc/exp/../src/shared/mpp/include/mpp_transmit_sma.h /home/wenchang/tigress/CM2.5_tiger2_wy_killtc/exp/../src/atmos_param/donner_deep/donner_nml.h /home/wenchang/tigress/CM2.5_tiger2_wy_killtc/exp/../src/land_lad2/vegetation/read_remap_cohort_data_new.inc /home/wenchang/tigress/CM2.5_tiger2_wy_killtc/exp/../src/shared/mpp/include/mpp_do_get_boundary.h /home/wenchang/tigress/CM2.5_tiger2_wy_killtc/exp/../src/shared/mpp/include/mpp_util_nocomm.inc /home/wenchang/tigress/CM2.5_tiger2_wy_killtc/exp/../src/land_lad2/shared/nf_utils/getput.inc /home/wenchang/tigress/CM2.5_tiger2_wy_killtc/exp/../src/shared/mpp/include/mpp_read_compressed.h /home/wenchang/tigress/CM2.5_tiger2_wy_killtc/exp/../src/shared/mpp/include/mpp_io_read.inc /home/wenchang/tigress/CM2.5_tiger2_wy_killtc/exp/../src/atmos_param/monin_obukhov/monin_obukhov_interfaces.h /home/wenchang/tigress/CM2.5_tiger2_wy_killtc/exp/../src/shared/mpp/include/mpp_domains_util.inc /home/wenchang/tigress/CM2.5_tiger2_wy_killtc/exp/../src/shared/mpp/include/mpp_comm_sma.inc /home/wenchang/tigress/CM2.5_tiger2_wy_killtc/exp/../src/shared/mpp/include/mpp_global_reduce.h /home/wenchang/tigress/CM2.5_tiger2_wy_killtc/exp/../src/shared/mpp/include/mpp_update_domains2D.h /home/wenchang/tigress/CM2.5_tiger2_wy_killtc/exp/../src/shared/mpp/include/mpp_util_sma.inc /home/wenchang/tigress/CM2.5_tiger2_wy_killtc/exp/../src/atmos_param/cosp/icarus-scops-3.7/congvec.H /home/wenchang/tigress/CM2.5_tiger2_wy_killtc/exp/../src/land_lad2/canopy_air/../shared/debug.inc /home/wenchang/tigress/CM2.5_tiger2_wy_killtc/exp/../src/shared/mpp/include/mpp_comm_nocomm.inc /home/wenchang/tigress/CM2.5_tiger2_wy_killtc/exp/../src/shared/mpp/include/mpp_chksum.h /home/wenchang/tigress/CM2.5_tiger2_wy_killtc/exp/../src/atmos_param/clubb/CLUBB_core/recl.inc /home/wenchang/tigress/CM2.5_tiger2_wy_killtc/exp/../src/shared/mpp/include/mpp_io_write.inc /home/wenchang/tigress/CM2.5_tiger2_wy_killtc/exp/../src/shared/fms/read_data_4d.inc /home/wenchang/tigress/CM2.5_tiger2_wy_killtc/exp/../src/shared/mpp/include/mpp_do_update_nonblock.h /home/wenchang/tigress/CM2.5_tiger2_wy_killtc/exp/../src/shared/mpp/include/mpp_comm_mpi.inc /home/wenchang/tigress/CM2.5_tiger2_wy_killtc/exp/../src/shared/mosaic/interp.h /home/wenchang/tigress/CM2.5_tiger2_wy_killtc/exp/../src/shared/mpp/include/mpp_update_nest_domains.h /home/wenchang/tigress/CM2.5_tiger2_wy_killtc/exp/../src/shared/mpp/include/mpp_data_nocomm.inc /home/wenchang/tigress/CM2.5_tiger2_wy_killtc/exp/../src/shared/mpp/include/mpp_data_mpi.inc /home/wenchang/tigress/CM2.5_tiger2_wy_killtc/exp/../src/shared/mpp/include/mpp_transmit_nocomm.h /home/wenchang/tigress/CM2.5_tiger2_wy_killtc/exp/../src/shared/mpp/include/mpp_chksum_scalar.h /home/wenchang/tigress/CM2.5_tiger2_wy_killtc/exp/../src/mom5/ocean_param/gotm-4.0/include/version.h /home/wenchang/tigress/CM2.5_tiger2_wy_killtc/exp/../src/shared/mpp/include/mpp_define_nest_domains.inc /home/wenchang/tigress/CM2.5_tiger2_wy_killtc/exp/../src/shared/mpp/include/mpp_domains_comm.inc /home/wenchang/tigress/CM2.5_tiger2_wy_killtc/exp/../src/shared/mpp/include/mpp_data_sma.inc /home/wenchang/tigress/CM2.5_tiger2_wy_killtc/exp/../src/shared/field_manager/parse.inc /home/wenchang/tigress/CM2.5_tiger2_wy_killtc/exp/../src/shared/mosaic/create_xgrid.h /home/wenchang/tigress/CM2.5_tiger2_wy_killtc/exp/../src/shared/mpp/include/mpp_gather.h /home/wenchang/tigress/CM2.5_tiger2_wy_killtc/exp/../src/mom5/ocean_core/ocean_memory.h /home/wenchang/tigress/CM2.5_tiger2_wy_killtc/exp/../src/shared/mpp/include/mpp_io_misc.inc /home/wenchang/tigress/CM2.5_tiger2_wy_killtc/exp/../src/shared/mpp/include/mpp_global_field.h /home/wenchang/tigress/CM2.5_tiger2_wy_killtc/exp/../src/shared/mpp/include/mpp_domains_define.inc /home/wenchang/tigress/CM2.5_tiger2_wy_killtc/exp/../src/shared/fms/read_data_2d.inc /home/wenchang/tigress/CM2.5_tiger2_wy_killtc/exp/../src/shared/mpp/include/mpp_do_global_field.h /home/wenchang/tigress/CM2.5_tiger2_wy_killtc/exp/../src/shared/mpp/include/system_clock.h /home/wenchang/tigress/CM2.5_tiger2_wy_killtc/exp/../src/shared/mpp/include/mpp_group_update.h /home/wenchang/tigress/CM2.5_tiger2_wy_killtc/exp/../src/shared/mpp/include/mpp_transmit.inc /home/wenchang/tigress/CM2.5_tiger2_wy_killtc/exp/../src/land_lad2/vegetation/vegn_cohort_io.inc /home/wenchang/tigress/CM2.5_tiger2_wy_killtc/exp/../src/shared/include/fms_platform.h /home/wenchang/tigress/CM2.5_tiger2_wy_killtc/exp/../src/shared/mpp/include/mpp_do_redistribute.h /home/wenchang/tigress/CM2.5_tiger2_wy_killtc/exp/../src/shared/mpp/include/mpp_do_updateV.h /home/wenchang/tigress/CM2.5_tiger2_wy_killtc/exp/../src/shared/fms/read_data_3d.inc /home/wenchang/tigress/CM2.5_tiger2_wy_killtc/exp/../src/shared/mpp/include/mpp_error_a_s.h /home/wenchang/tigress/CM2.5_tiger2_wy_killtc/exp/../src/shared/mpp/include/mpp_sum_nocomm.h /home/wenchang/tigress/CM2.5_tiger2_wy_killtc/exp/../src/shared/mpp/include/mpp_scatter.h /home/wenchang/tigress/CM2.5_tiger2_wy_killtc/exp/../src/shared/mosaic/mosaic_util.h /home/wenchang/tigress/CM2.5_tiger2_wy_killtc/exp/../src/shared/mpp/include/mpp_reduce_mpi.h /home/wenchang/tigress/CM2.5_tiger2_wy_killtc/exp/../src/shared/mpp/include/mpp_do_check.h /home/wenchang/tigress/CM2.5_tiger2_wy_killtc/exp/../src/shared/mpp/include/mpp_write_compressed.h /home/wenchang/tigress/CM2.5_tiger2_wy_killtc/exp/../src/shared/drifters/fms_switches.h /home/wenchang/tigress/CM2.5_tiger2_wy_killtc/exp/../src/shared/mpp/include/mpp_update_domains2D_nonblock.h /home/wenchang/tigress/CM2.5_tiger2_wy_killtc/exp/../src/atmos_param/donner_deep/donner_types.h /home/wenchang/tigress/CM2.5_tiger2_wy_killtc/exp/../src/shared/mosaic/constant.h /home/wenchang/tigress/CM2.5_tiger2_wy_killtc/exp/../src/shared/mpp/include/mpp_error_a_a.h /home/wenchang/tigress/CM2.5_tiger2_wy_killtc/exp/../src/shared/drifters/drifters_set_field.h /home/wenchang/tigress/CM2.5_tiger2_wy_killtc/exp/../src/shared/mpp/include/mpp_do_updateV_nonblock.h /home/wenchang/tigress/CM2.5_tiger2_wy_killtc/exp/../src/shared/mpp/include/mpp_read_2Ddecomp.h /home/wenchang/tigress/CM2.5_tiger2_wy_killtc/exp/../src/shared/mpp/include/mpp_global_sum_tl.h /home/wenchang/tigress/CM2.5_tiger2_wy_killtc/exp/../src/land_lad2/shared/debug.inc /home/wenchang/tigress/CM2.5_tiger2_wy_killtc/exp/../src/shared/mpp/include/mpp_sum.inc /home/wenchang/tigress/CM2.5_tiger2_wy_killtc/exp/../src/shared/mpp/include/mpp_global_sum.h /home/wenchang/tigress/CM2.5_tiger2_wy_killtc/exp/../src/atmos_param/strat_cloud/strat_nml.h /home/wenchang/tigress/CM2.5_tiger2_wy_killtc/exp/../src/shared/mpp/include/mpp_reduce_nocomm.h /home/wenchang/tigress/CM2.5_tiger2_wy_killtc/exp/../src/shared/mpp/include/mpp_chksum_int.h /home/wenchang/tigress/CM2.5_tiger2_wy_killtc/exp/../src/shared/mpp/include/mpp_io_util.inc /home/wenchang/tigress/CM2.5_tiger2_wy_killtc/exp/../src/shared/drifters/drifters_compute_k.h /home/wenchang/tigress/CM2.5_tiger2_wy_killtc/exp/../src/shared/mpp/include/mpp_sum_mpi.h /home/wenchang/tigress/CM2.5_tiger2_wy_killtc/exp/../src/atmos_param/cosp/cosp_defs.H /home/wenchang/tigress/CM2.5_tiger2_wy_killtc/exp/../src/shared/mpp/include/mpp_util_mpi.inc /home/wenchang/tigress/CM2.5_tiger2_wy_killtc/exp/../src/shared/mpp/include/mpp_get_boundary.h /home/wenchang/tigress/CM2.5_tiger2_wy_killtc/exp/../src/shared/drifters/drifters_push.h /home/wenchang/tigress/CM2.5_tiger2_wy_killtc/exp/../src/shared/mpp/include/mpp_error_s_a.h /home/wenchang/tigress/CM2.5_tiger2_wy_killtc/exp/../src/shared/mosaic/gradient_c2l.h /home/wenchang/tigress/CM2.5_tiger2_wy_killtc/exp/../src/shared/mpp/include/mpp_comm.inc /home/wenchang/tigress/CM2.5_tiger2_wy_killtc/exp/../src/shared/mpp/include/mpp_do_update_nest.h /home/wenchang/tigress/CM2.5_tiger2_wy_killtc/exp/../src/shared/mpp/include/mpp_do_checkV.h /home/wenchang/tigress/CM2.5_tiger2_wy_killtc/exp/../src/shared/mpp/include/mpp_reduce_sma.h /home/wenchang/tigress/CM2.5_tiger2_wy_killtc/exp/../src/land_lad2/vegetation/read_remap_cohort_data.inc /home/wenchang/tigress/CM2.5_tiger2_wy_killtc/exp/../src/shared/mpp/include/mpp_do_update.h /home/wenchang/tigress/CM2.5_tiger2_wy_killtc/exp/../src/mom5/ocean_param/gotm-4.0/include/cppdefs.h /home/wenchang/tigress/CM2.5_tiger2_wy_killtc/exp/../src/shared/mpp/include/mpp_sum_sma.h /home/wenchang/tigress/CM2.5_tiger2_wy_killtc/exp/../src/shared/mpp/include/mpp_write_2Ddecomp.h /home/wenchang/tigress/CM2.5_tiger2_wy_killtc/exp/../src/shared/fms/write_data.inc /home/wenchang/tigress/CM2.5_tiger2_wy_killtc/exp/../src/shared/mosaic/read_mosaic.h /home/wenchang/tigress/CM2.5_tiger2_wy_killtc/exp/../src/shared/mpp/include/mpp_util.inc /home/wenchang/tigress/CM2.5_tiger2_wy_killtc/exp/../src/shared/mpp/include/mpp_io_connect.inc /home/wenchang/tigress/CM2.5_tiger2_wy_killtc/exp/../src/shared/mpp/include/mpp_write.h /home/wenchang/tigress/CM2.5_tiger2_wy_killtc/exp/../src/shared/mpp/include/mpp_transmit_mpi.h
OBJ = rad_output_file.o ocean_nphysics.o time_interp_external.o atmos_tracer_driver.o ocean_blob_diag.o ocean_velocity_diag.o xgrid.o diffusion.o kbalgebraic.o cloud_spec.o cosp_rttov_simulator.o ocean_xlandmix.o ocean_residency.o parameters_tunable.o output_netcdf.o moist_processes_utils.o ocean_age_tracer.o ocean_overflow_OFP.o r_ratio.o diag_cloud_rad.o mpp_parameter.o cosp_modis_simulator.o mpp_domains.o cmue_c.o ice_nucl.o ocean_bih_tracer.o input_reader.o ocean_vert_const.o alpha_mnb.o ocean_workspace.o atmos_lib.o sorted_index.o potentialml.o sim_nc_mod.o standalone_clouds.o vegetation.o land_tile_diag_buff.o moistproc_kernels.o coupler_main.o sea_esf_rad.o fv_mp_mod.o diag_data.o mo_fphoto.o lhsw_driver.o get_cal_time.o format_input.o mixing_length.o generic_CFC.o llnl_stats.o test_fms_io.o endian.o ocean_nphysics_skew.o vert_turb_driver.o matrix_operations.o parameter_indices.o atmos_tracer_utilities.o physics_driver.o gridinterpol.o gases.o check_nan.o pdf_parameter_module.o adv_center.o vegn_harvesting.o diag_clouds_W.o ocean_lapcgrid_friction.o mask.o tp_core.o dsd.o atmos_age_tracer.o conv_utilities.o donner_lite_k.o mo_setsox.o donner_types.o atmos_sea_salt.o polysvp.o mg_drag.o ocean_velocity.o ocean_topog.o lengthscaleeq.o hconst.o mcm_swtbls.o astronomy.o ocean_lap_friction.o write_ocean_data.o cmue_b.o vegn_cohort.o vert_diff.o betaDistribution.o cloud_generator.o variables_prognostic_module.o fv_grid_utils.o fs_profile.o ocean_vert_kpp_mom4p0.o ocean_blob_dynamic_bottom.o mpp_pset.o oda_core_ecda.o fv_nudge.o donner_deep_miz.o stats_rad_zm.o dissipationeq.o gamma_mg.o bulkphys_rad.o parameters_radiation.o stats_zm.o rad_utilities.o cosp.o land_tile_io.o ocean_mixdownslope.o atmos_ch3i.o deep_conv.o fv_arrays.o sphum.o ocean_blob_static_free.o ocean_bihcst_friction.o original_fms_rad.o ocean_tracer_util.o mgrp_prscr_clds.o land_model.o horiz_interp.o moz.mods.o conv_plumes_k.o cosp_utils.o cloud_zonal.o moz.mat.o ocean_tracer_advect.o ocmip2_co2calc.o microphys_cloud.o sat_vapor_pres.o microphysics.o turbulence.o cldwat2m_micro.o stats_zt.o simple_pdf.o variances.o fk_craig.o table_printer.o generic_miniBLING.o mono_flux_limiter.o lidar_simulator.o cosp_constants.o scops.o cosp_types.o sw_core.o atmosphere.o mo_usrrxt.o radar_simulator.o strat_chem_utilities.o wet_deposition_0D.o coupler_types.o mpp_memutils.o external_ic.o cosp_diagnostics.o random_numbers.o memuse.o cmue_rf.o gmres_wrap.o util.o cosp_stats.o mo_setinv.o create_xgrid.o land_constants.o axis_utils.o river_physics.o atmos_ocean_fluxes.o nc_cond.o morrison_gettelman_microp.o amip_interp.o vegn_disturbance.o ocean_velocity_advect.o longwave_driver.o diag_grid.o ocean_tracer_diag.o ocean_nphysics_diff.o fv_sg.o parameters_microphys.o rdparm.o array_index.o fv_nesting.o model_flags.o tridiagonal.o fv_ada_nudge.o river.o input_names.o fv_control.o land_tile.o fv_timing.o test_horiz_interp.o fv_physics.o ocean_bihcgrid_friction.o ode_solvers.o extrapolation.o ocean_momentum_source.o oda_types.o topo_rough.o ocean_pert_co2.o ocean_pressure.o ocean_overflow.o clubb_core.o test_mpp_pset.o stats_subs.o CLUBB_3D_var.o cmue_d.o atmos_nudge.o moist_processes.o tkealgebraic.o ocean_tpm_util.o shortwave_driver.o ocean_rough.o mo_photo.o mpp_data.o rad_diag.o compute_rist.o isccp_cloud_types.o ocean_nphysics_new.o mcm_swnew.o ocean_bgc_restore.o strat_cloud.o mo_chem_utls.o ocean_domains.o stats_type.o ocmip2_abiotic.o advance_wp2_wp3_module.o threadloc.o donner_cape_k.o atmos_convection_tracer.o ocean_util.o ocean_vert_pp.o cloud_rad.o ice_albedo.o longwave.o ocean_tracer.o nfc.o ice_thm.o fv_update_phys.o ocean_residency_ml.o ocean_obc_barotrop.o ice_spec.o optical_path.o grey_radiation.o CLUBB_driver_SCM.o ocean_convect.o m_tracname.o cg_drag.o MISR_simulator.o diag_table.o time_manager.o cmue_sg.o sealw99.o stock_constants.o vert_advection.o strat_cloud_utilities.o cmue_a.o clouds.o convert_fluxes.o donner_rad_k.o horiz_interp_spherical.o topography.o pos_definite_module.o monin_obukhov_kernel.o strat_chem_model.o pf_to_mr.o strat_netcdf.o fv_diagnostics.o nf_utils.o nh_core.o ocean_rivermix.o ice_type.o generic_tracer.o betts_miller.o donner_deep_k.o ocean_bbc.o tropchem_driver.o mo_hook.o mt95.o fms_donner.o rh_clouds.o vegn_static_override.o interpolator.o lscale_cond.o cloud_obs.o advance_xp2_xpyp_module.o nonfms_donner.o vegn_photosynthesis.o clip_explicit.o atmos_sulfate.o test_mpp_domains.o river_type.o sigma_sqd_w_module.o cumulus_closure_k.o atmos_soa.o Skw_module.o specified_clouds_W.o ispralength.o xbt_drop_rate_adjust.o fms_io.o interpolation.o output_grads.o ocean_model.o ice_grid.o hyper_diffusion_4th_ord.o diag_util.o kpp.o generic_COBALT.o calendar.o ocean_shortwave_gfdl.o time_interp.o snow.o memutils.o time.o gaussian_topog.o aerosol_params.o diag_axis.o strat_chem_driver.o lagrange.o eqstate.o variables_diagnostic_module.o ocean_form_drag.o ocean_nphysicsC.o ocean_vert_chen.o shortwave.o interp.o ocean_tempsalt.o fv_io.o horiz_interp_conserve.o kbeq.o generic_ERGOM.o external_sst.o platform.o surface_varnce_module.o land_data.o ocean_xlandinsert.o ocean_thickness.o oda_driver.o ocean_vert_kpp_test.o test_cases.o drifters.o ocean_density.o micro_mg.o diag_cloud.o mo_exp_slv.o vegn_radiation.o anl_erf.o sat_vapor_pres_k.o fv_tracer2d.o monin_obukhov.o gradient.o lake_tile.o ocean_drifters.o dyn_core.o fft.o moist_conv.o generic_tracer_utils.o vert_diff_driver.o ocean_nphysics_tensor.o constants.o bm_massflux.o ocean_lapgen_friction.o advance_sclrm_Nd_module.o column_diagnostics.o clip_semi_implicit.o internal_wave.o ocean_sponges_tracer.o nsclock.o land_io.o ocean_increment_eta.o diff_center.o lin_cloud_microphys.o cosp_lidar.o uw_conv.o ocean_wave_type.o zetac_clouds_W.o radiation_diag.o land_numerics.o read_climate_nudge_data.o moz.subs.o error_code.o soil.o tracer_manager.o ocean_wave.o MG_microp_3D.o isccp_clouds.o a2b_edge.o ocean_shortwave.o lmd_ipsl_stats.o mo_rodas_slv.o aerosol.o ocmip2_he.o generic_BLING.o ocean_overexchange.o csr_matrix_class_3array.o mpp_utilities.o nfu.o advance_helper_module.o ocean_sponges_eta.o fill_holes.o numerical_check.o compute_cpsi3.o mcm_sw_driver.o diag_output.o hydrostatic_module.o modis_simulator.o array_lib.o uw_clouds_W.o mean_adv.o fv_dynamics.o fv_mapz.o ocean_sbc.o cosp_simulator.o test_xgrid.o atmos_model.o stable_bl_turb.o glacier.o donner_deep.o fv_restart.o strat_clouds_W.o co2int.o aer_ccn_act_k.o ocean_bih_friction.o grid_class.o esfsw_parameters.o aerosol_cloud.o ocean_frazil.o ocean_bihgen_friction.o gotm_lib_version.o ocean_diagnostics.o mo_fastjx.o drifters_io.o grid.o cmue_ma.o ocean_submesoscale.o conv_closures.o strat_cloud_legacy.o vegn_tile.o aer_in_act.o fv_grid_tools.o fv_fill.o vegn_data.o field_manager.o canopy_air.o q2over2eq.o transport_matrix.o stats_LH_zt.o ocean_nphysics_util.o gmres_cache.o atmos_carbon_aerosol.o ras.o load_hydrometeor_classes.o co2_data.o ocean_advection_velocity.o horiz_interp_type.o fv_surf_map.o test_mpp.o ocean_residency_range.o prec_scops.o ocean_nphysics_util_new.o zeff.o detr_ice_num.o pdf_closure_module.o mo_chemdr.o affinity.o gradient_c2l.o ocean_blob_dynamic_free.o ocean_sigma_transport.o cana_tile.o lake.o edt.o conv_utilities_k.o ocean_vert_gotm.o uptake.o ocean_residency_integrand.o ocean_increment_tracer.o glac_tile.o ocmip2_biotic.o cloud_interpolator.o saturation.o alt_cloud.o FMS_coupler_util.o ocean_albedo.o tkeeq.o variables_radiation_module.o atmos_dust.o genericeq.o mosaic.o mo_jpl.o horiz_interp_bicubic.o ocean_lap_tracer.o mrgrnk.o ocean_vert_kpp_mom4p1.o donner_utilities_k.o oda_core.o diag_manager.o land_tile_diag.o mg_const.o mo_read_sim_chm.o ocean_lapcst_friction.o stat_file_module.o ice_bergs.o ice_model.o quicksort.o sponge_layer_damping.o ocean_nphysics_flux.o surface_flux.o donner_meso_k.o cosp_isccp_simulator.o math_lib.o ocean_sponges_velocity.o diffusivity.o ocean_residency_meta.o ocean_operators.o FMS_ocmip2_co2calc.o mpp_io.o ocean_coriolis.o vegn_cohort_io.o lapack_wrap.o init_hydro.o fm_util.o ozone.o bm_omp.o ice_dyn.o file_functions.o ocean_adv_vel_diag.o read_mosaic.o advance_xm_wpxp_module.o ocean_nphysicsB.o my25_turb.o generic_TOPAZ.o radiation_driver.o rot.o ocean_vert_util.o shallow_conv.o transitions.o donner_cloud_model_k.o mcm_lw.o production.o donner_deep_clouds_W.o fft99.o ocean_vert_tidal_test.o donner_lscloud_k.o ensemble_manager.o land_utils.o diag_integral.o cosp_driver.o soil_tile.o cu_mo_trans.o fv_climate_nudge.o corr_matrix_module.o tke_turb.o ocean_passive.o mo_imp_slv.o data_override.o ocean_obc.o land_debug.o ocean_increment_velocity.o ocean_vert_mix.o mosaic_util.o ocean_grids.o rotstayn_klein_mp.o longwave_fluxes.o parameters_model.o diff_face.o ocean_shortwave_csiro.o vegn_dynamics.o ocean_generic_tracer.o ocmip2_cfc.o mo_chemini.o ocean_parameters.o radar_simulator_types.o ocean_shortwave_jerlov.o aerosolrad_package.o cloudrad_package.o fms.o tridiagonal_gotm.o test_mpp_io.o ocean_blob_util.o cosp_io.o ocean_nphysicsA.o constants_clubb.o drifters_comm.o icarus.o ocean_riverspread.o microphys_rad.o stats_variables.o rh_based_clouds.o cloudrad_diagnostics.o stats_sfc.o boundary.o longwave_tables.o damping_driver.o advance_windm_edsclrm_module.o ocean_barotropic.o esfsw_driver.o gas_tf.o dry_adj.o topo_drag.o ocean_po4_pre.o ocean_tpm.o cosp_misr_simulator.o atmos_sulfur_hex.o ocean_blob.o algebraiclength.o drifters_input.o entrain.o T_in_K_module.o ocean_blob_static_bottom.o horiz_interp_bilinear.o radiative_gases.o longwave_clouds.o longwave_params.o epsbalgebraic.o clubb_precision.o flux_exchange.o ocean_vert_tidal.o fv_eta.o snow_tile.o atmos_radon.o aer_ccn_act.o ocean_types.o land_tile_diag_sel.o cosp_radar.o fsrad.o stats_rad_zt.o ocean_ibgc.o atmos_co2.o lw_gases_stdtf.o mpp.o MersenneTwister.o drifters_core.o optics_lib.o load_mie_table.o
OFF = $(SRCROOT)mom5/ocean_param/gotm-4.0/util/adv_center.F90 $(SRCROOT)mom5/ocean_param/gotm-4.0/turbulence/kbalgebraic.F90 $(SRCROOT)atmos_param/clubb/CLUBB_core/stats_zt.F90 $(SRCROOT)mom5/ocean_bgc/ocmip2_abiotic.F90 $(SRCROOT)atmos_param/cosp/quickbeam/radar_simulator_types.F90 $(SRCROOT)atmos_shared/tracer_driver/tropchem/m_tracname.F90 $(SRCROOT)atmos_param/fsrad/co2_data.F90 $(SRCROOT)mom5/ocean_param/gotm-4.0/turbulence/epsbalgebraic.F90 $(SRCROOT)atmos_param/strat_cloud/cldwat2m_micro.F90 $(SRCROOT)atmos_shared/tracer_driver/atmos_dust.F90 $(SRCROOT)mom5/ocean_tracers/ocean_tpm_util.F90 $(SRCROOT)mom5/ocean_param/sources/ocean_shortwave_gfdl.F90 /home/wenchang/tigress/CM2.5_tiger2_wy_killtc/exp/../src/atmos_param/monin_obukhov/monin_obukhov_interfaces.h $(SRCROOT)shared/mpp/mpp_pset.F90 $(SRCROOT)shared/drifters/drifters_io.F90 $(SRCROOT)atmos_param/clubb/CLUBB_core/array_index.F90 $(SRCROOT)mom5/ocean_param/gotm-4.0/turbulence/algebraiclength.F90 $(SRCROOT)mom5/ocean_param/gotm-4.0/util/time.F90 $(SRCROOT)atmos_param/sea_esf_rad/standalone_clouds.F90 $(SRCROOT)atmos_cubed_sphere/tools/fv_grid_tools.F90 $(SRCROOT)shared/exchange/test_xgrid.F90 $(SRCROOT)mom5/ocean_core/ocean_grids.F90 $(SRCROOT)atmos_param/betts_miller/bm_massflux.F90 /home/wenchang/tigress/CM2.5_tiger2_wy_killtc/exp/../src/shared/mpp/include/mpp_comm_mpi.inc $(SRCROOT)atmos_param/clubb/CLUBB_core/hydrostatic_module.F90 $(SRCROOT)atmos_shared/tracer_driver/atmos_age_tracer.F90 $(SRCROOT)atmos_param/cosp/quickbeam/optics_lib.F90 $(SRCROOT)atmos_cubed_sphere/model/fv_dynamics.F90 $(SRCROOT)shared/drifters/cloud_interpolator.F90 $(SRCROOT)atmos_param/clubb/CLUBB_core/variables_radiation_module.F90 $(SRCROOT)atmos_param/clubb/CLUBB_core/interpolation.F90 $(SRCROOT)atmos_param/betts_miller/betts_miller.F90 /home/wenchang/tigress/CM2.5_tiger2_wy_killtc/exp/../src/shared/mpp/include/mpp_group_update.h $(SRCROOT)atmos_cubed_sphere/tools/fv_timing.F90 $(SRCROOT)shared/coupler/atmos_ocean_fluxes.F90 $(SRCROOT)shared/mpp/mpp_data.F90 $(SRCROOT)land_lad2/canopy_air/canopy_air.F90 $(SRCROOT)mom5/ocean_param/neutral/ocean_nphysics_skew.F90 $(SRCROOT)shared/amip_interp/amip_interp.F90 $(SRCROOT)atmos_coupled/atmos_model.F90 $(SRCROOT)atmos_param/clubb/CLUBB_core/sigma_sqd_w_module.F90 $(SRCROOT)atmos_param/clubb/CLUBB_core/Skw_module.F90 $(SRCROOT)shared/mosaic/gradient.F90 $(SRCROOT)atmos_param/cosp/cosp_modis_simulator.F90 $(SRCROOT)mom5/ocean_core/ocean_advection_velocity.F90 $(SRCROOT)atmos_cubed_sphere/tools/fv_surf_map.F90 $(SRCROOT)atmos_param/clubb/CLUBB_core/advance_xm_wpxp_module.F90 /home/wenchang/tigress/CM2.5_tiger2_wy_killtc/exp/../src/shared/drifters/fms_switches.h $(SRCROOT)mom5/ocean_param/gotm-4.0/util/tridiagonal_gotm.F90 $(SRCROOT)mom5/ocean_core/ocean_topog.F90 $(SRCROOT)mom5/ocean_param/sources/ocean_increment_eta.F90 $(SRCROOT)ocean_shared/generic_tracers/FMS_coupler_util.F90 $(SRCROOT)shared/drifters/drifters_comm.F90 $(SRCROOT)mom5/ocean_blobs/ocean_blob_static_bottom.F90 $(SRCROOT)atmos_param/monin_obukhov/monin_obukhov_kernel.F90 $(SRCROOT)mom5/ocean_param/gotm-4.0/turbulence/lengthscaleeq.F90 /home/wenchang/tigress/CM2.5_tiger2_wy_killtc/exp/../src/shared/mpp/include/mpp_global_sum_tl.h /home/wenchang/tigress/CM2.5_tiger2_wy_killtc/exp/../src/land_lad2/shared/debug.inc /home/wenchang/tigress/CM2.5_tiger2_wy_killtc/exp/../src/shared/mpp/include/mpp_chksum_int.h $(SRCROOT)atmos_param/cosp/quickbeam/format_input.F90 /home/wenchang/tigress/CM2.5_tiger2_wy_killtc/exp/../src/shared/drifters/drifters_compute_k.h $(SRCROOT)shared/mpp/test_mpp_domains.F90 $(SRCROOT)mom5/ocean_param/gotm-4.0/turbulence/alpha_mnb.F90 $(SRCROOT)atmos_shared/tracer_driver/atmos_sea_salt.F90 $(SRCROOT)land_lad2/vegetation/vegetation.F90 /home/wenchang/tigress/CM2.5_tiger2_wy_killtc/exp/../src/shared/mpp/include/mpp_reduce_sma.h $(SRCROOT)atmos_param/cosp/actsim/lidar_simulator.F90 /home/wenchang/tigress/CM2.5_tiger2_wy_killtc/exp/../src/land_lad2/vegetation/read_remap_cohort_data.inc $(SRCROOT)mom5/ocean_core/ocean_pressure.F90 $(SRCROOT)ice_sis/ice_dyn.F90 $(SRCROOT)mom5/ocean_param/gotm-4.0/util/convert_fluxes.F90 /home/wenchang/tigress/CM2.5_tiger2_wy_killtc/exp/../src/shared/mpp/include/mpp_do_update.h $(SRCROOT)mom5/ocean_param/neutral/ocean_nphysics_new.F90 $(SRCROOT)mom5/ocean_param/vertical/ocean_vert_tidal.F90 $(SRCROOT)mom5/ocean_param/vertical/ocean_vert_gotm.F90 $(SRCROOT)ice_sis/ice_grid.F90 $(SRCROOT)atmos_param/clubb/CLUBB_core/mt95.f90 $(SRCROOT)shared/fft/fft99.F90 $(SRCROOT)atmos_shared/tracer_driver/tropchem/strat_chem_utilities.F90 $(SRCROOT)shared/oda_tools/oda_core_ecda.F90 $(SRCROOT)ocean_shared/generic_tracers/FMS_ocmip2_co2calc.F90 $(SRCROOT)land_lad2/vegetation/vegn_radiation.F90 $(SRCROOT)mom5/ocean_param/neutral/ocean_nphysicsA.F90 $(SRCROOT)mom5/ocean_param/neutral/ocean_nphysics_util_new.F90 $(SRCROOT)atmos_param/sea_esf_rad/shortwave_driver.F90 /home/wenchang/tigress/CM2.5_tiger2_wy_killtc/exp/../src/land_lad2/shared/nf_utils/getput_compressed.inc /home/wenchang/tigress/CM2.5_tiger2_wy_killtc/exp/../src/shared/mpp/include/mpp_domains_misc.inc $(SRCROOT)coupler/flux_exchange.F90 $(SRCROOT)atmos_param/cosp/cosp_types.F90 $(SRCROOT)atmos_param/cosp/quickbeam/dsd.F90 $(SRCROOT)mom5/ocean_core/ocean_thickness.F90 $(SRCROOT)atmos_param/fsrad/longwave.F90 $(SRCROOT)shared/field_manager/fm_util.F90 $(SRCROOT)atmos_param/strat_cloud/mg_const.F90 $(SRCROOT)atmos_param/clubb/CLUBB_core/parameters_radiation.F90 $(SRCROOT)shared/coupler/ensemble_manager.F90 $(SRCROOT)atmos_param/monin_obukhov/monin_obukhov.F90 $(SRCROOT)mom5/ocean_param/gotm-4.0/turbulence/gotm_lib_version.F90 $(SRCROOT)atmos_param/diag_integral/diag_integral.F90 $(SRCROOT)mom5/ocean_param/gotm-4.0/turbulence/dissipationeq.F90 $(SRCROOT)atmos_param/my25_turb/my25_turb.F90 $(SRCROOT)atmos_param/fsrad/mcm_lw.F90 /home/wenchang/tigress/CM2.5_tiger2_wy_killtc/exp/../src/atmos_param/clubb/CLUBB_core/recl.inc $(SRCROOT)mom5/ocean_param/lateral/ocean_submesoscale.F90 $(SRCROOT)atmos_param/clubb/CLUBB_driver_SCM.F90 $(SRCROOT)atmos_param/moist_conv/moist_conv.F90 $(SRCROOT)atmos_param/clubb/CLUBB_core/input_reader.F90 $(SRCROOT)land_lad2/shared/land_tile_diag.F90 $(SRCROOT)mom5/ocean_param/gotm-4.0/turbulence/turbulence.F90 $(SRCROOT)shared/mpp/test_mpp_io.F90 $(SRCROOT)mom5/ocean_param/gotm-4.0/util/ode_solvers.F90 $(SRCROOT)shared/diag_manager/diag_grid.F90 $(SRCROOT)atmos_param/sea_esf_rad/lhsw_driver.F90 $(SRCROOT)atmos_param/clubb/CLUBB_core/advance_helper_module.F90 $(SRCROOT)land_lad2/transitions/transitions.F90 $(SRCROOT)mom5/ocean_diag/ocean_tracer_diag.F90 $(SRCROOT)ice_sis/ice_type.F90 $(SRCROOT)atmos_param/cosp/cosp_constants.F90 $(SRCROOT)atmos_param/physics_driver/physics_driver.F90 $(SRCROOT)atmos_cubed_sphere/model/a2b_edge.F90 $(SRCROOT)atmos_param/clubb/CLUBB_core/clubb_core.F90 $(SRCROOT)land_lad2/vegetation/vegn_cohort.F90 $(SRCROOT)atmos_param/clubb/CLUBB_core/variables_prognostic_module.F90 /home/wenchang/tigress/CM2.5_tiger2_wy_killtc/exp/../src/shared/mpp/include/mpp_do_global_field.h $(SRCROOT)atmos_shared/tracer_driver/tropchem/moz.subs.F90 $(SRCROOT)shared/constants/constants.F90 $(SRCROOT)mom5/ocean_param/neutral/ocean_nphysics_flux.F90 $(SRCROOT)mom5/ocean_bgc/ocean_ibgc.F90 $(SRCROOT)mom5/ocean_param/gotm-4.0/turbulence/variances.F90 $(SRCROOT)atmos_shared/tracer_driver/tropchem/mo_photo.F90 $(SRCROOT)atmos_param/sea_esf_rad/microphys_cloud.F90 $(SRCROOT)land_lad2/shared/land_tile_diag_sel.F90 $(SRCROOT)mom5/ocean_param/lateral/ocean_bihcgrid_friction.F90 /home/wenchang/tigress/CM2.5_tiger2_wy_killtc/exp/../src/shared/mpp/include/mpp_do_check.h $(SRCROOT)shared/mosaic/read_mosaic.c $(SRCROOT)mom5/ocean_param/sources/ocean_overexchange.F90 $(SRCROOT)shared/mosaic/interp.c /home/wenchang/tigress/CM2.5_tiger2_wy_killtc/exp/../src/shared/mpp/include/mpp_error_a_a.h $(SRCROOT)atmos_param/fsrad/hconst.F90 $(SRCROOT)mom5/ocean_param/gotm-4.0/turbulence/genericeq.F90 $(SRCROOT)atmos_param/cosp/cosp_misr_simulator.F90 $(SRCROOT)mom5/ocean_param/gotm-4.0/turbulence/r_ratio.F90 $(SRCROOT)atmos_param/cosp/cosp_stats.F90 $(SRCROOT)coupler/surface_flux.F90 $(SRCROOT)mom5/ocean_param/vertical/ocean_vert_kpp_test.F90 $(SRCROOT)mom5/ocean_core/ocean_sbc.F90 $(SRCROOT)atmos_param/shallow_cu/conv_utilities_k.F90 $(SRCROOT)ocean_shared/generic_tracers/generic_miniBLING.F90 /home/wenchang/tigress/CM2.5_tiger2_wy_killtc/exp/../src/shared/mpp/include/mpp_get_boundary.h /home/wenchang/tigress/CM2.5_tiger2_wy_killtc/exp/../src/shared/mpp/include/mpp_do_update_nest.h $(SRCROOT)atmos_cubed_sphere/model/fv_control.F90 $(SRCROOT)shared/drifters/drifters.F90 $(SRCROOT)atmos_param/clubb/CLUBB_core/constants_clubb.F90 $(SRCROOT)atmos_param/fsrad/mcm_sw_driver.F90 $(SRCROOT)shared/mosaic/grid.F90 $(SRCROOT)land_lad2/shared/land_tile_diag_buff.F90 $(SRCROOT)shared/coupler/coupler_types.F90 $(SRCROOT)atmos_param/moist_processes/moist_processes_utils.F90 $(SRCROOT)land_lad2/snow/snow_tile.F90 $(SRCROOT)atmos_shared/tracer_driver/tropchem/mo_chemini.F90 $(SRCROOT)atmos_cubed_sphere/tools/fv_mp_mod.F90 $(SRCROOT)shared/time_interp/time_interp_external.F90 /home/wenchang/tigress/CM2.5_tiger2_wy_killtc/exp/../src/shared/mpp/include/mpp_domains_reduce.inc /home/wenchang/tigress/CM2.5_tiger2_wy_killtc/exp/../src/shared/mpp/include/mpp_error_s_s.h $(SRCROOT)mom5/ocean_param/sources/ocean_increment_velocity.F90 $(SRCROOT)land_lad2/vegetation/vegn_disturbance.F90 $(SRCROOT)atmos_param/strat_cloud/strat_netcdf.F90 /home/wenchang/tigress/CM2.5_tiger2_wy_killtc/exp/../src/shared/mpp/include/mpp_do_get_boundary.h $(SRCROOT)mom5/ocean_bgc/ocean_pert_co2.F90 $(SRCROOT)atmos_param/clubb/CLUBB_core/mono_flux_limiter.F90 $(SRCROOT)mom5/ocean_bgc/ocean_po4_pre.F90 $(SRCROOT)atmos_cubed_sphere/model/fv_grid_utils.F90 /home/wenchang/tigress/CM2.5_tiger2_wy_killtc/exp/../src/shared/mpp/include/mpp_read_compressed.h $(SRCROOT)mom5/ocean_param/lateral/ocean_bih_tracer.F90 $(SRCROOT)atmos_cubed_sphere/model/fv_update_phys.F90 $(SRCROOT)atmos_param/entrain/entrain.F90 $(SRCROOT)shared/horiz_interp/test_horiz_interp.F90 $(SRCROOT)atmos_param/donner_deep/fms_donner.F90 $(SRCROOT)atmos_param/clubb/MG_microp_3D.F90 $(SRCROOT)mom5/ocean_tracers/ocean_residency_meta.F90 /home/wenchang/tigress/CM2.5_tiger2_wy_killtc/exp/../src/shared/mpp/include/mpp_comm_sma.inc $(SRCROOT)atmos_cubed_sphere/driver/coupled/fv_physics.F90 $(SRCROOT)mom5/ocean_param/vertical/ocean_vert_mix.F90 /home/wenchang/tigress/CM2.5_tiger2_wy_killtc/exp/../src/shared/mpp/include/mpp_util_sma.inc $(SRCROOT)atmos_param/clubb/CLUBB_core/clubb_precision.F90 $(SRCROOT)mom5/ocean_wave/ocean_wave.F90 $(SRCROOT)atmos_param/sea_esf_rad/gas_tf.F90 $(SRCROOT)land_lad2/shared/table_printer.F90 $(SRCROOT)atmos_cubed_sphere/tools/fv_io.F90 $(SRCROOT)atmos_param/cosp/llnl/prec_scops.F $(SRCROOT)shared/platform/platform.F90 $(SRCROOT)atmos_shared/tracer_driver/tropchem/moz.mods.F90 $(SRCROOT)land_lad2/shared/land_utils.F90 $(SRCROOT)atmos_param/clubb/CLUBB_core/endian.F90 $(SRCROOT)mom5/ocean_bgc/ocmip2_co2calc.F90 $(SRCROOT)atmos_param/vert_diff/vert_diff.F90 /home/wenchang/tigress/CM2.5_tiger2_wy_killtc/exp/../src/shared/mpp/include/mpp_data_mpi.inc /home/wenchang/tigress/CM2.5_tiger2_wy_killtc/exp/../src/shared/mpp/include/mpp_transmit_nocomm.h $(SRCROOT)atmos_param/clubb/CLUBB_core/sponge_layer_damping.F90 $(SRCROOT)shared/fft/fft.F90 $(SRCROOT)atmos_param/clubb/CLUBB_core/pdf_closure_module.F90 $(SRCROOT)atmos_param/mg_drag/mg_drag.F90 $(SRCROOT)mom5/ocean_param/gotm-4.0/turbulence/kbeq.F90 $(SRCROOT)atmos_shared/tracer_driver/atmos_ch3i.F90 /home/wenchang/tigress/CM2.5_tiger2_wy_killtc/exp/../src/shared/field_manager/parse.inc $(SRCROOT)atmos_param/vert_turb_driver/vert_turb_driver.F90 $(SRCROOT)atmos_param/donner_deep/nonfms_donner.F90 /home/wenchang/tigress/CM2.5_tiger2_wy_killtc/exp/../src/shared/mpp/include/mpp_io_misc.inc /home/wenchang/tigress/CM2.5_tiger2_wy_killtc/exp/../src/shared/mpp/include/mpp_global_field.h $(SRCROOT)atmos_param/strat_cloud/gamma_mg.F90 $(SRCROOT)atmos_cubed_sphere/model/fv_arrays.F90 $(SRCROOT)atmos_param/clubb/CLUBB_core/stat_file_module.F90 $(SRCROOT)shared/oda_tools/oda_types.F90 $(SRCROOT)ocean_shared/generic_tracers/generic_COBALT.F90 $(SRCROOT)atmos_param/cosp/cosp_lidar.F90 $(SRCROOT)atmos_shared/tracer_driver/atmos_radon.F90 $(SRCROOT)atmos_param/clubb/CLUBB_core/gmres_wrap.F90 $(SRCROOT)land_lad2/vegetation/vegn_photosynthesis.F90 $(SRCROOT)atmos_param/clubb/CLUBB_core/output_netcdf.F90 $(SRCROOT)atmos_param/sea_esf_rad/aerosolrad_package.F90 /home/wenchang/tigress/CM2.5_tiger2_wy_killtc/exp/../src/shared/include/fms_platform.h $(SRCROOT)atmos_param/sea_esf_rad/rh_based_clouds.F90 $(SRCROOT)mom5/ocean_core/ocean_obc.F90 $(SRCROOT)atmos_shared/tracer_driver/tropchem/mo_chem_utls.F90 $(SRCROOT)shared/horiz_interp/horiz_interp_type.F90 $(SRCROOT)atmos_param/stable_bl_turb/stable_bl_turb.F90 /home/wenchang/tigress/CM2.5_tiger2_wy_killtc/exp/../src/shared/mpp/include/mpp_sum_nocomm.h $(SRCROOT)land_lad2/river/river.F90 $(SRCROOT)mom5/ocean_diag/ocean_velocity_diag.F90 $(SRCROOT)mom5/ocean_tracers/ocean_residency_integrand.F90 $(SRCROOT)mom5/ocean_param/lateral/ocean_mixdownslope.F90 $(SRCROOT)ice_sis/ice_model.F90 $(SRCROOT)atmos_param/clubb/CLUBB_core/parameter_indices.F90 $(SRCROOT)mom5/ocean_tracers/ocean_tracer.F90 $(SRCROOT)shared/mpp/mpp_utilities.F90 /home/wenchang/tigress/CM2.5_tiger2_wy_killtc/exp/../src/shared/mpp/include/mpp_do_updateV_nonblock.h $(SRCROOT)land_lad2/land_data.F90 $(SRCROOT)mom5/ocean_param/lateral/ocean_sigma_transport.F90 $(SRCROOT)atmos_param/cosp/cosp_diagnostics.F90 $(SRCROOT)shared/data_override/data_override.F90 $(SRCROOT)atmos_param/ras/ras.F90 /home/wenchang/tigress/CM2.5_tiger2_wy_killtc/exp/../src/shared/mpp/include/mpp_global_sum.h $(SRCROOT)shared/time_manager/time_manager.F90 $(SRCROOT)land_lad2/soil/soil.F90 $(SRCROOT)atmos_param/sea_esf_rad/isccp_clouds.F90 $(SRCROOT)shared/drifters/quicksort.F90 $(SRCROOT)atmos_param/sea_esf_rad/microphys_rad.F90 $(SRCROOT)atmos_param/clubb/CLUBB_core/stats_LH_zt.F90 $(SRCROOT)mom5/ocean_bgc/ocmip2_he.F90 $(SRCROOT)mom5/ocean_param/gotm-4.0/util/diff_center.F90 $(SRCROOT)mom5/ocean_param/sources/ocean_momentum_source.F90 $(SRCROOT)atmos_param/strat_cloud/nc_cond.F90 /home/wenchang/tigress/CM2.5_tiger2_wy_killtc/exp/../src/mom5/ocean_param/gotm-4.0/include/cppdefs.h $(SRCROOT)atmos_cubed_sphere/model/dyn_core.F90 $(SRCROOT)atmos_shared/tracer_driver/atmos_co2.F90 /home/wenchang/tigress/CM2.5_tiger2_wy_killtc/exp/../src/shared/mosaic/read_mosaic.h $(SRCROOT)mom5/ocean_param/sources/ocean_sponges_eta.F90 $(SRCROOT)mom5/ocean_core/ocean_obc_barotrop.F90 $(SRCROOT)atmos_param/rh_clouds/rh_clouds.F90 $(SRCROOT)shared/diag_manager/diag_util.F90 $(SRCROOT)atmos_param/cosp/quickbeam/zeff.F90 $(SRCROOT)shared/column_diagnostics/column_diagnostics.F90 $(SRCROOT)atmos_param/cosp/MISR_simulator/MISR_simulator.F $(SRCROOT)atmos_param/fsrad/mcm_swnew.F90 $(SRCROOT)atmos_cubed_sphere/tools/fv_climate_nudge.F90 $(SRCROOT)atmos_shared/tracer_driver/atmos_carbon_aerosol.F90 /home/wenchang/tigress/CM2.5_tiger2_wy_killtc/exp/../src/shared/mpp/include/mpp_util_nocomm.inc $(SRCROOT)land_lad2/shared/land_tile_io.F90 $(SRCROOT)atmos_param/clubb/CLUBB_core/T_in_K_module.F90 $(SRCROOT)shared/field_manager/field_manager.F90 $(SRCROOT)atmos_param/sea_esf_rad/longwave_params.F90 $(SRCROOT)atmos_param/clubb/CLUBB_core/input_names.F90 $(SRCROOT)atmos_cubed_sphere/tools/sim_nc_mod.F90 $(SRCROOT)atmos_param/sea_esf_rad/sea_esf_rad.F90 $(SRCROOT)land_lad2/land_constants.F90 $(SRCROOT)mom5/ocean_param/gotm-4.0/turbulence/cmue_d.F90 $(SRCROOT)shared/mpp/mpp.F90 $(SRCROOT)shared/tridiagonal/tridiagonal.F90 $(SRCROOT)atmos_param/sea_esf_rad/optical_path.F90 $(SRCROOT)mom5/ocean_tracers/ocean_tempsalt.F90 $(SRCROOT)atmos_param/clubb/CLUBB_core/advance_windm_edsclrm_module.F90 $(SRCROOT)atmos_param/clubb/CLUBB_core/file_functions.F90 $(SRCROOT)mom5/ocean_param/lateral/ocean_lap_tracer.F90 $(SRCROOT)mom5/ocean_param/gotm-4.0/turbulence/cmue_rf.F90 $(SRCROOT)mom5/ocean_tracers/ocean_residency.F90 $(SRCROOT)atmos_param/clubb/CLUBB_core/variables_diagnostic_module.F90 /home/wenchang/tigress/CM2.5_tiger2_wy_killtc/exp/../src/shared/mpp/include/mpp_io_write.inc /home/wenchang/tigress/CM2.5_tiger2_wy_killtc/exp/../src/shared/fms/read_data_4d.inc $(SRCROOT)mom5/ocean_tracers/ocean_residency_range.F90 $(SRCROOT)mom5/ocean_blobs/ocean_blob_dynamic_free.F90 $(SRCROOT)atmos_cubed_sphere/tools/read_climate_nudge_data.F90 $(SRCROOT)atmos_param/cloud_generator/cloud_generator.F90 $(SRCROOT)shared/sat_vapor_pres/sat_vapor_pres_k.F90 $(SRCROOT)shared/mosaic/gradient_c2l.c $(SRCROOT)atmos_param/sea_esf_rad/zetac_clouds_W.F90 $(SRCROOT)ocean_shared/generic_tracers/generic_CFC.F90 $(SRCROOT)mom5/ocean_param/sources/ocean_xlandinsert.F90 $(SRCROOT)mom5/ocean_param/gotm-4.0/turbulence/cmue_b.F90 $(SRCROOT)atmos_param/clubb/CLUBB_core/stats_subs.F90 $(SRCROOT)mom5/ocean_param/sources/ocean_overflow_OFP.F90 $(SRCROOT)shared/mpp/mpp_parameter.F90 $(SRCROOT)atmos_param/clubb/CLUBB_core/stats_variables.F90 $(SRCROOT)atmos_param/lin_cloud_microphys/lin_cloud_microphys.F90 $(SRCROOT)atmos_param/sea_esf_rad/cloudrad_diagnostics.F90 /home/wenchang/tigress/CM2.5_tiger2_wy_killtc/exp/../src/shared/mpp/include/system_clock.h $(SRCROOT)atmos_param/cosp/icarus-scops-3.7/isccp_cloud_types.F $(SRCROOT)shared/random_numbers/MersenneTwister.F90 $(SRCROOT)land_lad2/shared/nf_utils/nfc.F90 $(SRCROOT)atmos_shared/tracer_driver/tropchem/mo_usrrxt.F90 $(SRCROOT)atmos_shared/vert_advection/vert_advection.F90 $(SRCROOT)shared/topography/topography.F90 $(SRCROOT)shared/mpp/test_mpp.F90 /home/wenchang/tigress/CM2.5_tiger2_wy_killtc/exp/../src/shared/mpp/include/mpp_do_updateV.h /home/wenchang/tigress/CM2.5_tiger2_wy_killtc/exp/../src/shared/mpp/include/mpp_error_a_s.h $(SRCROOT)mom5/ocean_core/ocean_coriolis.F90 $(SRCROOT)shared/mpp/affinity.c $(SRCROOT)mom5/ocean_tracers/transport_matrix.F90 /home/wenchang/tigress/CM2.5_tiger2_wy_killtc/exp/../src/shared/mpp/include/mpp_write_compressed.h $(SRCROOT)atmos_cubed_sphere/model/fv_mapz.F90 $(SRCROOT)atmos_shared/tracer_driver/tropchem/mo_chemdr.F90 $(SRCROOT)atmos_param/shallow_cu/conv_closures.F90 /home/wenchang/tigress/CM2.5_tiger2_wy_killtc/exp/../src/atmos_param/donner_deep/donner_types.h $(SRCROOT)atmos_param/clubb/CLUBB_core/advance_wp2_wp3_module.F90 $(SRCROOT)atmos_param/clubb/CLUBB_core/output_grads.F90 $(SRCROOT)atmos_shared/tracer_driver/tropchem/mo_fastjx.F90 $(SRCROOT)atmos_param/cosp/quickbeam/load_hydrometeor_classes.F90 $(SRCROOT)mom5/ocean_core/ocean_workspace.F90 $(SRCROOT)atmos_param/cosp/cosp_io.F90 $(SRCROOT)atmos_param/cloud_zonal/cloud_zonal.F90 $(SRCROOT)atmos_shared/tracer_driver/atmos_sulfur_hex.F90 $(SRCROOT)mom5/ocean_param/sources/ocean_shortwave_csiro.F90 $(SRCROOT)mom5/ocean_bgc/ocmip2_cfc.F90 $(SRCROOT)land_lad2/glacier/glacier.F90 /home/wenchang/tigress/CM2.5_tiger2_wy_killtc/exp/../src/shared/mpp/include/mpp_sum_mpi.h $(SRCROOT)atmos_shared/tracer_driver/tropchem/mo_rodas_slv.F90 $(SRCROOT)land_lad2/land_model.F90 $(SRCROOT)shared/oda_tools/oda_core.F90 $(SRCROOT)atmos_param/shallow_cu/conv_utilities.F90 $(SRCROOT)shared/time_manager/get_cal_time.F90 $(SRCROOT)atmos_param/clubb/CLUBB_core/hyper_diffusion_4th_ord.F90 $(SRCROOT)atmos_cubed_sphere/model/fv_nesting.F90 $(SRCROOT)mom5/ocean_param/gotm-4.0/turbulence/compute_cpsi3.F90 $(SRCROOT)mom5/ocean_param/vertical/ocean_form_drag.F90 $(SRCROOT)atmos_param/cu_mo_trans/cu_mo_trans.F90 $(SRCROOT)atmos_param/sea_esf_rad/cloud_spec.F90 $(SRCROOT)land_lad2/lake/lake_tile.F90 $(SRCROOT)shared/diag_manager/diag_output.F90 $(SRCROOT)atmos_param/clubb/CLUBB_core/clip_explicit.F90 $(SRCROOT)mom5/ocean_param/vertical/ocean_vert_tidal_test.F90 $(SRCROOT)mom5/ocean_param/sources/ocean_shortwave_jerlov.F90 $(SRCROOT)atmos_param/clubb/CLUBB_core/stats_sfc.F90 $(SRCROOT)atmos_param/cosp/MODIS_simulator/modis_simulator.F90 $(SRCROOT)atmos_param/clubb/CLUBB_core/saturation.F90 $(SRCROOT)atmos_cubed_sphere/tools/external_ic.F90 /home/wenchang/tigress/CM2.5_tiger2_wy_killtc/exp/../src/atmos_param/donner_deep/donner_nml.h /home/wenchang/tigress/CM2.5_tiger2_wy_killtc/exp/../src/land_lad2/vegetation/read_remap_cohort_data_new.inc $(SRCROOT)atmos_param/donner_deep/donner_cloud_model_k.F90 $(SRCROOT)atmos_param/cosp/llnl/llnl_stats.F90 $(SRCROOT)atmos_param/cosp/cosp_driver.F90 $(SRCROOT)land_lad2/shared/land_debug.F90 $(SRCROOT)land_lad2/land_tile.F90 $(SRCROOT)mom5/ocean_param/gotm-4.0/turbulence/cmue_a.F90 $(SRCROOT)atmos_cubed_sphere/model/fv_sg.F90 $(SRCROOT)atmos_param/cosp/icarus-scops-3.7/icarus.F90 $(SRCROOT)atmos_shared/tracer_driver/atmos_sulfate.F90 $(SRCROOT)atmos_param/damping_driver/damping_driver.F90 $(SRCROOT)atmos_param/betts_miller/bm_omp.F90 /home/wenchang/tigress/CM2.5_tiger2_wy_killtc/exp/../src/land_lad2/shared/nf_utils/getput.inc $(SRCROOT)mom5/ocean_param/gotm-4.0/util/eqstate.F90 $(SRCROOT)atmos_param/cosp/actsim/lmd_ipsl_stats.F90 $(SRCROOT)atmos_shared/tracer_driver/tropchem/moz.mat.F90 $(SRCROOT)atmos_param/clubb/CLUBB_core/clip_semi_implicit.F90 $(SRCROOT)atmos_param/clubb/CLUBB_core/parameters_microphys.F90 $(SRCROOT)atmos_param/clouds/clouds.F90 /home/wenchang/tigress/CM2.5_tiger2_wy_killtc/exp/../src/shared/mpp/include/mpp_domains_util.inc $(SRCROOT)mom5/ocean_param/sources/ocean_rivermix.F90 /home/wenchang/tigress/CM2.5_tiger2_wy_killtc/exp/../src/shared/mpp/include/mpp_update_domains2D.h $(SRCROOT)atmos_param/clubb/CLUBB_core/mean_adv.F90 $(SRCROOT)shared/mpp/test_mpp_pset.F90 $(SRCROOT)shared/fms/fms_io.F90 $(SRCROOT)atmos_shared/tracer_driver/stratchem/strat_chem_driver.F90 $(SRCROOT)mom5/ocean_param/vertical/ocean_convect.F90 $(SRCROOT)atmos_param/sea_esf_rad/donner_deep_clouds_W.F90 /home/wenchang/tigress/CM2.5_tiger2_wy_killtc/exp/../src/atmos_param/cosp/icarus-scops-3.7/congvec.H $(SRCROOT)atmos_param/cg_drag/cg_drag.F90 $(SRCROOT)atmos_param/donner_deep/wet_deposition_0D.F90 $(SRCROOT)ice_sis/ice_bergs.F90 $(SRCROOT)atmos_shared/tracer_driver/tropchem/mo_jpl.F90 $(SRCROOT)atmos_param/fsrad/rdparm.F90 $(SRCROOT)mom5/ocean_core/ocean_barotropic.F90 $(SRCROOT)ocean_shared/generic_tracers/generic_tracer_utils.F90 /home/wenchang/tigress/CM2.5_tiger2_wy_killtc/exp/../src/shared/mosaic/interp.h $(SRCROOT)mom5/ocean_bgc/ocean_bgc_restore.F90 $(SRCROOT)atmos_shared/tracer_driver/aer_ccn_act/aerosol_params.F90 $(SRCROOT)mom5/ocean_param/lateral/ocean_bih_friction.F90 $(SRCROOT)atmos_param/cosp/llnl/pf_to_mr.F $(SRCROOT)land_lad2/vegetation/vegn_cohort_io.F90 $(SRCROOT)land_lad2/soil/soil_tile.F90 $(SRCROOT)mom5/ocean_wave/ocean_wave_type.F90 $(SRCROOT)atmos_param/cosp/quickbeam/load_mie_table.F90 $(SRCROOT)ice_sis/mask.F90 $(SRCROOT)atmos_param/clubb/CLUBB_3D_var.F90 /home/wenchang/tigress/CM2.5_tiger2_wy_killtc/exp/../src/shared/mpp/include/mpp_data_sma.inc /home/wenchang/tigress/CM2.5_tiger2_wy_killtc/exp/../src/shared/mosaic/create_xgrid.h $(SRCROOT)atmos_param/sea_esf_rad/lw_gases_stdtf.F90 $(SRCROOT)mom5/ocean_param/neutral/ocean_nphysicsC.F90 $(SRCROOT)mom5/ocean_param/sources/ocean_overflow.F90 $(SRCROOT)mom5/ocean_param/gotm-4.0/turbulence/compute_rist.F90 $(SRCROOT)atmos_cubed_sphere/tools/test_cases.F90 $(SRCROOT)atmos_param/strat_cloud/aerosol_cloud.F90 $(SRCROOT)atmos_param/clubb/CLUBB_core/parameters_tunable.F90 $(SRCROOT)mom5/ocean_core/ocean_bbc.F90 $(SRCROOT)atmos_param/lscale_cond/lscale_cond.F90 $(SRCROOT)atmos_param/fsrad/shortwave.F90 /home/wenchang/tigress/CM2.5_tiger2_wy_killtc/exp/../src/shared/mpp/include/mpp_do_redistribute.h $(SRCROOT)atmos_cubed_sphere/model/fv_tracer2d.F90 $(SRCROOT)atmos_shared/tracer_driver/tropchem/mo_read_sim_chm.F90 $(SRCROOT)mom5/ocean_param/lateral/ocean_lap_friction.F90 $(SRCROOT)mom5/ocean_core/ocean_density.F90 $(SRCROOT)land_lad2/vegetation/vegn_tile.F90 $(SRCROOT)atmos_param/sea_esf_rad/original_fms_rad.F90 $(SRCROOT)mom5/ocean_param/gotm-4.0/turbulence/cmue_sg.F90 $(SRCROOT)mom5/ocean_bgc/ocean_generic_tracer.F90 $(SRCROOT)atmos_param/donner_deep/cumulus_closure_k.F90 $(SRCROOT)mom5/ocean_param/sources/ocean_sponges_velocity.F90 $(SRCROOT)mom5/ocean_tracers/ocean_passive.F90 $(SRCROOT)atmos_param/sea_esf_rad/longwave_fluxes.F90 $(SRCROOT)land_lad2/shared/sphum.F90 $(SRCROOT)mom5/ocean_param/vertical/ocean_vert_const.F90 $(SRCROOT)shared/axis_utils/axis_utils.F90 /home/wenchang/tigress/CM2.5_tiger2_wy_killtc/exp/../src/shared/mpp/include/mpp_sum.inc $(SRCROOT)atmos_param/shallow_cu/conv_plumes_k.F90 $(SRCROOT)ice_sis/ice_spec.F90 $(SRCROOT)atmos_param/donner_deep/donner_rad_k.F90 /home/wenchang/tigress/CM2.5_tiger2_wy_killtc/exp/../src/shared/mpp/include/mpp_util_mpi.inc /home/wenchang/tigress/CM2.5_tiger2_wy_killtc/exp/../src/shared/drifters/drifters_push.h $(SRCROOT)atmos_shared/tracer_driver/tropchem/mo_setsox.F90 $(SRCROOT)mom5/ocean_param/gotm-4.0/util/lagrange.F90 $(SRCROOT)atmos_shared/tracer_driver/aer_ccn_act/aer_ccn_act_k.F90 /home/wenchang/tigress/CM2.5_tiger2_wy_killtc/exp/../src/shared/mpp/include/mpp_do_checkV.h $(SRCROOT)mom5/ocean_param/gotm-4.0/turbulence/cmue_c.F90 $(SRCROOT)atmos_param/clubb/CLUBB_core/diffusion.F90 $(SRCROOT)atmos_param/radiation_driver/radiation_driver.F90 $(SRCROOT)atmos_param/clubb/CLUBB_core/pdf_parameter_module.F90 $(SRCROOT)shared/horiz_interp/horiz_interp_conserve.F90 /home/wenchang/tigress/CM2.5_tiger2_wy_killtc/exp/../src/shared/mpp/include/mpp_sum_sma.h $(SRCROOT)atmos_param/strat_cloud/strat_cloud_legacy.F90 $(SRCROOT)atmos_cubed_sphere/tools/init_hydro.F90 $(SRCROOT)atmos_shared/tracer_driver/tropchem/tropchem_driver.F90 $(SRCROOT)atmos_param/sea_esf_rad/longwave_tables.F90 $(SRCROOT)atmos_cubed_sphere/driver/coupled/atmosphere.F90 /home/wenchang/tigress/CM2.5_tiger2_wy_killtc/exp/../src/shared/mpp/include/mpp_write.h $(SRCROOT)atmos_param/grey_radiation/grey_radiation.F90 $(SRCROOT)shared/diag_manager/diag_manager.F90 $(SRCROOT)mom5/ocean_core/ocean_parameters.F90 $(SRCROOT)mom5/ocean_diag/ocean_adv_vel_diag.F90 $(SRCROOT)mom5/ocean_core/ocean_domains.F90 /home/wenchang/tigress/CM2.5_tiger2_wy_killtc/exp/../src/shared/mpp/include/mpp_transmit_sma.h $(SRCROOT)mom5/ocean_core/ocean_util.F90 $(SRCROOT)shared/fms/fms.F90 $(SRCROOT)land_lad2/snow/snow.F90 $(SRCROOT)mom5/ocean_param/sources/ocean_sponges_tracer.F90 $(SRCROOT)atmos_param/clubb/alt_cloud.F90 $(SRCROOT)atmos_cubed_sphere/tools/fv_restart.F90 $(SRCROOT)shared/astronomy/astronomy.F90 $(SRCROOT)land_lad2/river/river_physics.F90 $(SRCROOT)mom5/ocean_param/neutral/ocean_nphysics_util.F90 $(SRCROOT)shared/oda_tools/xbt_drop_rate_adjust.f90 $(SRCROOT)mom5/ocean_param/vertical/ocean_vert_pp.F90 $(SRCROOT)mom5/ocean_param/gotm-4.0/turbulence/ispralength.F90 $(SRCROOT)atmos_shared/tracer_driver/tropchem/mo_exp_slv.F90 /home/wenchang/tigress/CM2.5_tiger2_wy_killtc/exp/../src/shared/mpp/include/mpp_global_reduce.h $(SRCROOT)atmos_param/clubb/CLUBB_core/error_code.F90 $(SRCROOT)atmos_shared/tracer_driver/aer_ccn_act/aer_in_act.F90 $(SRCROOT)atmos_param/fsrad/fsrad.F90 $(SRCROOT)atmos_param/cosp/icarus-scops-3.7/scops.F $(SRCROOT)atmos_param/fsrad/rad_diag.F90 $(SRCROOT)shared/memutils/memuse.c $(SRCROOT)atmos_param/cosp/llnl/cosp_radar.F90 /home/wenchang/tigress/CM2.5_tiger2_wy_killtc/exp/../src/land_lad2/canopy_air/../shared/debug.inc $(SRCROOT)mom5/ocean_param/lateral/ocean_lapcst_friction.F90 /home/wenchang/tigress/CM2.5_tiger2_wy_killtc/exp/../src/shared/mpp/include/mpp_comm_nocomm.inc $(SRCROOT)mom5/ocean_tracers/ocean_frazil.F90 $(SRCROOT)atmos_param/fsrad/mcm_swtbls.F90 /home/wenchang/tigress/CM2.5_tiger2_wy_killtc/exp/../src/shared/mpp/include/mpp_chksum.h $(SRCROOT)atmos_param/strat_cloud/strat_cloud_utilities.F90 $(SRCROOT)atmos_param/clubb/CLUBB_core/grid_class.F90 $(SRCROOT)atmos_param/donner_deep/donner_deep_k.F90 $(SRCROOT)ocean_shared/generic_tracers/generic_tracer.F90 $(SRCROOT)atmos_param/strat_cloud/microphysics.F90 $(SRCROOT)atmos_shared/tracer_driver/stratchem/strat_chem_model.F90 $(SRCROOT)atmos_shared/tracer_driver/tropchem/mo_hook.F90 $(SRCROOT)atmos_param/tke_turb/tke_turb.F90 $(SRCROOT)atmos_param/moist_processes/moistproc_kernels.F90 $(SRCROOT)shared/diag_manager/diag_data.F90 $(SRCROOT)mom5/ocean_param/gotm-4.0/turbulence/kpp.F90 $(SRCROOT)shared/mpp/mpp_memutils.F90 $(SRCROOT)mom5/ocean_param/neutral/ocean_nphysicsB.F90 $(SRCROOT)atmos_param/sea_esf_rad/diag_clouds_W.F90 $(SRCROOT)atmos_param/donner_deep/donner_lscloud_k.F90 $(SRCROOT)shared/fms/test_fms_io.F90 $(SRCROOT)shared/horiz_interp/horiz_interp_spherical.F90 $(SRCROOT)atmos_param/strat_cloud/simple_pdf.F90 $(SRCROOT)land_lad2/glacier/glac_tile.F90 $(SRCROOT)mom5/ocean_param/vertical/ocean_vert_kpp_mom4p1.F90 $(SRCROOT)atmos_shared/tracer_driver/atmos_soa.F90 $(SRCROOT)mom5/ocean_tracers/ocean_tpm.F90 $(SRCROOT)atmos_param/cosp/cosp_isccp_simulator.F90 $(SRCROOT)shared/tracer_manager/tracer_manager.F90 $(SRCROOT)atmos_param/cosp/quickbeam/array_lib.F90 $(SRCROOT)atmos_param/clubb/CLUBB_core/gmres_cache.F90 $(SRCROOT)atmos_param/cosp/quickbeam/radar_simulator.F90 $(SRCROOT)atmos_param/cosp/cosp.F90 $(SRCROOT)land_lad2/soil/uptake.F90 $(SRCROOT)mom5/ocean_param/sources/ocean_riverspread.F90 $(SRCROOT)shared/horiz_interp/horiz_interp_bilinear.F90 $(SRCROOT)land_lad2/vegetation/vegn_data.F90 $(SRCROOT)land_lad2/topo_rough/topo_rough.F90 $(SRCROOT)shared/mosaic/mosaic_util.c $(SRCROOT)shared/exchange/xgrid.F90 $(SRCROOT)atmos_param/donner_deep/donner_cape_k.F90 $(SRCROOT)atmos_param/clubb/advance_sclrm_Nd_module.F90 $(SRCROOT)atmos_cubed_sphere/tools/fv_nudge.F90 /home/wenchang/tigress/CM2.5_tiger2_wy_killtc/exp/../src/shared/mpp/include/mpp_update_domains2D_nonblock.h $(SRCROOT)mom5/ocean_bgc/ocmip2_biotic.F90 $(SRCROOT)land_lad2/vegetation/vegn_dynamics.F90 /home/wenchang/tigress/CM2.5_tiger2_wy_killtc/exp/../src/shared/mosaic/constant.h $(SRCROOT)land_lad2/shared/land_io.F90 $(SRCROOT)atmos_param/cloud_rad/cloud_rad.F90 $(SRCROOT)shared/drifters/drifters_input.F90 $(SRCROOT)mom5/ocean_core/ocean_types.F90 $(SRCROOT)atmos_param/clubb/CLUBB_core/stats_zm.F90 $(SRCROOT)ice_param/ocean_rough.F90 /home/wenchang/tigress/CM2.5_tiger2_wy_killtc/exp/../src/atmos_param/strat_cloud/strat_nml.h $(SRCROOT)ice_sis/ice_thm.F90 /home/wenchang/tigress/CM2.5_tiger2_wy_killtc/exp/../src/shared/mpp/include/mpp_reduce_nocomm.h /home/wenchang/tigress/CM2.5_tiger2_wy_killtc/exp/../src/shared/mpp/include/mpp_io_util.inc $(SRCROOT)coupler/coupler_main.F90 $(SRCROOT)atmos_param/sea_esf_rad/cloudrad_package.F90 $(SRCROOT)atmos_param/topo_drag/topo_drag.F90 /home/wenchang/tigress/CM2.5_tiger2_wy_killtc/exp/../src/shared/mpp/include/mpp_error_s_a.h $(SRCROOT)atmos_shared/tracer_driver/atmos_convection_tracer.F90 $(SRCROOT)shared/oda_tools/write_ocean_data.F90 $(SRCROOT)atmos_param/fsrad/co2int.F90 $(SRCROOT)atmos_param/sea_esf_rad/uw_clouds_W.F90 $(SRCROOT)atmos_param/strat_cloud/strat_cloud.F90 $(SRCROOT)atmos_param/cosp/cosp_simulator.F90 $(SRCROOT)atmos_cubed_sphere/model/boundary.F90 $(SRCROOT)shared/time_interp/time_interp.F90 $(SRCROOT)atmos_param/clubb/CLUBB_core/mixing_length.F90 $(SRCROOT)shared/horiz_interp/horiz_interp.F90 $(SRCROOT)atmos_param/cosp/quickbeam/mrgrnk.F90 $(SRCROOT)atmos_param/sea_esf_rad/specified_clouds_W.F90 $(SRCROOT)land_lad2/river/river_type.F90 /home/wenchang/tigress/CM2.5_tiger2_wy_killtc/exp/../src/shared/mpp/include/mpp_util.inc $(SRCROOT)atmos_param/clubb/CLUBB_core/stats_rad_zt.F90 $(SRCROOT)atmos_shared/atmos_nudge/atmos_nudge.F90 $(SRCROOT)atmos_shared/tracer_driver/aer_ccn_act/ice_nucl.F90 $(SRCROOT)mom5/ocean_param/neutral/ocean_nphysics_diff.F90 $(SRCROOT)mom5/ocean_blobs/ocean_blob_static_free.F90 $(SRCROOT)mom5/ocean_param/neutral/ocean_nphysics.F90 $(SRCROOT)shared/mosaic/create_xgrid.c /home/wenchang/tigress/CM2.5_tiger2_wy_killtc/exp/../src/shared/mpp/include/mpp_transmit_mpi.h $(SRCROOT)land_lad2/shared/land_numerics.F90 $(SRCROOT)atmos_param/strat_cloud/rotstayn_klein_mp.F90 $(SRCROOT)atmos_param/shallow_cu/deep_conv.F90 $(SRCROOT)atmos_cubed_sphere/model/tp_core.F90 $(SRCROOT)atmos_param/clubb/CLUBB_core/model_flags.F90 $(SRCROOT)mom5/ocean_param/gotm-4.0/turbulence/tkeeq.F90 $(SRCROOT)atmos_shared/interpolator/interpolator.F90 /home/wenchang/tigress/CM2.5_tiger2_wy_killtc/exp/../src/shared/mpp/include/mpp_io_read.inc $(SRCROOT)land_lad2/lake/lake.F90 $(SRCROOT)atmos_param/sea_esf_rad/ozone.F90 $(SRCROOT)atmos_param/cosp/quickbeam/math_lib.F90 $(SRCROOT)mom5/ocean_param/lateral/ocean_bihcst_friction.F90 $(SRCROOT)atmos_shared/tracer_driver/atmos_tracer_utilities.F90 $(SRCROOT)mom5/ocean_param/gotm-4.0/turbulence/production.F90 $(SRCROOT)atmos_shared/tracer_driver/tropchem/mo_imp_slv.F90 $(SRCROOT)atmos_param/donner_deep/donner_types.F90 $(SRCROOT)atmos_param/donner_deep/donner_lite_k.F90 $(SRCROOT)atmos_param/sea_esf_rad/longwave_clouds.F90 $(SRCROOT)atmos_param/sea_esf_rad/radiation_diag.F90 $(SRCROOT)ice_param/ocean_albedo.F90 $(SRCROOT)ocean_shared/generic_tracers/generic_ERGOM.F90 $(SRCROOT)shared/mosaic/mosaic.F90 $(SRCROOT)mom5/ocean_param/gotm-4.0/turbulence/internal_wave.F90 $(SRCROOT)atmos_param/clubb/CLUBB_core/stats_type.F90 $(SRCROOT)atmos_param/clubb/CLUBB_core/lapack_wrap.F90 $(SRCROOT)shared/topography/gaussian_topog.F90 /home/wenchang/tigress/CM2.5_tiger2_wy_killtc/exp/../src/mom5/ocean_param/gotm-4.0/include/version.h $(SRCROOT)atmos_shared/tracer_driver/atmos_tracer_driver.F90 $(SRCROOT)mom5/ocean_param/gotm-4.0/util/diff_face.F90 $(SRCROOT)atmos_param/cosp/quickbeam/atmos_lib.F90 $(SRCROOT)land_lad2/vegetation/vegn_static_override.F90 /home/wenchang/tigress/CM2.5_tiger2_wy_killtc/exp/../src/shared/mpp/include/mpp_define_nest_domains.inc $(SRCROOT)atmos_param/strat_cloud/polysvp.F90 $(SRCROOT)atmos_param/donner_deep/donner_utilities_k.F90 $(SRCROOT)atmos_param/sea_esf_rad/esfsw_driver.F90 /home/wenchang/tigress/CM2.5_tiger2_wy_killtc/exp/../src/shared/mpp/include/mpp_domains_comm.inc $(SRCROOT)mom5/ocean_param/sources/ocean_xlandmix.F90 $(SRCROOT)atmos_param/cloud_obs/cloud_obs.F90 $(SRCROOT)atmos_param/cosp/cosp_rttov_simulator.F90 /home/wenchang/tigress/CM2.5_tiger2_wy_killtc/exp/../src/shared/mpp/include/mpp_gather.h $(SRCROOT)atmos_param/cosp/cosp_utils.F90 /home/wenchang/tigress/CM2.5_tiger2_wy_killtc/exp/../src/shared/mpp/include/mpp_transmit.inc $(SRCROOT)mom5/ocean_tracers/ocean_age_tracer.F90 $(SRCROOT)land_lad2/vegetation/vegn_harvesting.F90 $(SRCROOT)atmos_param/moist_processes/detr_ice_num.F90 $(SRCROOT)mom5/ocean_core/ocean_model.F90 $(SRCROOT)mom5/ocean_param/gotm-4.0/util/util.F90 $(SRCROOT)land_lad2/canopy_air/cana_tile.F90 $(SRCROOT)mom5/ocean_param/lateral/ocean_lapgen_friction.F90 $(SRCROOT)shared/random_numbers/random_numbers.F90 $(SRCROOT)mom5/ocean_tracers/ocean_residency_ml.F90 /home/wenchang/tigress/CM2.5_tiger2_wy_killtc/exp/../src/land_lad2/vegetation/vegn_cohort_io.inc $(SRCROOT)atmos_cubed_sphere/model/sw_core.F90 $(SRCROOT)mom5/ocean_param/vertical/ocean_vert_chen.F90 /home/wenchang/tigress/CM2.5_tiger2_wy_killtc/exp/../src/shared/fms/read_data_3d.inc $(SRCROOT)mom5/ocean_param/lateral/ocean_bihgen_friction.F90 $(SRCROOT)atmos_param/sea_esf_rad/sealw99.F90 /home/wenchang/tigress/CM2.5_tiger2_wy_killtc/exp/../src/shared/mosaic/mosaic_util.h $(SRCROOT)atmos_shared/tracer_driver/tropchem/mo_fphoto.F90 $(SRCROOT)atmos_param/sea_esf_rad/mgrp_prscr_clds.F90 $(SRCROOT)atmos_param/sea_esf_rad/radiative_gases.F90 $(SRCROOT)mom5/ocean_param/gotm-4.0/turbulence/tkealgebraic.F90 $(SRCROOT)mom5/ocean_core/oda_driver.F90 /home/wenchang/tigress/CM2.5_tiger2_wy_killtc/exp/../src/shared/drifters/drifters_set_field.h $(SRCROOT)mom5/ocean_blobs/ocean_blob_dynamic_bottom.F90 /home/wenchang/tigress/CM2.5_tiger2_wy_killtc/exp/../src/shared/mpp/include/mpp_read_2Ddecomp.h $(SRCROOT)mom5/ocean_param/sources/ocean_shortwave.F90 $(SRCROOT)atmos_param/sea_esf_rad/bulkphys_rad.F90 $(SRCROOT)atmos_cubed_sphere/tools/fv_ada_nudge.F90 $(SRCROOT)mom5/ocean_param/vertical/ocean_vert_kpp_mom4p0.F90 $(SRCROOT)atmos_param/edt/edt.F90 $(SRCROOT)atmos_param/sea_esf_rad/strat_clouds_W.F90 /home/wenchang/tigress/CM2.5_tiger2_wy_killtc/exp/../src/atmos_param/cosp/cosp_defs.H $(SRCROOT)mom5/ocean_blobs/ocean_blob_util.F90 $(SRCROOT)atmos_param/vert_diff_driver/vert_diff_driver.F90 /home/wenchang/tigress/CM2.5_tiger2_wy_killtc/exp/../src/shared/mpp/include/mpp_comm.inc $(SRCROOT)ice_sis/rot.F90 $(SRCROOT)atmos_param/cloud_generator/betaDistribution.F90 $(SRCROOT)atmos_param/donner_deep/donner_deep_miz.F90 $(SRCROOT)shared/diag_manager/diag_axis.F90 $(SRCROOT)atmos_param/clubb/CLUBB_core/surface_varnce_module.F90 $(SRCROOT)mom5/ocean_core/ocean_velocity_advect.F90 $(SRCROOT)mom5/ocean_param/lateral/ocean_lapcgrid_friction.F90 $(SRCROOT)atmos_param/strat_cloud/micro_mg.F90 $(SRCROOT)land_lad2/shared/nf_utils/nfu.F90 /home/wenchang/tigress/CM2.5_tiger2_wy_killtc/exp/../src/shared/fms/write_data.inc $(SRCROOT)atmos_param/sea_esf_rad/longwave_driver.F90 $(SRCROOT)atmos_param/clubb/CLUBB_core/csr_matrix_class_3array.F90 $(SRCROOT)mom5/ocean_param/gotm-4.0/turbulence/q2over2eq.F90 $(SRCROOT)atmos_param/clubb/CLUBB_core/numerical_check.F90 $(SRCROOT)mom5/ocean_diag/ocean_tracer_util.F90 $(SRCROOT)mom5/ocean_core/ocean_velocity.F90 $(SRCROOT)atmos_cubed_sphere/model/fv_fill.F90 $(SRCROOT)mom5/ocean_tracers/ocean_tracer_advect.F90 $(SRCROOT)shared/horiz_interp/horiz_interp_bicubic.F90 $(SRCROOT)atmos_param/donner_deep/donner_meso_k.F90 $(SRCROOT)atmos_param/sea_esf_rad/esfsw_parameters.F90 $(SRCROOT)mom5/ocean_core/ocean_operators.F90 $(SRCROOT)atmos_param/sea_esf_rad/aerosol.F90 $(SRCROOT)atmos_param/clubb/CLUBB_core/matrix_operations.F90 $(SRCROOT)shared/sat_vapor_pres/sat_vapor_pres.F90 $(SRCROOT)mom5/ocean_param/neutral/ocean_nphysics_tensor.F90 $(SRCROOT)shared/mpp/mpp_io.F90 $(SRCROOT)atmos_shared/tracer_driver/aer_ccn_act/aer_ccn_act.F90 $(SRCROOT)mom5/ocean_param/vertical/ocean_vert_util.F90 $(SRCROOT)atmos_param/fsrad/fs_profile.F90 $(SRCROOT)shared/drifters/drifters_core.F90 $(SRCROOT)atmos_param/strat_cloud/morrison_gettelman_microp.F90 $(SRCROOT)land_lad2/shared/nf_utils/nf_utils.F90 $(SRCROOT)atmos_param/donner_deep/donner_deep.F90 $(SRCROOT)shared/diag_manager/diag_table.F90 $(SRCROOT)mom5/ocean_blobs/ocean_blob.F90 $(SRCROOT)atmos_param/sea_esf_rad/rad_utilities.F90 /home/wenchang/tigress/CM2.5_tiger2_wy_killtc/exp/../src/shared/mpp/include/mpp_do_update_nonblock.h $(SRCROOT)ocean_shared/generic_tracers/generic_BLING.F90 /home/wenchang/tigress/CM2.5_tiger2_wy_killtc/exp/../src/shared/mpp/include/mpp_data_nocomm.inc /home/wenchang/tigress/CM2.5_tiger2_wy_killtc/exp/../src/shared/mpp/include/mpp_update_nest_domains.h $(SRCROOT)shared/mpp/mpp_domains.F90 $(SRCROOT)mom5/ocean_blobs/ocean_blob_diag.F90 /home/wenchang/tigress/CM2.5_tiger2_wy_killtc/exp/../src/shared/mpp/include/mpp_chksum_scalar.h $(SRCROOT)atmos_param/dry_adj/dry_adj.F90 $(SRCROOT)ocean_shared/generic_tracers/generic_TOPAZ.F90 $(SRCROOT)atmos_param/clubb/CLUBB_core/parameters_model.F90 $(SRCROOT)atmos_param/diag_cloud/diag_cloud.F90 $(SRCROOT)atmos_param/clubb/CLUBB_core/stats_rad_zm.F90 $(SRCROOT)atmos_param/clubb/CLUBB_core/advance_xp2_xpyp_module.F90 /home/wenchang/tigress/CM2.5_tiger2_wy_killtc/exp/../src/shared/fms/read_data_2d.inc /home/wenchang/tigress/CM2.5_tiger2_wy_killtc/exp/../src/shared/mpp/include/mpp_domains_define.inc /home/wenchang/tigress/CM2.5_tiger2_wy_killtc/exp/../src/mom5/ocean_core/ocean_memory.h $(SRCROOT)mom5/ocean_diag/ocean_diagnostics.F90 $(SRCROOT)mom5/ocean_param/gotm-4.0/turbulence/potentialml.F90 $(SRCROOT)atmos_cubed_sphere/tools/fv_diagnostics.F90 $(SRCROOT)shared/memutils/memutils.F90 $(SRCROOT)atmos_shared/tracer_driver/tropchem/mo_setinv.F90 $(SRCROOT)mom5/ocean_param/gotm-4.0/turbulence/fk_craig.F90 /home/wenchang/tigress/CM2.5_tiger2_wy_killtc/exp/../src/shared/mpp/include/mpp_scatter.h $(SRCROOT)mom5/ocean_param/gotm-4.0/util/gridinterpol.F90 /home/wenchang/tigress/CM2.5_tiger2_wy_killtc/exp/../src/shared/mpp/include/mpp_reduce_mpi.h $(SRCROOT)atmos_param/diffusivity/diffusivity.F90 $(SRCROOT)mom5/ocean_diag/ocean_drifters.F90 $(SRCROOT)atmos_param/clubb/CLUBB_core/fill_holes.F90 $(SRCROOT)atmos_param/strat_cloud/check_nan.F90 $(SRCROOT)atmos_cubed_sphere/tools/sorted_index.F90 $(SRCROOT)atmos_param/cosp/quickbeam/gases.F90 $(SRCROOT)shared/mpp/threadloc.c $(SRCROOT)atmos_param/shallow_conv/shallow_conv.F90 /home/wenchang/tigress/CM2.5_tiger2_wy_killtc/exp/../src/shared/mosaic/gradient_c2l.h $(SRCROOT)atmos_cubed_sphere/model/null/nh_core.F90 $(SRCROOT)atmos_param/shallow_cu/uw_conv.F90 $(SRCROOT)atmos_param/clubb/CLUBB_core/corr_matrix_module.F90 $(SRCROOT)atmos_param/clubb/CLUBB_core/calendar.F90 $(SRCROOT)atmos_param/moist_processes/moist_processes.F90 $(SRCROOT)atmos_cubed_sphere/tools/fv_eta.F90 $(SRCROOT)atmos_param/sea_esf_rad/rad_output_file.F90 /home/wenchang/tigress/CM2.5_tiger2_wy_killtc/exp/../src/shared/mpp/include/mpp_write_2Ddecomp.h $(SRCROOT)shared/mpp/nsclock.c $(SRCROOT)atmos_param/clubb/CLUBB_core/pos_definite_module.F90 $(SRCROOT)atmos_param/clubb/CLUBB_core/extrapolation.F90 $(SRCROOT)atmos_param/clubb/CLUBB_core/anl_erf.F90 $(SRCROOT)ice_param/ice_albedo.F90 $(SRCROOT)shared/exchange/stock_constants.F90 $(SRCROOT)atmos_cubed_sphere/tools/external_sst.F90 /home/wenchang/tigress/CM2.5_tiger2_wy_killtc/exp/../src/shared/mpp/include/mpp_io_connect.inc $(SRCROOT)mom5/ocean_param/sources/ocean_increment_tracer.F90 $(SRCROOT)atmos_param/diag_cloud_rad/diag_cloud_rad.F90 $(SRCROOT)mom5/ocean_param/gotm-4.0/turbulence/cmue_ma.F90
clean: neat
	-rm -f .CM2.5.x.cppdefs $(OBJ) CM2.5.x
neat:
	-rm -f $(TMPFILES)
localize: $(OFF)
	cp $(OFF) .
TAGS: $(SRC)
	etags $(SRC)
tags: $(SRC)
	ctags $(SRC)
CM2.5.x: $(OBJ) 
	$(LD) $(OBJ) -o CM2.5.x  $(LDFLAGS)
