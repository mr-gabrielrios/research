## SWISHE small-scale manuscript - figure guide
This `readme` is intended to provide a summary of what each figure is and how to reproduce it.

- `TC_density_TS.png` and `TC_density_C15w.png`
	- Summary: Plot of number of TC days per N-degree x N-degree bin for each model over the first 50 years of each model run.
	- Data: run `tc_analysis.tc_track_data(models=['AM2.5', 'FLOR', 'HIRAM'], experiments=['control', 'swishe'], storm_type='C15w' OR 'TS', year_range=(101, 150))`
	- Figure generation: run `visualization.density_grid(data, model_names=['AM2.5', 'FLOR', 'HIRAM'], bin_resolution=5)`

- `TC_max_wind-TS.png`
	- Summary: distribution of maximum winds for all storms for all models over first 50 years of each model run.
	- Data: run `tc_analysis.tc_track_data(models=['AM2.5', 'FLOR', 'HIRAM'], experiments=['control', 'swish    e'], storm_type='C15w' OR 'TS', year_range=(101, 150))`
	- Figure generation: run `visualization.pdf(track_data, 'max_wnd', num_bins=40)`

- `TC_min_slp-TS-TS.png`
	- Summary: distribution of minimum SLPs for all storms for all models over first 50 years of each model run.
	- Data: run `tc_analysis.tc_track_data(models=['AM2.5', 'FLOR', 'HIRAM'], experiments=['control', 'swish    e'], storm_type='C15w' OR 'TS', year_range=(101, 150))`
	- Figure generation: run `visualization.pdf(track_data, 'min_slp', num_bins=40)`

- `TC_center_lat-TS-genesis.png`
	- Summary: distribution of TC centers for all storms at LMI for all models over first 50 years of each model run.
	- Data: run `tc_analysis.tc_track_data(models=['AM2.5', 'FLOR', 'HIRAM'], experiments=['control', 'swish    e'], storm_type='C15w' OR 'TS', year_range=(101, 150))`
	- Figure generation: run `visualization.pdf(track_data, 'center_lat', num_bins=40)`

- `TC-field_slp-bin_b3.png'
	- Summary: planar composite of sea-level pressure at the maximum intensity bin for all storms.
	- Data: run `tc_analysis.tc_model_data(['AM2.5', 'FLOR', 'HIRAM'], ['control', 'swishe'], num_storms=50)
	- Figure generation: run `visualization.planar_composite(data, ['AM2.5', 'FLOR', 'HIRAM'], ['b3'], 'slp', None,  ['control', 'swishe', 'swishe-control'])`
	- Statistics:
		- Number of records used for the composite for AM2.5, control, b3: 84
		- Number of records used for the composite for AM2.5, swishe, b3: 5
		- Number of records used for the composite for FLOR, control, b3: 322
		- Number of records used for the composite for FLOR, swishe, b3: 20
		- Number of records used for the composite for HIRAM, control, b3: 86
		- Number of records used for the composite for HIRAM, swishe, b3: 24

- `TC-field_WVP-bin_b3.png'
	- Summary: planar composite of column-integ. water vapor at the maximum intensity bin for all storms.
	- Data: run `tc_analysis.tc_model_data(['AM2.5', 'FLOR', 'HIRAM'], ['control', 'swishe'], num_storms=50)
	- Figure generation: run `visualization.planar_composite(data, ['AM2.5', 'FLOR', 'HIRAM'], ['b3'], 'WVP', None,  ['control', 'swishe', 'swishe-control'])`
	- Statistics:
		- Number of records used for the composite for AM2.5, control, b3: 84
		- Number of records used for the composite for AM2.5, swishe, b3: 5
		- Number of records used for the composite for FLOR, control, b3: 322
		- Number of records used for the composite for FLOR, swishe, b3: 20
		- Number of records used for the composite for HIRAM, control, b3: 86
		- Number of records used for the composite for HIRAM, swishe, b3: 24

- `TC-field_shflx-bin_b3.png'
	- Summary: planar composite of sensible heat flux at the maximum intensity bin for all storms.
	- Data: run `tc_analysis.tc_model_data(['AM2.5', 'FLOR', 'HIRAM'], ['control', 'swishe'], num_storms=50)
	- Figure generation: run `visualization.planar_composite(data, ['AM2.5', 'FLOR', 'HIRAM'], ['b3'], 'shflx', None,  ['control', 'swishe', 'swishe-control'])`
	- Statistics:
		- Number of records used for the composite for AM2.5, control, b3: 84
		- Number of records used for the composite for AM2.5, swishe, b3: 5
		- Number of records used for the composite for FLOR, control, b3: 322
		- Number of records used for the composite for FLOR, swishe, b3: 20
		- Number of records used for the composite for HIRAM, control, b3: 86
		- Number of records used for the composite for HIRAM, swishe, b3: 24

- `TC-field_lhflx-bin_b3.png'
	- Summary: planar composite of latent heat flux at the maximum intensity bin for all storms.
	- Data: run `tc_analysis.tc_model_data(['AM2.5', 'FLOR', 'HIRAM'], ['control', 'swishe'], num_storms=50)
	- Figure generation: run `visualization.planar_composite(data, ['AM2.5', 'FLOR', 'HIRAM'], ['b3'], 'shflx', None,  ['control', 'swishe', 'swishe-control'])`
	- Statistics:
		- Number of records used for the composite for AM2.5, control, b3: 84
		- Number of records used for the composite for AM2.5, swishe, b3: 5
		- Number of records used for the composite for FLOR, control, b3: 322
		- Number of records used for the composite for FLOR, swishe, b3: 20
		- Number of records used for the composite for HIRAM, control, b3: 86
		- Number of records used for the composite for HIRAM, swishe, b3: 24


