## SWISHE small-scale manuscript - figure guide
This `readme` is intended to provide a summary of what each figure is and how to reproduce it.

- `TC_density_C15w.png`
	- Summary: Plot of number of TC days per N-degree x N-degree bin for each model over the first 50 years of each model run.
	- Data: run `tc_analysis.tc_track_data(models=['AM2.5', 'FLOR', 'HIRAM'], experiments=['control', 'swishe'], storm_type='C15w', year_range=(101, 150))`
	- Figure generation: run `visualization.density_grid(data, model_names=['AM2.5', 'FLOR', 'HIRAM'], bin_resolution=5)`

