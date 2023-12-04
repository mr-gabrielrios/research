import numpy as np
import pandas as pd
import os
import xarray as xr

class TC:
    def __init__(self, storm_id, max_wind, min_slp):
        self.storm_id = storm_id
        self.max_wind = max_wind
        self.min_slp = min_slp