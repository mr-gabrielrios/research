import cdsapi
import os
import datetime
import argparse

# Parse input arguments
parser = argparse.ArgumentParser()

parser.add_argument('--varname')
parser.add_argument('--year_range')
parser.add_argument('--month_range', nargs='?', const='1:12')

args = parser.parse_args()
varname = args.varname
year_range = args.year_range
month_range = args.month_range

# Format arguments
year_range = [int(year) for year in year_range.split(':')]
month_range = [int(month) for month in month_range.split(':')] if month_range else [1, 12]

# Perform checks
assert isinstance(min(year_range), int) & isinstance(max(year_range), int), 'Year range entries must both be integers.'
# Generate list of strings from year range
year_range_str = [str(year) for year in range(min(year_range), max(year_range))]

assert isinstance(min(month_range), int) & isinstance(max(month_range), int), 'Month range entries must both be integers.'
# Generate list of strings from month range
month_range_str = [str(month) for month in range(min(month_range), max(month_range))]

# Generate list of strings for day range
day_range_str = [f'{day:02d}' for day in range(1, 32)]

dataset = "reanalysis-era5-single-levels"
request = {
    "product_type": ["reanalysis"],
    "variable": [varname],
    "year": year_range_str,
    "month": month_range_str,
    "day": day_range_str,
    "time": ["00:00", "06:00", "12:00", "18:00"],
    "data_format": "netcdf",
    "download_format": "unarchived"
}

print(request)

# Save data
dirname = "/scratch/gpfs/GEOCLIM/gr7610/tiger3/reference/datasets"
filename = f"ERA5_request.varname-{varname}.year-{min(year_range)}_{max(year_range)}.month-{min(month_range)}_{max(month_range)}.{datetime.datetime.now().strftime('%Y%m%d%H%M')}.nc"
target = os.path.join(dirname, filename)

client = cdsapi.Client()
client.retrieve(dataset, request, target)