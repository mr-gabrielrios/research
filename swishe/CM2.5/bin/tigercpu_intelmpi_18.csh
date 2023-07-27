source /usr/share/Modules/init/csh

module purge
#module load intel/18.0/64/18.0.1.163
module load intel/18.0/64/18.0.3.222 # WY: 18.0.1.163 became unavailable on Jun 13, 2018
#module load intel-mpi/intel/2018.1/64
module load intel-mpi/intel/2018.3/64 # WY: 2018.1/64 became unavailable on Jun 13, 2018
module load hdf5/intel-16.0/intel-mpi/1.8.16
module load netcdf/intel-16.0/hdf5-1.8.16/intel-mpi/4.4.0
which mpif90
which mpicc
