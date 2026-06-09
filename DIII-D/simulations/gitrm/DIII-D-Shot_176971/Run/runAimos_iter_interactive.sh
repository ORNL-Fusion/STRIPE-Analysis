bin=/gpfs/u/home/MPFS/MPFSnthd/barn/GITRm-MultiSpecies/build-dcs-rhel8-gcc84-gitrm/GITRm
td=/gpfs/u/home/MPFS/MPFSnthd/scratch/DIII-D-176971_data
prof_solps=profiles_solps_final_new.nc
prof_flux=DIII-D_flux.nc
solps_matrix=lu.nc
bfile=bField.nc
mesh=176971_2.osh

mpirun --bind-to core -np 1 $bin \
 --kokkos-ndevices=1 \
 ${td}/${mesh}\
 ${td}/gitrm_sep_1.ptn \
 ${td}/$prof_solps \
 ${td}/$prof_flux \
 ${td}/$bfile \
 ${td}/$solps_matrix \
 ${td}/ADAS_Rates_C.nc \
 ${td}/ADAS_Rates_W.nc \
 ${td}/ftridynSelf.nc \
 ${td}/surface_model_GITRm_rustbca_C_W_d.nc \
 -
rm /tmp/hosts.$SLURM_JOB_ID

