NP_CALC = 4
NP_UPDATE = 4

enkf: enkf-prep enkf-calc enkf-update

enoi: clean enoi-prep enoi-calc enoi-update

stats: clean stats-prep stats-calc

singleob: clean singleob-calc singleob-update

pointlogs: pointlogs-calc pointlogs-update

global: clean enkf-prep enkf-calc-global enkf-update-global

enkf-prep: enkf_prep obs enkf.prm
	./enkf_prep --no-superobing enkf.prm 2>&1 | tee prep.out

enkf-calc: enkf_calc observations.nc prior
	nice -n +19 mpirun -np $(NP_CALC) ./enkf_calc --use-rmsd-for-obsstats --ignore-no-obs enkf.prm 2>&1 | tee calc.out

enkf-update: enkf_update X5.nc prior
	nice -n +19 mpirun -np $(NP_UPDATE) ./enkf_update --calculate-spread enkf.prm 2>&1 | tee update.out

enoi-prep: enkf_prep conf obs enoi.prm
	./enkf_prep enoi.prm 2>&1 | tee prep.out

enoi-calc: enkf_calc observations.nc conf background_6565 ensemble_static
	nice -n +19 mpirun -np $(NP_CALC) ./enkf_calc enoi.prm 2>&1 | tee calc.out

enoi-update: enkf_update w.nc conf background_6565 ensemble_static
	nice -n +19 mpirun -np $(NP_UPDATE) ./enkf_update --calculate-spread --output-increment enoi.prm 2>&1 | tee update.out

stats-prep: enkf_prep conf obs stats.prm
	./enkf_prep stats.prm 2>&1 | tee prep.out

stats-calc: enkf_calc observations.nc conf background_6565
	nice -n +19 mpirun -np $(NP_CALC) ./enkf_calc --forecast-stats-only stats.prm 2>&1 | tee calc.out

singleob-calc: clean enkf_calc prior
	nice -n +19 mpirun -np $(NP_CALC) ./enkf_calc --single-observation-ijk 350 360 0 sea_ice_thickness 2 .003 singleob.prm 2>&1 | tee calc.out

#singleob-calc: clean enkf_calc prior
#	nice -n +19 mpirun -np $(NP_CALC) ./enkf_calc --single-observation-xyz .5 .5 0 sea_ice_thickness 2 .003 singleob.prm 2>&1 | tee calc.out

singleob-update: enkf_update prior
	nice -n +19 mpirun -np $(NP_UPDATE) ./enkf_update --output-increment singleob.prm 2>&1 | tee update.out

pointlogs-calc: enkf_calc observations.nc conf prior X5.nc
	cp enkf.prm pointlogs.prm
	echo "POINTLOG 50 50" >> pointlogs.prm
	echo "POINTLOG 70 50" >> pointlogs.prm
	nice -n +19 mpirun -np $(NP_CALC) ./enkf_calc --point-logs-only pointlogs.prm 2>&1 | tee calc-pointlogs.out

pointlogs-update: enkf_update X5.nc conf prior
	nice -n +19 mpirun -np $(NP_UPDATE) ./enkf_update --no-fields-write pointlogs.prm 2>&1 | tee update-pointlogs.out

enkf-calc-global: enkf_calc observations.nc conf prior
	nice -n +19 ./enkf_calc enkf-global.prm 2>&1 | tee calc.out

enkf-update-global: enkf_update X5.nc conf prior
	nice -n +19 mpirun -np $(NP_UPDATE) ./enkf_update --calculate-spread --write-inflation enkf-global.prm 2>&1 | tee update.out
	rm -f prior/mem*_*-???.nc

#enkf_prep:
#	cp ../../bin/enkf_prep .
#
#enkf_calc:
#	cp ../../bin/enkf_calc .
#
#enkf_update:
#	cp ../../bin/enkf_update .

clean:
	@rm -f X5*.nc w*.nc enkf_diag*.nc badbatches.txt observations*.nc *.out prior/*.analysis* prior/*.increment* prior/*-???.nc background_6565/*.analysis* background_6565/*.increment* prior/ens_*-???.nc background_6565/bg_*-???.nc pointlogs.prm
	@rm -rf .enkftmp
