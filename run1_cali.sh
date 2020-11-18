#!/bin/bash

case_set=cali
case_name=growth_test_aklan

# Set soil salinity
#i1=20
#i2=25
#i3=30
i4=25

rm seib.exe
rm ./input_case_${case_set}_*
rm ./input/parameter_coupled_module_${case_set}_*

SRC_DIR1=./seib_coupled_code
SRC_DIR2=./stand_level_coupled_code

gfortran -fbounds-check ${SRC_DIR1}/modules.f90 ${SRC_DIR2}/mod_param.f90 ${SRC_DIR2}/mod_water_vapor.f90 ${SRC_DIR2}/mod_math_tools.f90 ${SRC_DIR2}/mod_soil_char.f90 ${SRC_DIR2}/mod_plant_hydraulics.f90 ${SRC_DIR2}/mod_nitrogen_profile.f90 ${SRC_DIR2}/mod_photosynthesis.f90 ${SRC_DIR2}/mod_leaf_boundary_layer.f90 ${SRC_DIR2}/mod_energy_water_balance.f90 ${SRC_DIR2}/mod_stomatal_conductance.f90 ${SRC_DIR2}/mod_tree_allometry.f90 ${SRC_DIR2}/mod_spac_photosynthesis.f90 ${SRC_DIR2}/mod_metabolism.f90 ${SRC_DIR2}/mod_crown_morphology.f90 ${SRC_DIR2}/mod_growth.f90 ${SRC_DIR1}/mod_radiation.f90 ${SRC_DIR1}/mod_soil_water_flux.f90 ${SRC_DIR1}/mod_monitoring.f90 ${SRC_DIR1}/start_point.f90 -O2 -o seib.exe

rm *.mod

mkdir ./output/${case_set}_${case_name}
cp ./input/parameter_coupled_module_${case_set}.txt ./output/${case_set}_${case_name}/parameter_coupled_module_${case_set}.txt
cp ./input/parameter_mangrove.txt ./output/${case_set}_${case_name}/parameter_mangrove.txt

#sed -e "s/SALT/$i1/g" -e "s/CALIB/$case_name/g" -e "s/CASE/$case_set/g" ./input/parameter_coupled_module_${case_set}.txt > ./input/parameter_coupled_module_${case_set}_${i1}.txt
#sed -e "s/SALT/$i1/g" ./input_case_${case_set}.in > ./input_case_${case_set}_${i1}.in
#sed -e "s/SALT/$i2/g" -e "s/CALIB/$case_name/g" -e "s/CASE/$case_set/g" ./input/parameter_coupled_module_${case_set}.txt > ./input/parameter_coupled_module_${case_set}_${i2}.txt
#sed -e "s/SALT/$i2/g" ./input_case_${case_set}.in > ./input_case_${case_set}_${i2}.in
#sed -e "s/SALT/$i3/g" -e "s/CALIB/$case_name/g" -e "s/CASE/$case_set/g" ./input/parameter_coupled_module_${case_set}.txt > ./input/parameter_coupled_module_${case_set}_${i3}.txt
#sed -e "s/SALT/$i3/g" ./input_case_${case_set}.in > ./input_case_${case_set}_${i3}.in
sed -e "s/SALT/$i4/g" -e "s/CALIB/$case_name/g" -e "s/CASE/$case_set/g" ./input/parameter_coupled_module_${case_set}.txt > ./input/parameter_coupled_module_${case_set}_${i4}.txt
sed -e "s/SALT/$i4/g" ./input_case_${case_set}.in > ./input_case_${case_set}_${i4}.in

#./seib.exe < input_case_${case_set}_${i1}.in &
#./seib.exe < input_case_${case_set}_${i2}.in &
#./seib.exe < input_case_${case_set}_${i3}.in &
#./seib.exe < input_case_${case_set}_${i4}.in &
#wait

./seib.exe < input_case_${case_set}_${i4}.in

rm ./input_case_${case_set}_*
rm ./input/parameter_coupled_module_${case_set}_*

