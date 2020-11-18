!*************************************************************************************************
! Start up procedures for selected sites
! (read common parameters and climate data)
!*************************************************************************************************

!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> MY:rm
!   Include 'modules.f90'
!<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< MY:rm
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> MY:add
   Include 'mod_meteorology.f90'
!<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< MY:add
   Include 'main.f90'
   Include 'initialize.f90'
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> MY:rm
!   Include 'metabolic.f90'
!<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< MY:rm
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> MY:add
   Include 'metabolic_mod.f90'
!<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< MY:add
   Include 'output.f90'
   Include 'physics.f90'
   Include 'population_regu.f90'
   Include 'spatial_calc.f90'
   Include 'etc.f90'

PROGRAM start_point

!_____________ Set Variables
!Namespace
   USE data_structure
!!!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>TN:Add
   USE vegi_status_current1
   USE vegi_status_current2
   USE grid_status_current2
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> MY:rm
!   USE mod_grid
!<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< MY:rm
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> MY:add
   USE mod_meteorology
!<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< MY:add
!!!<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<TN:Add
   implicit none

!Climate data
   integer YearMaxClimate !year length of the climate data
   integer YearMaxCO2     !year length of the CO2 data

   real,allocatable,dimension(:,:)  ::&
   tmp_air       ,& !1. Surface air temperature (Celcius)
   tmp_soil1     ,& !2. Soil temperature   0- 10cm depth (Celcius)
   tmp_soil2     ,& !3. Soil temperature  10-200cm depth (Celcius)
   tmp_soil3     ,& !4. Soil temperature 200-300cm depth (Celcius)
   prec          ,& !5. Precipitation (mm day-1)
   cloud         ,& !6. Cloudness (fraction)
   wind          ,& !7. Wind velocity (m s-1)
   humid         ,& !8. Specific humidity (kg kg-1)
   tmp_air_range    !9. Daily range of tmp_air (Celcius)

   real,allocatable,dimension(:,:,:)::&
    tmp_soil      !Soil temperature for each layers (Celcius)

!Atomospheric CO2 time-series @ ppm
   real,allocatable,dimension(:)::aco2_annual

!Location data
   integer Mask         !Land ocean mask (1:land, 0:ocean)
   real    ALT          !altitude (m above MSL)
   real    Albedo_soil0 !albedo, default
   real    W_fi         !filed capacity   (m3/m3, 0.0 -> 1.0)
   real    W_wilt       !wilting point    (m3/m3, 0.0 -> 1.0)
   real    W_sat        !saturate point   (m3/m3, 0.0 -> 1.0)
   real    W_mat        !matrix potential (m, -0.0001 -> -3.0)

!Others
   real    LAT, LON       !latitude and logitude for simulate
   integer GlobalZone     !ID number of global zone
   integer point          !
   integer i, j, count    !for general usage
   real    x              !for general usage

!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> MY:add
   character(128) :: input_file_name(2)
   character*128  :: Fn_rad_out
   integer, parameter :: Fn_rad = 50
!<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< MY:add

!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> MY:moved to below
!!_____________ Set location
!!  LAT    north:+, south:- (decimalized)
!!  LON    east:+, west:-   (decimalized)
!!!!! Fukido@ Ishigaki Island, Japan <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<TN changed
!    LAT    =  24.48
!    LON    = 124.18
!
!!Set cordinate variables
!   point = (90-int(LAT)-1)*360 + int(LON+180) + 1
!
!!GlobalZone: Set Location category
!   if (LON>=-20 .and. 60>=LON .and. 35.0>=LAT) then
!      !African continent
!      GlobalZone = 1
!   elseif (LON>=100 .and. 170>=LON .and. 50.0<=LAT) then
!      !Eastern Siberia
!      GlobalZone = 2
!   else
!      !Default
!      GlobalZone = 0
!   endif
!<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< MY:moved to below

!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> MY:add
   ! Read input file names

   read(5,*) input_file_name
   write(*,*) 'input_file_name for SEIB-DGVM', input_file_name(1)
   write(*,*) 'input_file_name for coupled module', input_file_name(2)
!<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< MY:add

!_____________ Read Parameters
!Read Parameter files
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> MY:rm
!   open (1, file='./input/parameter_mangrove.txt', action='READ', status='OLD')
!<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< MY:rm
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> MY:add
   open (1, file = input_file_name(1), action='READ', status='OLD')
!<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< MY:add
      read ( unit=1, nml=Control)
      read ( unit=1, nml=PFT_type)
      read ( unit=1, nml=Respiration)
      read ( unit=1, nml=Turnover_n)
      read ( unit=1, nml=Metabolic)
      read ( unit=1, nml=Assimilation)
      read ( unit=1, nml=Dynamics)
      read ( unit=1, nml=Disturbance)
      read ( unit=1, nml=Soil_resp)
   close (1)
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> MY:add

    ! Check the inputs

    write(*,*) 'Simulation_year', Simulation_year
    write(*,*) 'Max_loc, STEP, C_in_drymass', Max_loc, STEP, C_in_drymass
    write(*,*) 'ALM5 Rh, Br', ALM5(1), ALM5(2)
    write(*,*) 'LA_max Rh, Br', LA_max(1), LA_max(2)
    write(*,*) 'DBH_limit Rh, Br', DBH_limit(1), DBH_limit(2)
    write(*,*) 'HGT_max Rh, Br', HGT_max(1), HGT_max(2)
    write(*,*) 'SLA Rh, Br', SLA(1), SLA(2)
    write(*,*) 'M4 Rh, Br', M4(1), M4(2)
    write(*,*) 'P_establish Rh, Br', P_establish(1), P_establish(2)
    write(*,*) 'PAR_min Rh, Br', PAR_min(1), PAR_min(2)
    write(*,*) 'AGE_max Rh, Br', AGE_max(1), AGE_max(2)
    write(*,*) 'Est_scenario, Est_year_change, Est_frac_random', Est_scenario, Est_year_change, Est_frac_random

    ! Read parameters for coupled module

    open (2, file = input_file_name(2), action='READ', status='OLD')
       read (unit = 2, nml = Calc_settings)
       read (unit = 2, nml = File_names)
       read (unit = 2, nml = Atmospheric_parameters)
       read (unit = 2, nml = Soil_parameters)
       read (unit = 2, nml = Waterflux_parameters)
       read (unit = 2, nml = PFT_parameters)
       read (unit = 2, nml = Root_parameters)
       read (unit = 2, nml = Stem_parameters)
       read (unit = 2, nml = Leaf_parameters)
       read (unit = 2, nml = Metabolic_parameters)
       read (unit = 2, nml = Mortality_parameters)
       read (unit = 2, nml = Allometric_parameters)
    close (2)

    ! Check the inputs

    write(*,*) 'growth_calc', growth_calc
    write(*,*) 'jra55_use_switch', jra55_use_switch
    write(*,*) 'flux_calculation', flux_calculation
    write(*,*) 'n_spe', n_spe
    write(*,*) 'mort_switch', mort_switch
    write(*,*) 'jra55_timezone', jra55_timezone
    write(*,*) 'jra55_lat', jra55_lat
    write(*,*) 'Fn_jra55', Fn_jra55
    write(*,*) 'Fn_forest_out', Fn_forest_out
    write(*,*) 'Fn_forest2_out', Fn_forest2_out
    write(*,*) 'Fn_monitor_tree_out', Fn_monitor_tree_out
    write(*,*) 'Fn_monitor_biomass_out', Fn_monitor_biomass_out
    write(*,*) 'Fn_monitor_plot_out', Fn_monitor_plot_out
    write(*,*) 'o2air, co2air', o2air, co2air
    write(*,*) 'ground_elev, soilresis', ground_elev, soilresis
    write(*,*) 'albsoi_vis, albsoi_nir', albsoi_vis, albsoi_nir
    write(*,*) 'hksat, moist, bsw', hksat, moist, bsw
    write(*,*) 'psisat, soil_t', psisat, soil_t
    write(*,*) 'sal_ini, din_ini, dip_ini', sal_ini, din_ini, dip_ini
    write(*,*) 'flux_fw, din_fw, dip_fw', flux_fw, din_fw, dip_fw
    write(*,*) 'sal_sw, din_sw, dip_sw', sal_sw, din_sw, dip_sw
    write(*,*) 'rhol_vis Rh, Br', rhol_vis(1), rhol_vis(2)
    write(*,*) 'rhol_nir Rh, Br', rhol_nir(1), rhol_nir(2)
    write(*,*) 'taul_vis Rh, Br', taul_vis(1), taul_vis(2)
    write(*,*) 'taul_nir Rh, Br', taul_nir(1), taul_nir(2)
    write(*,*) 'emleaf', emleaf
    write(*,*) 'xl Rh, Br', xl(1), xl(2)
    write(*,*) 'root_filter Rh, Br', root_filter(1), root_filter(2)
    write(*,*) 'root_resist Rh, Br', root_resist(1), root_resist(2)
    write(*,*) 'root_density Rh, Br', root_density(1), root_density(2)
    write(*,*) 'root_radius Rh, Br', root_radius(1), root_radius(2)
    write(*,*) 'root_depth', root_depth
    write(*,*) 'fine_root_ratio Rh, Br', fine_root_ratio(1), fine_root_ratio(2)
    write(*,*) 'k_sap Rh, Br', k_sap(1), k_sap(2)
    write(*,*) 'wood_rho Rh, Br', wood_rho(1), wood_rho(2)
    write(*,*) 'c_n_leaf Rh, Br', c_n_leaf(1), c_n_leaf(2)
    write(*,*) 'c_n_stem Rh, Br', c_n_stem(1), c_n_stem(2)
    write(*,*) 'c_n_root Rh, Br', c_n_root(1), c_n_root(2)
    write(*,*) 't_acclim Rh, Br', t_acclim(1), t_acclim(2)
    write(*,*) 't_growth Rh, Br', t_growth(1), t_growth(2)
    write(*,*) 't_home Rh, Br', t_home(1), t_home(2)
    write(*,*) 'optimum_dpai Rh, Br', optimum_dpai(1), optimum_dpai(2)
    write(*,*) 'minlp Rh, Br', minlp(1), minlp(2)
    write(*,*) 'dleaf Rh, Br', dleaf(1), dleaf(2)
    write(*,*) 'gs_max Rh, Br', gs_max(1), gs_max(2)
    write(*,*) 'vcmaxpft Rh, Br', vcmaxpft(1), vcmaxpft(2)
    write(*,*) 'iota Rh, Br', iota(1), iota(2)
    write(*,*) 'grow_resp Rh, Br', grow_resp(1), grow_resp(2)
    write(*,*) 'main_resp_stem Rh, Br', main_resp_stem(1), main_resp_stem(2)
    write(*,*) 'main_resp_root Rh, Br', main_resp_root(1), main_resp_root(2)
    write(*,*) 'root_turn Rh, Br', root_turn(1), root_turn(2)
    write(*,*) 'coarse_root_turn Rh, Br', coarse_root_turn(1), coarse_root_turn(2)
    write(*,*) 'leaf_turn Rh, Br', leaf_turn(1), leaf_turn(2)
    write(*,*) 'leaf_resorp Rh, Br', leaf_resorp(1), leaf_resorp(2)
    write(*,*) 'stock_trunk_ratio Rh, Br', stock_trunk_ratio(1), stock_trunk_ratio(2)
    write(*,*) 'mort_d1 Rh, Br', mort_d1(1), mort_d1(2)
    write(*,*) 'mort_c_d2 Rh, Br', mort_c_d2(1), mort_c_d2(2)
    write(*,*) 'mort_c_d3 Rh, Br', mort_c_d3(1), mort_c_d3(2)
    write(*,*) 'crown_a Rh, Br', crown_a(1), crown_a(2)
    write(*,*) 'crown_b Rh, Br', crown_b(1), crown_b(2)
    write(*,*) 'tree_h_a Rh, Br', tree_h_a(1), tree_h_a(2)
    write(*,*) 'tree_h_b Rh, Br', tree_h_b(1), tree_h_b(2)
    write(*,*) 'pr_s_a Rh, Br', pr_s_a(1), pr_s_a(2)
    write(*,*) 'pr_s_b Rh, Br', pr_s_b(1), pr_s_b(2)
    write(*,*) 'pr_h_a Rh, Br', pr_h_a(1), pr_h_a(2)
    write(*,*) 'pr_h_b Rh, Br', pr_h_b(1), pr_h_b(2)
    write(*,*) 'pr_d Rh, Br', pr_d(1), pr_d(2)

    ! Soil salinity (psu) -> (mol/m3)

    sal_ini = sal_ini * 1000.0d0 / 58.44d0
    sal_sw = sal_sw * 1000.0d0 / 58.44d0

    ! Give initial values for soil water variables

    sal = sal_ini
    din = din_ini
    dip = dip_ini

    ! Read tide (predicted by naotide): Jan1 2014 00:00 - Dec31 2014 23:00

    open (81, file = './input/tide.txt'   , status = 'old')
    read (81,*) time_series_tide
    close (81)

    ! Read meteorological data (JMA hourly-data or JRA-55 3hourly-data) : Jan1 2014 00:00 - Dec31 2014 23:00

    Fn_rad_out = './output/hourly_rad.csv'
    open ( Fn_rad, file = Fn_rad_out, status = 'replace')

    call meteorology_read ( &
                                         ! *** From SEIB-DGVM ***
    jra55_use_switch                , &  ! Switch of JRA-55 data (1: use, 0: no use)
    jra55_timezone                  , &  ! Time zone of JRA-55 data (correction for local time)
    jra55_lat                       , &  ! Latitude of JRA-55 data
    Fn_jra55                        , &  ! File name of JRA-55 data
    Fn_rad                            &  ! File number for hourly solar radiation
    )

!<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< MY:add

!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> MY:moved from above
!_____________ Set location
!  LAT    north:+, south:- (decimalized)
!  LON    east:+, west:-   (decimalized)
!!!! Fukido@ Ishigaki Island, Japan <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<TN changed
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> MY:rm
!    LAT    =  24.48
!<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< MY:rm
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> MY:rm
    LAT    = jra55_lat
!<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< MY:rm
    LON    = 124.18

!Set cordinate variables
   point = (90-int(LAT)-1)*360 + int(LON+180) + 1

!GlobalZone: Set Location category
   if (LON>=-20 .and. 60>=LON .and. 35.0>=LAT) then
      !African continent
      GlobalZone = 1
   elseif (LON>=100 .and. 170>=LON .and. 50.0<=LAT) then
      !Eastern Siberia
      GlobalZone = 2
   else
      !Default
      GlobalZone = 0
   endif
!<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< MY:moved from above
!!!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>TN:rm
!!Read Location data
!   open ( File_no(1), file=Fn_location, status='OLD')
!   do i=1, point
!      read(File_no(1),*) Mask, ALT, Albedo_soil0, W_sat, W_fi, W_mat, W_wilt
!   end do
!   close( File_no(1) )
!!!<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<TN:rm
!!!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>TN:add >>>>>>>>> MY:rm
!!Grid data
!   Call read_ROMS_files('./Fukido_grd_v6.nc'           &
!   &                   ,'./ocean_his7.nc'    &
!   &                   , 4, 253  )
!
    Mask    = 1         !Land ocean mask (1:land, 0:ocean)
    ALT     = 0.0       !altitude (m above MSL)
    Albedo_soil0 = 0.1  !albedo, default                      !!!�l�Ă��Ƃ�
    W_sat   = 0.4565    !saturate point   (m3/m3, 0.0 -> 1.0) !!!�l�Ă��Ƃ�
    W_fi    = 0.41708   !filed capacity   (m3/m3, 0.0 -> 1.0) !!!�l�Ă��Ƃ�
    W_mat   =-0.3       !matrix potential (m, -0.0001 -> -3.0)!!!�l�Ă��Ƃ�
    W_wilt  = 0.16642   !wilting point    (m3/m3, 0.0 -> 1.0) !!!�l�Ă��Ƃ�
!    write(*,*) GRID%N_x,GRID%N_y
!    write(*,*) GRID%Max_x,GRID%Max_y
!
!   ! Allocate variables for each floor cell
!   allocate (gmass_leaf     (GRID%N_x, GRID%N_y)     )
!   allocate (gmass_root     (GRID%N_x, GRID%N_y)     )
!   allocate (gmass_available(GRID%N_x, GRID%N_y)     )
!   allocate (gmass_stock    (GRID%N_x, GRID%N_y)     )
!   allocate (lai_grass      (GRID%N_x, GRID%N_y)     )
!   allocate (lai_opt_grass_RunningRecord(20, GRID%N_x, GRID%N_y) )
!
!   allocate (patch_vacant   (GRID%N_x, GRID%N_y)     )
!
!   allocate (par_floor_rel  (GRID%N_x, GRID%N_y)     )
!   allocate (par_grass_rel  (GRID%N_x, GRID%N_y)     )
!   allocate (sum_par_floor  (GRID%N_x, GRID%N_y)     )
!!!<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<TN:add <<<<<<<<< MY:rm
   !data processing
   if (W_fi   > W_sat ) W_fi   = W_sat
   if (W_wilt > W_sat ) W_wilt = W_sat

!_____________ Read Climate Data
!Scan climate data length
   count = 0
   open (File_no(1), file=Fn_climate, status='OLD')
   do i=1,Day_in_Year * 10000
      count = count + 1
      read(File_no(1),*,end=100)
   enddo
   100 continue
   close (File_no(1))
   YearMaxClimate = int(count/Day_in_Year)

!Set sizes of allocatable climate data table
   allocate (tmp_air       (Day_in_Year, YearMaxClimate)     )
   allocate (tmp_soil1     (Day_in_Year, YearMaxClimate)     )
   allocate (tmp_soil2     (Day_in_Year, YearMaxClimate)     )
   allocate (tmp_soil3     (Day_in_Year, YearMaxClimate)     )
   allocate (cloud         (Day_in_Year, YearMaxClimate)     )
   allocate (prec          (Day_in_Year, YearMaxClimate)     )
   allocate (humid         (Day_in_Year, YearMaxClimate)     )
   allocate (wind          (Day_in_Year, YearMaxClimate)     )
   allocate (tmp_air_range (Day_in_Year, YearMaxClimate)     )

   allocate (tmp_soil      (Day_in_Year, YearMaxClimate, NumSoil) )

!Read Climatic data
   open (File_no(1), file=Fn_climate, status='OLD')
   do i=1, YearMaxClimate
   do j=1, Day_in_Year
      read (File_no(1),*) tmp_air(j,i), tmp_soil1(j,i), tmp_soil2(j,i), tmp_soil3(j,i), &
                          prec(j,i), cloud(j,i), wind(j,i), humid(j,i), tmp_air_range(j,i)

      !Give adhock values for soil temperature
      Call tmp_soil_interpolate (tmp_soil1(j,i), tmp_soil2(j,i), tmp_soil3(j,i), tmp_soil(j,i,:))
   enddo
   enddo
   close (File_no(1))

   !______________ Read time-series of atmospheric CO2
!Scan climate data length
   YearMaxCO2 = 0
   open (File_no(1), file=Fn_CO2, status='OLD')
   do i=1, 10000
      YearMaxCO2 = YearMaxCO2 + 1
      read(File_no(1),*,end=200)
   enddo
   200 continue
   close (File_no(1))
   YearMaxCO2 = YearMaxCO2 - 1

!Set sizes of allocatable CO2 data table
   allocate ( aco2_annual(YearMaxCO2) )

!Read CO2 data
   Open (File_no(1), file=Fn_CO2, status='OLD')
   do i = 1, YearMaxCO2
      read(File_no(1), *) aco2_annual(i)
   end do
   Close (File_no(1))

   !___________ For employing different random seed for each run (by Shigeki IKEDA @ Kyoto Univ.)
   IF (Flag_randomization) then
      call random_seed(size=seedsize)
      allocate(seed(seedsize))

      do size_count=1,seedsize
      call system_clock(count=clock)
      seed(size_count)=clock
      end do

      call random_seed(put=seed)
   EndIf

!_____________ Display Simulation Conditions
   !Print location properties
   write (*,*)
   write (*,*) '*********  Coodinate Configurations  *********'
   write (*,*) 'Latitude, Longtitude   :', LAT, LON
   write (*,*) 'Point nomal            :', point
   write (*,*)

   !Print soil properties
   write (*,*) '*********  Location properties  *********'
   write (*,*) 'Altitude    :', ALT
   write (*,*) 'Albedo_soil0:', Albedo_soil0
   write (*,*) 'W_sat       :', W_sat
   write (*,*) 'W_fi        :', W_fi
   write (*,*) 'W_wilt      :', W_wilt
   write (*,*)

   !Print climate properties
   x = real(Day_in_Year * YearMaxClimate)

   write (*,*) '*********  Wether statistics (annual mean of the all years)  *********'
   write (*,*) '2m air temperature  (Cecius)  :', sum(tmp_air  (:,:)) / x
   write (*,*) 'precipitation       (mm/year) :', sum(prec     (:,:)) / real(YearMaxClimate)
   write (*,*) 'humidity            (kg/kg)   :', sum(humid    (:,:)) / x
   write (*,*) 'wind                (m/s)     :', sum(wind     (:,:)) / x
   write (*,*)

!_____________ Simulation
   !Call simulation loop
   IF (Mask==0 .or. Albedo_soil0<=0.0 .or. W_sat<=0.0 .or. W_fi<=0.0 .or. W_wilt<=0.0 ) then
      write(*,*) 'Error: Invalid location properties'
   ELSE
      write (*,*) '*********  Now simulating  *********'
      Call main_loop ( &
      LAT, LON, GlobalZone, YearMaxClimate, YearMaxCO2, &
      tmp_air(:,:), prec(:,:), cloud(:,:), wind(:,:), humid(:,:), &
      tmp_air_range(:,:), tmp_soil(:,:,:), &
      aco2_annual, ALT, Albedo_soil0, W_fi, W_wilt, W_sat, W_mat )
      write (*,*) '*********  Done  *********'
   END IF

   STOP

END PROGRAM start_point
