module mod_monitoring

  !-----------------------------------------------------------------------
  ! !DESCRIPTION:
  !
  !
  ! !USES:
  !
  !
  ! !PUBLIC TYPES:
  implicit none
  !
  ! !PUBLIC MEMBER FUNCTIONS:
  public :: monitor_index
  public :: monitor_tree
  public :: monitor_plot
  !-----------------------------------------------------------------------

contains

  !-----------------------------------------------------------------------
  subroutine monitor_index ( &
                                       ! *** From SEIB-DGVM ***
    flag_tree_presence            , &  ! Flag of tree presence (0: No trees, 1: There are trees)
    tree_exist                    , &  ! Flag of tree presence
                                       !
                                       ! *** In/Output ***
    monitor                         &  ! Tree index for monitoring (0: No tree decided)
    )
    !
    ! !DESCRIPTION:
    ! Decide a tree to monitor.
    !
    ! !USES:
    use data_structure, only : Max_no, randf
    !
    ! !ARGUMENTS:
    implicit none
    integer, intent(in)    :: flag_tree_presence
    logical, intent(in)    :: tree_exist(:)
    integer, intent(inout) :: monitor
    !
    ! !LOCAL VARIABLES:
    integer, dimension(Max_no) :: tree_matrix  ! Matrix of tree index
    integer :: tree_n                          ! Number of trees
    integer :: no                              ! Tree index
    real :: rand_n                             ! Random value (0.0 - 1.0)
    integer :: loop                            ! Number of loops
    !---------------------------------------------------------------------

    ! Return if there are no trees.

    if (flag_tree_presence == 0) then
       return
    end if

    ! Check if the monitoring tree is still alive.
    ! If not, set a new tree for monitoring.

    if (monitor > 0) then
       if (.not. tree_exist(monitor)) then
          monitor = 0
       end if
    end if

    ! Return if the monitoring tree is decided.

    if (monitor > 0) then
       return
    end if

    !---------------------------------------------------------------------
    ! Making matrix of tree index
    !---------------------------------------------------------------------

    ! Zero out

    tree_matrix(:) = 0
    tree_n = 0

    do no = 1, Max_no

       if ( .not. tree_exist(no) ) cycle

       tree_n = tree_n + 1
       tree_matrix(tree_n) = no

    end do

    !---------------------------------------------------------------------
    ! Decide a tree to monitor
    !---------------------------------------------------------------------

    do loop = 1, 100
       do no = 1, tree_n
          rand_n = randf() ! 0.0 - 1.0
          if (rand_n < (1/real(tree_n))) exit
       end do
       if (rand_n < (1/real(tree_n))) exit
    end do

    monitor = tree_matrix(no)

    ! Error check
    if ( .not. tree_exist(monitor)) then
       write(*,*) 'ERROR: Monitoring tree is not existing.'
    end if

  end subroutine monitor_index

  !-----------------------------------------------------------------------
  subroutine monitor_tree ( &
                                       ! *** Input ***
                                       !
    time_series_rad_dir           , &  ! Time-series direct radiation at canopy top (W/m2)
    time_series_rad_dif           , &  ! Time-series diffused radiation at canopy top (W/m2)
                                       !
                                       ! *** From SEIB-DGVM ***
    Fn                            , &  ! File number for monitoring tree
    STEP                          , &  ! Canopy layer thickness (m)
    year                          , &  ! Year
    day_of_year                   , &  ! Day of the year (1 - 365)
    pft                           , &  ! Species index (1: Rh, 2: Br)
    bottom_layer_monitor          , &  ! Bottom layer of canopy
    top_layer_monitor             , &  ! Top layer of canopy
    par_direct_rel                , &  ! Profile of relative intensity of direct PAR within canopy compared to canopy top (fraction)
    par_diffuse_rel               , &  ! Profile of relative intensity of diffused PAR within canopy compared to canopy top (fraction)
    gs_monitor                    , &  ! Midday stand-level stomatal conductance (mol H2O/m2/s)
    an_top_monitor                , &  ! Midday top layer leaf net photosynthesis (umol CO2/m2 leaf/s)
    an_bot_monitor                , &  ! Midday bottom layer leaf net photosynthesis (umol CO2/m2 leaf/s)
    et_top_monitor                , &  ! Midday top layer leaf transpiration rate (mol H2O/m2 leaf/s)
    et_bot_monitor                , &  ! Midday bottom layer leaf transpiration rate (mol H2O/m2 leaf/s)
    c_uptake_bottom_day_monitor   , &  ! Daily bottom layer carbon uptake rate by photosynthesis (mol CO2/m2 leaf/day)
    n_uptake_bottom_day_monitor   , &  ! Daily bottom layer nitrogen uptake rate (mol N/m2 leaf/day)
    max_water_uptake_day_monitor  , &  ! Daily maximum stand-level sapflow rate (mol H2O/tree/day)
    max_stand_et_day_monitor      , &  ! Daily maximum transpiration rate at stand-level (mol H2O/tree/day)
    water_uptake_day_monitor      , &  ! Daily plant water uptake rate at stand-level (m3 H2O/tree/day)
    available_c_monitor           , &  ! Remaining carbon for tree growth after respiration and turnover (mol C/tree/day)
    available_n_monitor           , &  ! Remaining nitrogen for tree growth after turnover (mol N/tree/day)
    height_max_monitor            , &  ! Potential tree height based on allometric relation (m)
    height_limit                  , &  ! Tree height limitation based on proximate trees (unit is STEP!!)
    crown_d_max_monitor           , &  ! Potential crown diameter based on allometric relation (m)
    radius_limit                  , &  ! Crown radius limitation based on proximate trees (m)
    dbh_heartwood                 , &  ! Heartwood diameter (m)
    dbh_sapwood                   , &  ! Sapwood diameter (m)
    tree_h                        , &  ! Tree height (m)
    crown_diameter                , &  ! Crown diameter (m)
    seib_height                   , &  ! SEIB-DGVM specific tree height (unit is STEP!!)
    seib_bole                     , &  ! SEIB-DGVM specific bole height (unit is STEP!!)
    lai_monitor                   , &  ! Leaf area index (m2 leaf/m2 ground)
    mass_trunk                    , &  ! Stand-level trunk biomass (g trunk/tree)
    mass_coarse_root              , &  ! Stand-level coarse root biomass (g/tree)
    mass_above_root               , &  ! Stand-level above-ground root biomass (g/tree)
    mass_leaf                     , &  ! Stand-level leaf biomass (g leaf/tree)
    mass_root                     , &  ! Stand-level fine root biomass (g root/tree)
    mass_stock                      &  ! Stand-level stock biomass (g stock/tree)
    )
    !
    ! !DESCRIPTION:
    ! Daily output for monitoring tree
    !
    ! !USES:
    !
    !
    ! !ARGUMENTS:
    implicit none
    real(8), intent(in) :: time_series_rad_dir(:), time_series_rad_dif(:)
    integer, intent(in) :: Fn
    real, intent(in)    :: STEP
    integer, intent(in) :: year, day_of_year, pft, bottom_layer_monitor, top_layer_monitor
    real, intent(in)    :: par_direct_rel(:), par_diffuse_rel(:)
    real(8), intent(in) :: gs_monitor, an_top_monitor, an_bot_monitor
    real(8), intent(in) :: et_top_monitor, et_bot_monitor
    real(8), intent(in) :: c_uptake_bottom_day_monitor, n_uptake_bottom_day_monitor
    real(8), intent(in) :: max_water_uptake_day_monitor, max_stand_et_day_monitor
    real(8), intent(in) :: water_uptake_day_monitor
    real(8), intent(in) :: available_c_monitor, available_n_monitor
    real(8), intent(in) :: height_max_monitor
    integer, intent(in) :: height_limit
    real(8), intent(in) :: crown_d_max_monitor
    real, intent(in)    :: radius_limit, dbh_heartwood, dbh_sapwood
    real(8), intent(in) :: tree_h
    real, intent(in)    :: crown_diameter
    integer, intent(in) :: seib_height, seib_bole
    real(8), intent(in) :: lai_monitor
    real, intent(in)    :: mass_trunk, mass_coarse_root, mass_above_root
    real, intent(in)    :: mass_leaf, mass_root, mass_stock
    !
    ! !LOCAL VARIABLES:
    real(8), parameter :: mmh2o = 18.02d0 / 1000.0d0  ! Molecular mass of water (kg/mol)is:starred
    real(8), parameter :: denh2o = 1000.0d0           ! Density of liquid water (kg/m3)
    integer :: hour_of_year
    real(8) :: par_midday                             ! PAR at canopy top at midday (umol photon/m2 ground/s)
    real(8) :: par_top                                ! PAR at crown top layer at midday (umol photon/m2 ground/s)
    real(8) :: par_bot                                ! PAR at crown bottom layer at midday (umol photon/m2 ground/s)
    real(8) :: max_water_uptake_day_monitor_kg        ! Daily maximum stand-level water uptake potential (kg H2O/tree/day)
    real(8) :: max_stand_et_day_monitor_kg            ! Daily maximum transpiration rate at stand-level (kg H2O/tree/day)
    real(8) :: water_uptake_day_monitor_kg            ! Daily plant water uptake rate at stand-level (kg H2O/tree/day)
    real(8) :: height_limit_m                         ! Tree height limitation based on proximate trees (m)
    real(8) :: crown_dep                              ! Crown depth (m)
    real(8) :: dpai                                   ! Layer leaf area index (m2 leaf/m2 ground)
    !---------------------------------------------------------------------

    ! PAR at canopy top at midday (umol photon/m2 ground/s)

    hour_of_year = (day_of_year -1)*24 + 12 + 1
    par_midday = 4.6d0 * 0.43d0 * time_series_rad_dir(hour_of_year) &
                 + 4.2d0 * 0.57d0 * time_series_rad_dif(hour_of_year)

    ! PAR at crown top and bottom layer at midday (umol photon/m2 ground/s)

    par_top = par_direct_rel(top_layer_monitor) * 4.6d0 * 0.43d0 * time_series_rad_dir(hour_of_year) &
              + par_diffuse_rel(top_layer_monitor) * 4.2d0 * 0.57d0 * time_series_rad_dif(hour_of_year)

    par_bot = par_direct_rel(bottom_layer_monitor) * 4.6d0 * 0.43d0 * time_series_rad_dir(hour_of_year) &
              + par_diffuse_rel(bottom_layer_monitor) * 4.2d0 * 0.57d0 * time_series_rad_dif(hour_of_year)

    ! mol H2O/tree/day -> kg H2O/tree/day

    max_water_uptake_day_monitor_kg = max_water_uptake_day_monitor * mmh2o
    max_stand_et_day_monitor_kg = max_stand_et_day_monitor * mmh2o

    ! m3 H2O/tree/day -> kg H2O/tree/day

    water_uptake_day_monitor_kg = water_uptake_day_monitor * denh2o

    ! STEP -> m

    height_limit_m = real(height_limit) * STEP + 1.3d0

    ! Crown depth (m)

    crown_dep = real(seib_height - seib_bole) * STEP

    ! Layer leaf area index (m2 leaf/m2 ground)

    dpai = lai_monitor / real(seib_height - seib_bole)

    ! Write monitoring tree output

    write (Fn, '(3(i4,a), 31(f12.5,a))') &
    year                            , ',', &  ! Year
    day_of_year                     , ',', &  ! Day of the year (1 - 365)
    pft                             , ',', &  ! Species index (1: Rh, 2: Br)
    par_midday                      , ',', &  ! PAR at canopy top at midday (umol photon/m2 ground/s)
    par_top                         , ',', &  ! PAR at crown top layer at midday (umol photon/m2 ground/s)
    par_bot                         , ',', &  ! PAR at crown bottom layer at midday (umol photon/m2 ground/s)
    gs_monitor                      , ',', &  ! Midday stand-level stomatal conductance (mol H2O/m2/s)
    an_top_monitor                  , ',', &  ! Midday top layer leaf net photosynthesis (umol CO2/m2 leaf/s)
    an_bot_monitor                  , ',', &  ! Midday bottom layer leaf net photosynthesis (umol CO2/m2 leaf/s)
    et_top_monitor                  , ',', &  ! Midday top layer leaf transpiration rate (mol H2O/m2 leaf/s)
    et_bot_monitor                  , ',', &  ! Midday bottom layer leaf transpiration rate (mol H2O/m2 leaf/s)
    c_uptake_bottom_day_monitor     , ',', &  ! Daily bottom layer carbon uptake rate by photosynthesis (mol CO2/m2 leaf/day)
    n_uptake_bottom_day_monitor     , ',', &  ! Daily bottom layer nitrogen uptake rate (mol N/m2 leaf/day)
    max_water_uptake_day_monitor_kg , ',', &  ! Daily maximum stand-level water uptake potential (kg H2O/tree/day)
    max_stand_et_day_monitor_kg     , ',', &  ! Daily maximum transpiration rate at stand-level (kg H2O/tree/day)
    water_uptake_day_monitor_kg     , ',', &  ! Daily plant water uptake rate at stand-level (kg H2O/tree/day)
    available_c_monitor             , ',', &  ! Remaining carbon for tree growth after respiration and turnover (mol C/tree/day)
    available_n_monitor             , ',', &  ! Remaining nitrogen for tree growth after turnover (mol N/tree/day)
    height_max_monitor              , ',', &  ! Potential tree height based on allometric relation (m)
    height_limit_m                  , ',', &  ! Tree height limitation based on proximate trees (m)
    crown_d_max_monitor             , ',', &  ! Potential crown diameter based on allometric relation (m)
    radius_limit * 2.0d0            , ',', &  ! Crown diameter limitation based on proximate trees (m)
    dbh_heartwood + dbh_sapwood     , ',', &  ! DBH (m)
    tree_h                          , ',', &  ! Tree height (m)
    crown_diameter                  , ',', &  ! Crown diameter (m)
    crown_dep                       , ',', &  ! Crown depth (m)
    lai_monitor                     , ',', &  ! Leaf area index (m2 leaf/m2 ground)
    dpai                            , ',', &  ! Layer leaf area index (m2 leaf/m2 ground)
    mass_trunk / 1000.0d0           , ',', &  ! Stand-level trunk biomass (kg trunk/tree)
    mass_coarse_root / 1000.0d0     , ',', &  ! Stand-level coarse root biomass (g/tree)
    mass_above_root / 1000.0d0      , ',', &  ! Stand-level above-ground root biomass (g/tree)
    mass_leaf / 1000.0d0            , ',', &  ! Stand-level leaf biomass (kg leaf/tree)
    mass_root / 1000.0d0            , ',', &  ! Stand-level fine root biomass (kg root/tree)
    mass_stock / 1000.0d0                    ! Stand-level stock biomass (kg stock/tree)

  end subroutine monitor_tree

  !-----------------------------------------------------------------------
  subroutine monitor_plot ( &
                                       ! *** Input ***
                                       !
    time_series_rad_dir           , &  ! Time-series direct radiation at canopy top (W/m2)
    time_series_rad_dif           , &  ! Time-series diffused radiation at canopy top (W/m2)
                                       !
                                       ! *** From SEIB-DGVM ***
    Fn1                           , &  ! File number for tree biomass output
    Fn2                           , &  ! File number for plot output
    year                          , &  ! Year
    day_of_year                   , &  ! Day of the year (1 - 365)
    Max_loc                       , &  ! Dimension of the virtual forest (m)
    Max_hgt                       , &  ! Number of canopy layer
    Max_no                        , &  ! Maximum number of individual stands
    tree_exist                    , &  ! Flag of tree presence
    pft                           , &  ! Species index (1: Rh, 2: Br)
    age                           , &  ! Tree age (year: 1~)
    tree_h                        , &  ! Tree height (m)
    dbh_heartwood                 , &  ! Heartwood diameter (m)
    dbh_sapwood                   , &  ! Sapwood diameter (m)
    mass_trunk                    , &  ! Stand-level trunk biomass (g trunk/tree)
    mass_coarse_root              , &  ! Stand-level coarse root biomass (g/tree)
    mass_above_root               , &  ! Stand-level above-ground root biomass (g/tree)
    mass_leaf                     , &  ! Stand-level leaf biomass (g leaf/tree)
    mass_root                     , &  ! Stand-level fine root biomass (g root/tree)
    dpai_layer_sum                , &  ! Sum of layer leaf area index for each forest layer (m2 leaf/me ground)
    irveg                         , &  ! Absorbed longwave radiation, vegetation (W/m2)
    irsoi                         , &  ! Absorbed longwave radiation, ground (W/m2)
    ir_tree                       , &  ! VIS-interruption-coefficient by tree corwn (0.0 - 1.0). lower => more absorptance.
    rnsoi_day                     , &  ! Daily mean net radiation, ground (W/m2)
    rnveg_day                     , &  ! Daily mean net radiation, vegetation (W/m2)
    etsoi_day                     , &  ! Daily soil evaporation rate (mm/day)
    etveg_day                     , &  ! Daily transpiration rate (mm/day)
    sal                           , &  ! Pore-water salinity (mol/m3)
    din                           , &  ! DIN concentration in pore-water (mol N/m3)
    dip                             &  ! DIP concentration in pore-water (mol P/m3)
    )
    !
    ! !DESCRIPTION:
    ! Daily output for plot-scale variables
    !
    ! !USES:
    use data_structure, only : n_spe
    !
    ! !ARGUMENTS:
    implicit none
    real(8), intent(in) :: time_series_rad_dir(:), time_series_rad_dif(:)
    integer, intent(in) :: Fn1, Fn2, year, day_of_year, Max_loc, Max_hgt, Max_no
    logical, intent(in) :: tree_exist(:)
    integer, intent(in) :: pft(:), age(:)
    real(8), intent(in) :: tree_h(:)
    real, intent(in)    :: dbh_heartwood(:), dbh_sapwood(:), mass_trunk(:)
    real, intent(in)    :: mass_coarse_root(:), mass_above_root(:), mass_leaf(:), mass_root(:)
    real(8), intent(in) :: dpai_layer_sum(:)
    real(8), intent(in) :: irveg, irsoi
    real, intent(in)    :: ir_tree
    real(8), intent(in) :: rnsoi_day, rnveg_day, etsoi_day, etveg_day
    real(8), intent(in) :: sal, din, dip
    !
    ! !LOCAL VARIABLES:
    integer :: p                                        ! Species index (1: Rh, 2: Br)
    integer :: no                                       ! Tree index
    integer :: i                                        ! Layer index
    real, dimension(n_spe) :: tree_density_plot         ! Tree density for each species (tree/m2)
    real(8), dimension(n_spe) :: mean_tree_h            ! Mean tree height for each species (m)
    real, dimension(n_spe) :: mean_dbh                  ! Mean DBH for each species (m)
    real, dimension(n_spe) :: above_biomass_plot        ! Above-ground biomass for each species (Mg/ha)
    real, dimension(n_spe) :: leaf_biomass_plot         ! Leaf biomass for each species (Mg/ha)
    real, dimension(n_spe) :: root_biomass_plot         ! Root biomass for each species (Mg/ha)
    real, dimension(n_spe) :: coarse_root_biomass_plot  ! Coarse root biomass for each species (Mg/ha)
    real, dimension(n_spe) :: above_root_biomass_plot   ! Above-ground root biomass for each species (Mg/ha)
    real(8) :: lai_plot                                 ! Plot-level LAI (m2 leaf/m2 ground)
    !---------------------------------------------------------------------

    !---------------------------------------------------------------------
    ! Tree biomass output
    !---------------------------------------------------------------------

    ! Zero out

    tree_density_plot(:) = 0.0d0
    mean_tree_h(:) = 0.0d0
    mean_dbh(:) = 0.0d0
    above_biomass_plot(:) = 0.0d0
    leaf_biomass_plot(:) = 0.0d0
    root_biomass_plot(:) = 0.0d0
    coarse_root_biomass_plot(:) = 0.0d0
    above_root_biomass_plot(:) = 0.0d0

    ! Summing for all trees

    do no = 1, Max_no

       if ( .not. tree_exist(no) ) cycle
       if (age(no) <= 2) cycle  ! Ignore very small trees whose survival is not affected by tree density and production rate.
       if (dbh_heartwood(no) + dbh_sapwood(no) <= 0.02) cycle  ! Ignore very small trees

       p = pft(no)
       tree_density_plot(p) = tree_density_plot(p) + 1.0d0                                        ! tree number
       mean_tree_h(p) = mean_tree_h(p) + tree_h(no)                                               ! m
       mean_dbh(p) = mean_dbh(p) + dbh_heartwood(no) + dbh_sapwood(no)                            ! m
       above_biomass_plot(p) = above_biomass_plot(p) + mass_trunk(no) *1.e-03                     ! kg
       leaf_biomass_plot(p) = leaf_biomass_plot(p) + mass_leaf(no) *1.e-03                        ! kg
       root_biomass_plot(p) = root_biomass_plot(p) + mass_root(no) *1.e-03                        ! kg
       coarse_root_biomass_plot(p) = coarse_root_biomass_plot(p) + mass_coarse_root(no) *1.e-03   ! kg
       above_root_biomass_plot(p) = above_root_biomass_plot(p) + mass_above_root(no) *1.e-03      ! kg

    end do

    ! Averaging for each species

    do p = 1, n_spe
       mean_tree_h(p) = mean_tree_h(p) / max(1.0d0, tree_density_plot(p))      ! m
       mean_dbh(p) = mean_dbh(p) / max(1.0d0, tree_density_plot(p))            ! m
       above_biomass_plot(p) = (above_biomass_plot(p) / (real(Max_loc)**2.0d0)) * 10.0d0               ! kg/m2 -> Mg/ha
       leaf_biomass_plot(p) = (leaf_biomass_plot(p) / (real(Max_loc)**2.0d0)) * 10.0d0                 ! kg/m2 -> Mg/ha
       root_biomass_plot(p) = (root_biomass_plot(p) / (real(Max_loc)**2.0d0)) * 10.0d0                 ! kg/m2 -> Mg/ha
       coarse_root_biomass_plot(p) = (coarse_root_biomass_plot(p) / (real(Max_loc)**2.0d0)) * 10.0d0   ! kg/m2 -> Mg/ha
       above_root_biomass_plot(p) = (above_root_biomass_plot(p) / (real(Max_loc)**2.0d0)) * 10.0d0   ! kg/m2 -> Mg/ha
    end do

    ! Compute tree density

    do p = 1, n_spe
       tree_density_plot(p) = tree_density_plot(p) / (real(Max_loc)**2.0d0)   ! tree/m2
    end do

    ! Write tree biomass output

    write (Fn1, '(2(i4,a))',advance='no') year,',', day_of_year, ','
    do p = 1, n_spe
       write (Fn1, '(f12.5,a)',advance='no') tree_density_plot(p),','
    end do
    do p = 1, n_spe
       write (Fn1, '(f12.5,a)',advance='no') mean_tree_h(p),','
    end do
    do p = 1, n_spe
       write (Fn1, '(f12.5,a)',advance='no') mean_dbh(p),','
    end do
    do p = 1, n_spe
       write (Fn1, '(f12.5,a)',advance='no') above_biomass_plot(p),','
    end do
    do p = 1, n_spe
       write (Fn1, '(f12.5,a)',advance='no') leaf_biomass_plot(p),','
    end do
    do p = 1, n_spe
       write (Fn1, '(f12.5,a)',advance='no') coarse_root_biomass_plot(p),','
    end do
    do p = 1, n_spe
       write (Fn1, '(f12.5,a)',advance='no') above_root_biomass_plot(p),','
    end do
    do p = 1, n_spe
       if (p < n_spe) then
       write (Fn1, '(f12.5,a)',advance='no') root_biomass_plot(p),','
       else
       write (Fn1, '(f12.5,a)') root_biomass_plot(p)
       end if
    end do

write(*,*) 'Rh      year', '      day', '     tree number', '     height_mean', '     dbh_mean'
write(*,*) year, day_of_year, tree_density_plot(1) * 7.0d0 * 7.0d0 * 3.1415d0, mean_tree_h(1), mean_dbh(1)
write(*,*) 'Br      year', '      day', '     tree number', '     height_mean', '     dbh_mean'
write(*,*) year, day_of_year, tree_density_plot(2) * 7.0d0 * 7.0d0 * 3.1415d0, mean_tree_h(2), mean_dbh(2)

    !---------------------------------------------------------------------
    ! Plot output
    !---------------------------------------------------------------------

    ! Zero out

    lai_plot = 0.0d0

    ! Summing for all layers

    do i = 1, Max_hgt
       lai_plot = lai_plot + dpai_layer_sum(i)
    end do

    ! Write plot output

    write (Fn2, '(2(i4,a))',advance='no') year,',', day_of_year, ','
    write (Fn2, '(5(f12.5,a))',advance='no') irveg,',', irsoi, ',', lai_plot, ',', ir_tree, ',', rnsoi_day, ','
    write (Fn2, '(3(f12.5,a))',advance='no') rnveg_day,',', etsoi_day, ',', etveg_day, ','
    write (Fn2, '(2(f12.5,a))',advance='no') sal * 58.44d0 / 1000.0d0,',', din, ','
    write (Fn2, '(f12.5,a)') dip

  end subroutine monitor_plot

end module mod_monitoring
