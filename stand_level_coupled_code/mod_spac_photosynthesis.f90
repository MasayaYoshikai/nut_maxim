module mod_spac_photosynthesis

  !-----------------------------------------------------------------------
  ! !Description:
  ! SPAC-Photosynthesis model
  !
  ! !Uses:
  use mod_param
  use data_structure, only : Max_hgt, STEP, File_no, Fn_diurnal, growth_calc
  use time_counter, only : year
  use vegi_status_current1, only : monitor, bottom_layer_monitor, top_layer_monitor, &
                                   gs_monitor, an_top_monitor, an_bot_monitor, &
                                   et_top_monitor, et_bot_monitor
  !
  ! !Public types:
  implicit none
  !
  ! !Public member functions:
  public :: spac_photosynthesis
  !
  ! !Private member functions:
  !
  !
  !-----------------------------------------------------------------------

  contains

  !-----------------------------------------------------------------------
  subroutine spac_photosynthesis ( &
                                       ! *** Input ***
    o2air                         , &  ! Atmospheric O2 (mmol/mol)
    co2air                        , &  ! Atmospheric CO2 (umol/mol)
    hksat                         , &  ! Soil hydraulic conductivity at saturation (mm H2O/s)
    moist                         , &  ! Relative soil moisture (-)
    bsw                           , &  ! Soil layer Clapp and Hornberger "b" parameter (-)
    psisat                        , &  ! Soil matric potential at saturation (MPa)
    soil_t                        , &  ! Soil water temperature (K)
    sal                           , &  ! Pore-water salinity (mol/m3)
    din                           , &  ! DIN concentration in pore-water (mol N/m3)
    root_filter                   , &  ! Filtering rate of salt at root surface (-)
    root_resist                   , &  ! Hydraulic resistivity of root tissue (MPa.s.g root/mmol H2O)
    root_density                  , &  ! Specific root density (fine root) (g root/m3 root)
    root_radius                   , &  ! Fine root radius (m)
    root_depth                    , &  ! Rooting depth (m)
    k_sap                         , &  ! Stem hydraulic conductivity (kg H2O.m/m2 sapwood/s/MPa)
    wood_rho                      , &  ! Wood density (g/cm3)
    c_n_stem                      , &  ! C/N ratio in mol in stem
    c_n_root                      , &  ! C/N ratio in mol in root
    t_growth                      , &  ! Growth temperature (degree), mean temperature for 30 days before
    t_home                        , &  ! Home temperature (degree), mean maximum temperature of the warmest month
    minlp                         , &  ! Minimum leaf water potential (MPa)
    dleaf                         , &  ! Leaf dimension (m)
    gs_max                        , &  ! Maximum stomatal conductance (mol H2O/m2/s)
    vcmaxpft                      , &  ! Maximum carboxylation rate at 25C (umol/m2/s)
    iota                          , &  ! Stomatal water-use efficiency (umol CO2/ mol H2O)
    C_in_drymass                  , &  ! Proportion of carbon in biomass (g C/g DW)
    ALM5                          , &  ! Sapwood diameter proportion (m sapwood/m dbh)
                                       !
                                       ! *** From SEIB-DGVM ***
    day_of_year                   , &  ! Day of the year (1 - 365)
    no                            , &  ! Tree index
    p                             , &  ! Species index (1: Rh, 2: Br)
    rPAR_dir                      , &  ! Profile of relative intensity of direct PAR within canopy compared to canopy top (fraction)
    rPAR_dif                      , &  ! Profile of relative intensity of diffused PAR within canopy compared to canopy top (fraction)
    rNIR_dir                      , &  ! Profile of relative intensity of direct NIR within canopy compared to canopy top (fraction)
    rNIR_dif                      , &  ! Profile of relative intensity of diffused NIR within canopy compared to canopy top (fraction)
    irleaf                        , &  ! Leaf absorbed longwave radiation for canopy layer (W/m2 leaf)
    kbm_vis                       , &  ! Scattering adjusted light extinction coefficient for VIS (-)
    kbm_nir                       , &  ! Scattering adjusted light extinction coefficient for NIR (-)
    rwind_profile                 , &  ! Relative wind speed profile to the canopy top (-)
    tree_h                        , &  ! Tree height (m)
    root_biomass                  , &  ! Stand-level fine root biomass (g root/tree)
    lai                           , &  ! Leaf area index (m2 leaf/m2 ground)
    crown_area                    , &  ! Stand-level leaf crown area (m2 ground/tree)
    dbh_heartwood                 , &  ! Heartwood diameter (m)
    dbh_sapwood                   , &  ! Sapwood diameter (m)
    trunk_biomass                 , &  ! Stand-level trunk biomass (g/tree)
    seib_height                   , &  ! SEIB-DGVM specific tree height (the unit is STEP!!)
    seib_bole                     , &  ! SEIB-DGVM specific bole height (the unit is STEP!!)
    time_series_tair              , &  ! Time-series air temperature (K)
    time_series_pair              , &  ! Time-series air pressure (Pa)
    time_series_wind              , &  ! Time-series wind speed (m/s)
    time_series_eair              , &  ! Time-series vapor pressure in air (Pa)
    time_series_rad_dir           , &  ! Time-series direct radiation at canopy top (W/m2)
    time_series_rad_dif           , &  ! Time-series diffused radiation at canopy top (W/m2)
                                       !
                                       ! *** Output ***
    n_leaf_layer                  , &  ! Number of leaf layers
    max_water_uptake_day          , &  ! Maximum stand-level water uptake rate (mol H2O/tree/s)
    max_stand_et_day              , &  ! Maximum stand-level transpiration rate of the day (mol H2O/tree/s)
    sapflow_day                   , &  ! Daily sapflow rate at stand-level (m3 H2O/tree/day)
    an_mean_day_max               , &  ! Daily maximum canopy mean leaf net photosynthesis (umol CO2/m2 leaf/s)
    an_top_day_max                , &  ! Daily maximum canopy top leaf net photosynthesis (umol CO2/m2 leaf/s)
    c_uptake_day                  , &  ! Daily carbon uptake rate by photosynthesis (mol C/tree/day)
    c_uptake_bottom_day           , &  ! Daily bottom layer carbon uptake rate by photosynthesis (mol CO2/m2 leaf/day)
    n_uptake_bottom_day           , &  ! Daily bottom layer nitrogen uptake rate (mol N/m2 leaf/day)
    n_uptake_day                  , &  ! Daily nitrogen uptake rate (mol N/tree/day)
    r_whole_root_increment        , &  ! Whole plant resistance after root biomass increment (MPa.s.tree/mmol H2O)
    r_whole_d_increment           , &  ! Whole plant resistance after diameter increment (MPa.s.tree/mmol H2O)
    r_root_out                    , &  ! Root-stem resistance (MPa.s.tree/mmol H2O)
    r_sap_out                     , &  ! Whole-plant stem resistance (MPa.s.tree/mmol H2O)
    tleaf_out                       &  ! Leaf temperature profile (K)
    )
    !
    ! !Description:
    ! Compute transpiration and photosynthesis for a tree
    !
    ! !Uses:
    use mod_water_vapor
    use mod_math_tools
    use mod_soil_char
    use mod_plant_hydraulics
    use mod_nitrogen_profile
    use mod_photosynthesis
    use mod_leaf_boundary_layer
    use mod_energy_water_balance
    use mod_stomatal_conductance
    use mod_tree_allometry
    !
    ! !Argumetns:
    implicit none
    real(8), intent(in)  :: o2air, co2air
    real(8), intent(in)  :: hksat, moist, bsw, psisat, soil_t
    real(8), intent(in)  :: sal, din
    real(8), intent(in)  :: root_filter(:), root_resist(:)
    real(8), intent(in)  :: root_density(:), root_radius(:), root_depth
    real(8), intent(in)  :: k_sap(:), wood_rho(:), c_n_stem(:), c_n_root(:)
    real(8), intent(in)  :: t_growth(:), t_home(:)
    real(8), intent(in)  :: minlp(:), dleaf(:), gs_max(:), vcmaxpft(:), iota(:)
    real(8), intent(in)  :: C_in_drymass, ALM5(:)
    integer, intent(in)  :: day_of_year, no, p
    real(8), intent(in)  :: rPAR_dir(:), rPAR_dif(:), rNIR_dir(:), rNIR_dif(:), irleaf(:)
    real(8), intent(in)  :: kbm_vis(:), kbm_nir(:), rwind_profile(:)
    real(8), intent(in)  :: tree_h, root_biomass, lai, crown_area
    real(8), intent(in)  :: dbh_heartwood, dbh_sapwood, trunk_biomass
    integer, intent(in)  :: seib_height, seib_bole
    real(8), intent(in)  :: time_series_tair(:), time_series_pair(:)
    real(8), intent(in)  :: time_series_wind(:), time_series_eair(:)
    real(8), intent(in)  :: time_series_rad_dir(:), time_series_rad_dif(:)
    integer, intent(out) :: n_leaf_layer
    real(8), intent(out) :: max_water_uptake_day, max_stand_et_day, sapflow_day
    real(8), intent(out) :: an_mean_day_max, an_top_day_max
    real(8), intent(out) :: c_uptake_day, c_uptake_bottom_day, n_uptake_bottom_day
    real(8), intent(out) :: n_uptake_day
    real(8), intent(out) :: r_whole_root_increment, r_whole_d_increment
    real(8), intent(out) :: r_root_out, r_sap_out
    real(8), intent(out) :: tleaf_out(:)
    !
    ! !Local variables:
    type(universal_type) :: universal
    type(atmos_type)     :: atmos
    type(layer_type)     :: layer
    type(flux_type)      :: flux
    integer :: i
    integer :: count
    integer :: hour  ! 0 - 23
    integer :: canopy_index
    integer :: bottom_layer                    ! Bottom layer of canopy
    integer :: hour_of_year
    real(8) :: root_area                       ! Stand-level root coverage area (m2 ground/tree)
    real(8) :: bole_h                          ! Bole height (m)
    real(8) :: tot_psi                         ! Total soil water potential (MPa)
    real(8) :: hk                              ! Soil hydraulic conductance (mmol H2O/m/s/MPa)
    integer :: flag_energy                     ! Flag (1: for stomatal optimization, 0: other purposes)
    real(8) :: val                             ! Dummy for calling the function energy_water_balance
    real(8) :: max_stand_et                    ! Maximum total transpiration rate at stand-level (mol H2O/tree/s)
    real(8) :: stand_et                        ! Stand-level transpiration rate calculated from the determined gs (mol H2O/tree/s)
    real(8) :: swsky_vis_dir                   ! Atmospheric direct VIS solar radiation (W/m2)
    real(8) :: swsky_vis_dif                   ! Atmospheric diffused VIS solar radiation (W/m2)
    real(8) :: swsky_nir_dir                   ! Atmospheric direct NIR solar radiation (W/m2)
    real(8) :: swsky_nir_dif                   ! Atmospheric diffused NIR solar radiation (W/m2)
    real(8) :: par_dir                         ! Direct PAR at canopy top (umol photon/m2 ground/s)
    real(8) :: par_dif                         ! Diffused PAR at canopy top (umol photon/m2 ground/s)
    real(8) :: pre_rad                         ! Net radiation of one-hour ago (W/m2)
    real(8) :: now_rad                         ! Current net radiation (W/m2)
    real(8) :: swleaf                          ! Leaf absorbed shortwave radiation per ground (W/m2 ground)
    real(8), dimension(Max_hgt) :: save_an     ! Save value of an(:)
    real(8), dimension(Max_hgt) :: save_et     ! Save value of etflux(:)
    real(8) :: qair                            ! Specific humidity (kg/kg)
    real(8) :: rhoair                          ! Air density (kg/m3)
    real(8) :: mmair                           ! Molecular mass of air (kg/mol)
    real(8) :: sapflow                         ! Stand-level sapflow rate  (mol H2O/tree/s)
    real(8) :: sapflow_m3                      ! Stand-level sapflow rate  (m3 H2O/tree/s)
    real(8) :: an_mean                         ! Canopy mean leaf net photosynthesis (umol CO2/m2 leaf/s)
    real(8) :: an_top                          ! Canopy top leaf net photosynthesis (umol CO2/m2 leaf/s)
    real(8) :: c_uptake                        ! Carbon uptake rate by photosynthesis (mol C/tree/s)
    real(8) :: c_uptake_bottom                 ! Bottom layer carbon uptake rate by photosynthesis (mol CO2/m2 leaf/s)
    real(8) :: et_mean                         ! Canopy mean leaf transpiration rate (mol H2O/m2 leaf/s)
    real(8) :: c_eff                           ! Carbon uptake efficiency (Negative: photosynthesis limitation, Positive: nutrient uptake limitation)
    real(8) :: et_bottom                       ! Bottom layer transpiration rate (m3 H2O/m2 leaf/s)
    real(8) :: below_resource_day_stem         ! Stem biomass increment potential based on nutrient uptake (g DW/tree/day)
    real(8) :: below_resource_day_root         ! Root biomass increment potential based on nutrient uptake (g DW/tree/day)
    real(8) :: rd_stand                        ! Leaf respiration rate at stand-level (mol C/tree/s)
    real(8) :: biomass_increment               ! Increment of trunk biomass for calculating hydraulic bottleneck (g/tree)
    real(8) :: new_dbh_heartwood               ! New heartwood diameter after biomass increment (m)
    real(8) :: new_dbh_sapwood                 ! New sapwood diameter after biomass increment (m)
    real(8) :: new_tree_h                      ! New tree height after biomass increment (m)
    real(8) :: new_root_biomass                ! New fine root biomass after biomass increment (g root/tree)
    !---------------------------------------------------------------------

    ! Allocate model parameters

    call init_allocate ( &
                              ! *** In/Output ***
    universal            , &  ! Universal variables
    atmos                , &  ! Atmospheric variables
    layer                , &  ! Layer variables
    flux                   &  ! Flux variables
    )

    associate ( &
                                                     ! *** Input ***
    n_layer      => universal%n_layer           , &  ! Number of layers (-) (Max_hgt in SEIB-DGVM)
    rgas         => universal%rgas              , &  ! Universal gas constant (J/K/mol)
    denh2o       => universal%denh2o            , &  ! Water density (kg/m3)
    mmh2o        => universal%mmh2o             , &  ! Molecular mass of water (kg/mol)
    tfrz         => universal%tfrz              , &  ! Freezing point of water (K)
    mmdry        => universal%mmdry             , &  ! Molecular mass of dry air (kg/mol)
    cpd          => universal%cpd               , &  ! Specific heat of dry air at constant pressure (J/kg/K)
    cpw          => universal%cpw               , &  ! Specific heat of water vapor at constant pressure (J/kg/K)
                                                     !
                                                     ! *** Derived variables
    tair         => atmos%tair                  , &  ! Air temperature (K)
    pair         => atmos%pair                  , &  ! Air pressure (Pa)
    wind         => atmos%wind                  , &  ! Wind speed (m/s)
    eair         => atmos%eair                  , &  ! Vapor pressure in air (Pa)
    rhomol       => atmos%rhomol                , &  ! Molar density (mol/m3)
    cp           => atmos%cp                    , &  ! Specific heat of air at constant pressure (J/mol/K)
    rad_dir      => atmos%rad_dir               , &  ! Direct radiation at canopy top (W/m2)
    rad_dif      => atmos%rad_dif               , &  ! Diffused radiation at canopy top (W/m2)
    leaf_layer   => layer%leaf_layer            , &  ! 1: Leaf, 0: No leaf
    sumpai       => layer%sumpai                , &  ! Cumulative leaf area (m2 leaf/m2 ground)
    dpai         => layer%dpai                  , &  ! Layer leaf area index (m2 leaf/m2 ground)
    rn           => flux%rn                     , &  ! Leaf net radiation profile (W/m2 leaf)
    apar         => flux%apar                   , &  ! Leaf absorbed PAR Profile (umol photon/m2 leaf/s)
    etflux       => flux%etflux                 , &  ! Leaf transpiration rate (mol H2O/m2 leaf/s)
    tleaf        => flux%tleaf                  , &  ! Leaf temperature (K)
    an           => flux%an                     , &  ! Leaf net photosynthesis (umol CO2/m2 leaf/s)
    gs           => flux%gs                     , &  ! Stand-level stomatal conductance (mol H2O/m2/s)
    rd           => flux%rd                     , &  ! Leaf respiration rate (umol CO2/m2 leaf/s)
    sapmax       => flux%sapmax                 , &  ! Stand-level maximum water uptake rate (mol H2O/tree/s)
    r_root       => flux%r_root                 , &  ! Root-stem resistance (MPa.s.tree/mmol H2O)
    r_sap        => flux%r_sap                  , &  ! Whole-plant stem resistance (MPa.s.tree/mmol H2O)
    r_whole      => flux%r_whole                  &  ! Whole plant resistance (MPa.s.tree/mmol H2O)
    )

    ! Stand-level root coverage area (m2 ground/tree): assume root area is same as crown area.

    root_area = crown_area

    ! Bole height (m)

    bole_h = real(seib_bole) * STEP + 1.3d0

    ! Leaf layer index (1: Leaf, 0: No leaf)                                              ! <= to be returned.

    n_leaf_layer = 0
    do i = 1, n_layer
       if (i >= (seib_bole + 1) .and. i <= seib_height) then
          leaf_layer(i) = 1
          n_leaf_layer = n_leaf_layer + 1
       else
          leaf_layer(i) = 0
       end if
    end do

    ! Identifying canopy bottom layer

    canopy_index = 0
    do i = 1, n_layer
       if (leaf_layer(i) == 1 .and. canopy_index == 0) then
          bottom_layer = i
          canopy_index = 1
       end if
    end do

    ! Layer leaf area index (m2 leaf/m2 ground)

    dpai = lai / real(n_leaf_layer)

    ! Cumulative leaf area from canopy top (m2 leaf/m2 ground)

    count = 0
    do i = n_layer, 1, -1
       if (leaf_layer(i) == 1) then ! Leaf layer
          count = count + 1
          if (count == 1) then ! Canopy top
             sumpai(i) = dpai
          else
             sumpai(i) = sumpai(i+1) + dpai
          end if
       else ! non-leaf layer
          sumpai(i) = 0.0d0
       end if
    end do

    ! Set soil hydraulic parameters

    call soil_char (       &
                              ! *** Input ***
    hksat                , &  ! Soil hydraulic conductivity at saturation (mm H2O/s)
    moist                , &  ! Relative soil moisture (-)
    bsw                  , &  ! Soil layer Clapp and Hornberger "b" parameter (-)
    psisat               , &  ! Soil matric potential at saturation (MPa)
    soil_t               , &  ! Soil water temperature (K)
    sal                  , &  ! Pore-water salinity (mol/m3)
    root_filter(p)       , &  ! Filtering rate of salt at root surface (-)
                              !
                              ! *** Input type ***
    universal            , &  ! Universal variables
                              !
                              ! *** Output ***
    tot_psi              , &  ! Total soil water potential (MPa)
    hk                     &  ! Soil hydraulic conductance (mmol H2O/m/s/MPa)
    )

    ! Calculate plant hydraulic resistance and maximum water uptake rate

    call plant_hydraulics (&
                              ! *** Input ***
    root_radius(p)       , &  ! Fine root radius (m)
    root_biomass         , &  ! Stand-level fine root biomass (g root/tree)
    root_area            , &  ! Stand-level root coverage area (m2 ground/tree)
    root_depth           , &  ! Rooting depth (m)
    root_density(p)      , &  ! Specific root density (fine root) (g root/m3 root)
    hk                   , &  ! Soil hydraulic conductance (mmol H2O/m/s/MPa)
    root_resist(p)       , &  ! Hydraulic resistivity of root tissue (MPa.s.g root/mmol H2O)
    dbh_heartwood        , &  ! Heartwood diameter (m)
    dbh_sapwood          , &  ! Sapwood diameter (m)
    k_sap(p)             , &  ! Stem hydraulic conductivity (kg H2O.m/m2 sapwood/s/MPa)
    tree_h               , &  ! Tree height (m)
    bole_h               , &  ! Bole height (m)
    tot_psi              , &  ! Total soil water potential (MPa)
    minlp(p)             , &  ! Minimum leaf water potential (MPa)
                              !
                              ! *** Input type ***
    universal            , &  ! Universal variables
                              !
                              ! *** Input/Output ***
    flux                   &  ! Flux variables
    )
    r_root_out = r_root                                                                   ! <= to be returned.
    r_sap_out = r_sap                                                                     ! <= to be returned.

    ! Maximum stand-level water uptake rate (mol H2O/tree/s)

    max_water_uptake_day = sapmax                                                         ! <= to be returned.

!write(*,*) 'sapmax kg H2O/tree/s', sapmax*mmh2o

    ! Canopy profile of nitrogen and photosynthetic capacity

    call nitrogen_profile (&
                              ! *** Input ***
    t_growth(p)          , &  ! Growth temperature (degree), mean temperature for 30 days before
    t_home(p)            , &  ! Home temperature (degree), mean maximum temperature of the warmest month
    vcmaxpft(p)          , &  ! Maximum carboxylation rate at 25C (umol/m2/s)
                              !
                              ! *** Input type ***
    universal            , &  ! Universal variables
                              !
                              ! *** Input/Output ***
    layer                  &  ! Layer variables
    )

    ! Set leaf-level parameters for photosynthesis model

    call photosynthesis_param (&
                              ! *** Input ***
    o2air                , &  ! Atmospheric O2 (mmol/mol)
    t_growth(p)          , &  ! Growth temperature (degree), mean temperature for 30 days before
    t_home(p)            , &  ! Home temperature (degree), mean maximum temperature of the warmest month
                              !
                              ! *** Input type ***
    universal              &  ! Universal variables
    )

    ! Initialize daily stand-level variables to return

    max_stand_et_day = 0.0d0
    sapflow_day = 0.0d0
    an_mean_day_max = 0.0d0
    an_top_day_max = 0.0d0
    c_uptake_day = 0.0d0
    c_uptake_bottom_day = 0.0d0
    n_uptake_bottom_day = 0.0d0

    !---------------------------------------------------------------------
    ! Hourly computation
    !---------------------------------------------------------------------

    do hour = 0, 23

       ! Set atmospheric forcing parameters

       hour_of_year = (day_of_year -1)*24 + hour + 1
       tair = time_series_tair(hour_of_year)          ! Air temperature (K)
       pair = time_series_pair(hour_of_year)          ! Air pressure (Pa)
       wind = time_series_wind(hour_of_year)          ! Wind speed (m/s)
       eair = time_series_eair(hour_of_year)          ! Vapor pressure in air (Pa)
       rad_dir = time_series_rad_dir(hour_of_year)    ! Direct radiation at canopy top (W/m2)
       rad_dif = time_series_rad_dif(hour_of_year)    ! Diffused radiation at canopy top (W/m2)

       ! Atmospheric VIS and NIR solar radiation (W/m2)
       ! The proportion is same as that of SEIB-DGVM. Assuming 4% of solar radiation is UV.

       swsky_vis_dir = rad_dir * 0.43d0
       swsky_vis_dif = rad_dif * 0.57d0
       swsky_nir_dir = rad_dir * (1.0d0 - 0.43d0 - 0.04d0)
       swsky_nir_dif = rad_dif * (1.0d0 - 0.57d0 - 0.04d0)

       ! Direct and diffused PAR at canopy top (umol photon/m2 ground/s)
       ! The proportion is same as that of SEIB-DGVM.

       par_dir = 4.6d0 * swsky_vis_dir
       par_dif = 4.2d0 * swsky_vis_dif

       ! Radiation of one-hour ago and now (W/m2)

       if (hour > 0) then
          pre_rad =  time_series_rad_dir(hour_of_year-1) + time_series_rad_dif(hour_of_year-1)
          now_rad =  rad_dir + rad_dif
       end if

       !---------------------------------------------------------------------
       ! Set atmospheric parameters
       !---------------------------------------------------------------------

!       qair = 0.622d0 * eair / (pair - 0.378d0 * eair)
       qair = mmh2o / mmdry * eair / (pair - (1.0d0 - mmh2o/mmdry) * eair)       ! Specific humidity (kg/kg)
       rhomol = pair / (rgas * tair)                                             ! Molar density (mol/m3)
!       rhoair = rhomol * mmdry * (1.0d0 - 0.378d0 * eair/pair);
       rhoair = rhomol * mmdry * (1.0d0 - (1.0d0 - mmh2o/mmdry) * eair / pair)   ! Air density (kg/m3)
       mmair = rhoair / rhomol                                                   ! Molecular mass of air (kg/mol)
!       cp = cpd * (1.0d0 + 0.84d0 * qair) * mmair
       cp = cpd * (1.0d0 + (cpw/cpd - 1.0d0) * qair) * mmair;                    ! Specific heat of air at constant pressure (J/mol/K)

       !---------------------------------------------------------------------
       ! Leaf net radiation (W/m2 leaf) and absorbed PAR (umol photon/m2 leaf/s)
       !
       !          |                          |
       !          v       par_layer(i)       v
       !        ________________________________
       !       |
       !       |        Canopy layer (i) <= apar(i)
       !       |________________________________
       !          |                          |
       !          v      par_layer(i-1)      v
       !
       !
       !          |                          |
       !          v       par_layer(k)       v
       !        ________________________________
       !       |
       !       |     Bottom canopy layer (k) <= apar(k)
       !       |________________________________
       !          |                          |
       !          v    par_layer(k-1) = 0.0  v
       !        ________________________________
       !       |
       !       |       Non-leaf layer (k-1)
       !       |________________________________
       !          |                          |
       !          v    par_layer(k-1) = 0.0  v
       !
       !---------------------------------------------------------------------

       ! Absorbed radiation and PAR at each layer

       do i = n_layer, 1, -1
          if (leaf_layer(i) == 1) then ! Leaf layer

!!! ---------- Based on differences between upper and lower layers --------- !!!
!             ! Absorbed shortwave radiation and PAR per ground area
!
!             swleaf = (swsky_vis_dir * rPAR_dir(i) + swsky_vis_dif * rPAR_dif(i)) &
!                      * (1.0d0 - exp(-1.0d0 * kbm_vis(p) * dpai)) &
!                      + (swsky_nir_dir * rNIR_dir(i) + swsky_nir_dif * rNIR_dif(i)) &
!                      * (1.0d0 - exp(-1.0d0 * kbm_nir(p) * dpai))
!
!             apar(i) = (par_dir * rPAR_dir(i) + par_dif * rPAR_dif(i)) &
!                       * (1.0d0 - exp(-1.0d0 * kbm_vis(p) * dpai))
!
!             ! per leaf area
!
!             rn(i) = swleaf / dpai + irleaf(i)
!             rn(i) = max(rn(i), 0.0d0)          ! To prevent negative transpiration rate.
!             if (abs(irleaf(i)) < 1.e-06) then
!                write(*,*) 'ERROR: Zero longwave radiation in SPAC-Photosynthesis module', irleaf(i)
!             end if
!             apar(i) = apar(i) / dpai
!!! ----------------- Based on light absorption coefficient ---------------- !!!
             swleaf = (swsky_vis_dir * rPAR_dir(i) + swsky_vis_dif * rPAR_dif(i)) * kbm_vis(p) &
                      + (swsky_nir_dir * rNIR_dir(i) + swsky_nir_dif * rNIR_dif(i)) * kbm_nir(p)
             rn(i) = swleaf + irleaf(i)
             apar(i) = (par_dir * rPAR_dir(i) + par_dif * rPAR_dif(i)) * kbm_vis(p)
!!! ----------------------------------------------------------------------- !!!

          else ! non-leaf layer

             rn(i) = 0.0d0
             apar(i) = 0.0d0
          end if
       end do

       !---------------------------------------------------------------------
       ! Leaf boundary layer conductances
       !---------------------------------------------------------------------

       call leaf_boundary_layer (      &
                                          ! *** Input ***
       dleaf(p)                      , &  ! Leaf dimension (m)
       rwind_profile                 , &  ! Relative wind speed profile to the canopy top (-)
                                          !
                                          ! *** Input type ***
       universal                     , &  ! Universal variables
       atmos                         , &  ! Atmospheric variables
       layer                         , &  ! Layer variables
                                          !
                                          ! *** Input/Output ***
       flux                            &  ! Flux variables
       )

       !---------------------------------------------------------------------
       ! Maximum stand-level transpiration rate (mol H2O/tree/s)
       ! Calculating leaf temperature and transpiration profile with maximum gs.
       ! val: (dummy) sapmax - stand_et
       !---------------------------------------------------------------------

       flag_energy = 0
       val = energy_water_balance (    &
                                          ! *** Input ***
       flag_energy                   , &  ! Flag (1: for stomatal optimization, 0: other purposes)
       crown_area                    , &  ! Stand-level leaf crown area (m2 ground/tree)
       gs_max(p)                     , &  ! Value for gs to use in calculations
                                          !
                                          ! *** Input type ***
       universal                     , &  ! Universal variables
       atmos                         , &  ! Atmospheric variables
       layer                         , &  ! Layer variables
                                          !
                                          ! *** Input/Output ***
       flux                            &  ! Flux variables
       )

       max_stand_et = 0.0d0
       do i = 1, n_layer
          if (leaf_layer(i) == 1) then
             max_stand_et = max_stand_et + etflux(i) * dpai * crown_area
          end if
       end do

!write(*,*) 'Maximum transpiration done'
!write(*,*) 'max_stand_et kg H2O/tree/s', max_stand_et*mmh2o

       !---------------------------------------------------------------------
       ! Stomatal conductance, leaf transpiration, and photosynthesis.
       ! Skip computation during succesive night time (one-hour previous value of an is used)
       ! Stand-level leaf respiration rate is computed when hour = 0 or radiation is zero.
       ! <- Growth rate is negative during night time.
       !---------------------------------------------------------------------

       if (hour > 0 .and. pre_rad < 1.0d0 .and. now_rad < 1.0d0) then ! succesive night time

          an = save_an
          etflux = save_et

       else ! 0:00 AM or day-time

          !---------------------------------------------------------------------
          ! Stomatal conductance based on water uptake maximization
          !---------------------------------------------------------------------

          if (sapmax > max_stand_et) then

             ! Transpiration (leaf area) is constraining sapflow.
             ! Stomatal conductance is set to the maximum value.

             gs = gs_max(p)

          else

             ! Plant hydraulics is constraining sapflow.
             ! Plant needs to regulate gs not to lose water by transpiration.
             ! Optimizing stomatal conductance to maximize water uptake rate
             ! then compute leaf temperature and transpiration profile using
             ! the optimized gs.

             call stomatal_conductance1 (    &
                                               ! *** Input ***
             crown_area                   , &  ! Stand-level leaf crown area (m2 ground/tree)
             gs_max(p)                    , &  ! Maximum stomatal conductance (mol H2O/m2/s)
                                               !
                                               ! *** Input type ***
             universal                    , &  ! Universal variables
             atmos                        , &  ! Atmospheric variables
             layer                        , &  ! Layer variables
                                               !
                                               ! *** Input/Output ***
             flux                           &  ! Flux variables
             )

          end if

          ! Leaf photosynthesis for a specified stomatal conductance
          ! based on the water uptake maximization

          call leaf_photosynthesis (        &
                                               ! *** Input ***
          o2air                           , &  ! Atmospheric O2 (mmol/mol)
          co2air                          , &  ! Atmospheric CO2 (umol/mol)
                                               !
                                               ! *** Input type ***
          universal                       , &  ! Universal variables
          atmos                           , &  ! Atmospheric variables
          layer                           , &  ! Layer variables
                                               !
                                               ! *** Input/Output ***
          flux                              &  ! Flux variables
          )

          ! Check which factor photosynthesis or water uptake could growth rate at the time.
          ! Stem C/N is used for the stoichiometry here.
          !
          ! an_mean (temporal use): Canopy mean leaf net photosynthesis (umol CO2/m2 leaf/s)
          ! et_mean (temporal use): Canopy mean leaf transpiration rate (mol H2O/m2 leaf/s)
          ! c_eff: Carbon uptake efficiency
          !        (Negative: photosynthesis limitation, Positive: nutrient uptake limitation)

          an_mean = 0.0d0
          et_mean = 0.0d0
          do i = 1, n_layer
             if (leaf_layer(i) == 1) then ! Leaf layer
                an_mean = an_mean + an(i)
                et_mean = et_mean + etflux(i)
             end if
          end do
          an_mean = an_mean / real(n_leaf_layer)
          et_mean = et_mean / real(n_leaf_layer)

          c_eff = an_mean * 1.e-06 &                                ! mol C/m2 leaf/s
                  - (et_mean * mmh2o / denh2o) * din * c_n_root(p)  ! mol C/m2 leaf/s

          ! Case low photosynthetic rate (like due to dark condition) could limit the growth rate

          if (c_eff < 0.0d0) then

! write(*,*) 'Photosynthesis limitation: hour = ', hour

             !---------------------------------------------------------------------
             ! Stomatal conductance based on water use efficiency
             ! Leaf transpiration and photosynthesis with the new gs
             ! are calculated in this routine.
             !---------------------------------------------------------------------

             call stomatal_conductance2 (    &
                                               ! *** Input ***
             o2air                        , &  ! Atmospheric O2 (mmol/mol)
             co2air                       , &  ! Atmospheric CO2 (umol/mol)
             iota(p)                      , &  ! Stomatal water-use efficiency (umol CO2/ mol H2O)
             crown_area                   , &  ! Stand-level leaf crown area (m2 ground/tree)
             n_leaf_layer                 , &  ! Number of leaf layers
                                               !
                                               ! *** Input type ***
             universal                    , &  ! Universal variables
             atmos                        , &  ! Atmospheric variables
             layer                        , &  ! Layer variables
                                               !
                                               ! *** Input/Output ***
             flux                           &  ! Flux variables
             )

          end if

          ! Save assimilation (umol CO2/m2 leaf/s) and transpiration (mol H2O/m2 leaf/s)
          ! profile for successive night time

          save_an = an
          save_et = etflux

       end if ! End of stomatal conductance, leaf transpiration, and photosynthesis calculations.

!write(*,*) 'hour, gs, an(top), et(top)', hour, gs, an(bottom_layer+n_leaf_layer-1),etflux(bottom_layer+n_leaf_layer-1)

       !---------------------------------------------------------------------
       ! Stand-level sapflow rate (= transpiration rate)
       ! Plant hydraulics or maximum transpiration rate is limiting sapflow rate.
       !---------------------------------------------------------------------

       ! Stand-level transpiration rate calculated from the detemined gs.

       stand_et = 0.0d0
       do i = 1, n_layer
          if (leaf_layer(i) == 1) then
             stand_et = stand_et + etflux(i) * dpai * crown_area
          end if
       end do

       ! Stand-level sapflow rate

       sapflow = stand_et                     ! (mol H2O/tree/s)
       sapflow_m3 = sapflow * mmh2o / denh2o  ! (m3 H2O/tree/s)

       ! Maximum stand-level transpiration rate of the day (mol H2O/tree/s)

       max_stand_et_day = max(max_stand_et_day, stand_et)                                 ! <= to be returned.

       ! Bottom layer transpiration rate (m3 H2O/m2 leaf/s)
       ! Multiplying mmh2o/denh2o for (mol H2O/m2 leaf/s) -> (m3 H2O/m2 leaf/s)

       et_bottom = etflux(bottom_layer) * mmh2o / denh2o

       !---------------------------------------------------------------------
       ! Carbon uptake rate by photosynthesis and leaf respiration rate
       !---------------------------------------------------------------------

       ! Carbon uptake rate by photosynthesis (mol C/tree/s)
       ! Leaf respiration rate (mol C/tree/s)
       ! Multiplying 1.e-06 for unit conversion (umol CO2/tree/s) -> (mol CO2/tree/s)
       !
       ! an_mean: Canopy mean leaf net photosynthesis (umol CO2/m2 leaf/s)
       ! an_top: Canopy top leaf net photosynthesis (umol CO2/m2 leaf/s)

       c_uptake = 0.0d0
       rd_stand = 0.0d0
       an_mean = 0.0d0
       do i = 1, n_layer
          if (leaf_layer(i) == 1) then ! Leaf layer
             c_uptake = c_uptake + an(i) * dpai * crown_area
             rd_stand = rd_stand + rd(i) * dpai * crown_area
             an_mean = an_mean + an(i)
          end if
       end do
       c_uptake = c_uptake * 1.e-06
       rd_stand = rd_stand * 1.e-06
       an_mean = an_mean / real(n_leaf_layer)
       an_top = an(bottom_layer + n_leaf_layer - 1)
       if (an_top == 0.0d0) then
          write(*,*) 'ERROR: Zero an at canopy top.'
       end if

       ! Daily maximum canopy mean and canopy top leaf net photosynthesis (umol CO2/m2 leaf/s)

       an_mean_day_max = max(an_mean_day_max, an_mean)                                    ! <= to be returned.
       an_top_day_max = max(an_top_day_max, an_top)                                       ! <= to be returned.

       ! For bottom layer (mol CO2/m2 leaf/s)

       c_uptake_bottom = an(bottom_layer) * 1.e-06
       if (c_uptake_bottom == 0.0d0) then
          write(*,*) 'ERROR: zero bottom layer net assimilation rate'
          write(*,*) 'n_leaf_layer, an_top_day_max, tree_h, crown_area', n_leaf_layer, an_top_day_max, tree_h, crown_area
       end if

       !---------------------------------------------------------------------
       ! Daily sum for sap flow and carbon assimilation
       !---------------------------------------------------------------------

       ! Daily sapflow rate at stand-level (m3 H2O/tree/day)

       sapflow_day = sapflow_day + sapflow_m3 * 60.0d0 * 60.0d0                           ! <= to be returned.

       ! Daily carbon uptake rate by photosynthesis (mol C/tree/day)

       c_uptake_day = c_uptake_day + c_uptake * 60.0d0 * 60.0d0                           ! <= to be returned.

       ! Daily bottom layer carbon uptake rate by photosynthesis (mol CO2/m2 leaf/day)

       c_uptake_bottom_day = c_uptake_bottom_day + c_uptake_bottom * 60.0d0 * 60.0d0      ! <= to be returned.

       ! Daily bottom layer nitrogen uptake rate (mol N/m2 leaf/day)
       ! Multiplying din (mol N/m3) for (m3 H2O/m2 leaf/s) -> (mol N/m2 leaf/s)

       n_uptake_bottom_day = n_uptake_bottom_day + et_bottom * din * 60.0d0 * 60.0d0      ! <= to be returned.

       !---------------------------------------------------------------------
       ! For monitoring tree
       !---------------------------------------------------------------------

       if (no == monitor) then
          if (hour == 12) then ! at midday
             bottom_layer_monitor = bottom_layer
             top_layer_monitor = bottom_layer + n_leaf_layer - 1
             gs_monitor = gs
             an_top_monitor = an(top_layer_monitor)
             an_bot_monitor = an(bottom_layer_monitor)
             et_top_monitor = etflux(top_layer_monitor)
             et_bot_monitor = etflux(bottom_layer_monitor)
          end if
       end if

       !---------------------------------------------------------------------
       ! Diurnal dynamics (only for no-growth calculation)
       !---------------------------------------------------------------------

       if (growth_calc .eqv. .false.) then
          if (no == monitor .and. year == 3) then
             write(File_no(24), '(1(i4,a), 5(f12.5,a))') &
             hour_of_year                            , ',', &  ! Hour of year
             par_dir+par_dif                         , ',', &  ! PAR at canopy top (umol photon/m2 ground/s)
             tair-tfrz                               , ',', &  ! Air temperature (degree C)
             eair                                    , ',', &  ! Vapor pressure in air (Pa)
             gs                                      , ',', &  ! Stand-level stomatal conductance (mol H2O/m2/s)
             an(top_layer_monitor)                             ! Leaf net photosynthesis at canopy top (umol CO2/m2 leaf/s)
          end if
       end if

    end do ! End of hourly computation

    !---------------------------------------------------------------------
    ! Daily nitrogen uptake rate (mol N/tree/day)
    !---------------------------------------------------------------------

    n_uptake_day = sapflow_day * din                                                      ! <= to be returned.

    !---------------------------------------------------------------------
    ! Daily biomass increment (g DW/tree/day)
    ! c_uptake_day is not necessary here because
    ! when we need to consider plant hydraulics optimization,
    ! the limiting factor is always n_uptake_day.
    !---------------------------------------------------------------------

    ! Based on nutrient uptake.
    ! Multiplying c_n_stem/root for (mol N/tree/day) -> (mol C/tree/day)
    ! Multiplying 12 for (mol C/tree/day) -> (g C/tree/day)
    ! Deviding by C_in_drymass for (g C/tree/day) -> (g DW/tree/day)

    below_resource_day_stem = n_uptake_day * c_n_stem(p) * 12.0d0 / C_in_drymass
    below_resource_day_root = n_uptake_day * c_n_root(p) * 12.0d0 / C_in_drymass

    ! New DBH and tree height after stem biomass increment

    biomass_increment = max(below_resource_day_stem, 1.0d0)

    call tree_allometry (  &
                              ! *** Input ***
    dbh_heartwood        , &  ! Heartwood diameter (m)
    dbh_sapwood          , &  ! Sapwood diameter (m)
    tree_h               , &  ! Tree height (m)
    trunk_biomass        , &  ! Stand-level trunk biomass (g/tree)
    wood_rho(p)          , &  ! Wood density (g/cm3)
    ALM5(p)              , &  ! Sapwood diameter proportion (m sapwood/m dbh)
    biomass_increment    , &  ! Increment of trunk biomass (g/tree)
                              !
                              ! *** Output ***
    new_dbh_heartwood    , &  ! New heartwood diameter after biomass increment (m)
    new_dbh_sapwood      , &  ! New sapwood diameter after biomass increment (m)
    new_tree_h             &  ! New tree height after biomass increment (m)
    )

    !---------------------------------------------------------------------
    ! Hydraulic resitance after biomass increment
    !---------------------------------------------------------------------

    ! Case 1: Increment of fine root biomass (g root/tree)

    new_root_biomass = root_biomass + max(below_resource_day_root, 1.0d0*c_n_root(p)/c_n_stem(p))

    call plant_hydraulics (&
                              ! *** Input ***
    root_radius(p)       , &  ! Fine root radius (m)
    new_root_biomass     , &  ! New fine root biomass after biomass increment (g root/tree)
    root_area            , &  ! Stand-level root coverage area (m2 ground/tree)
    root_depth           , &  ! Rooting depth (m)
    root_density(p)      , &  ! Specific root density (fine root) (g root/m3 root)
    hk                   , &  ! Soil hydraulic conductance (mmol H2O/m/s/MPa)
    root_resist(p)       , &  ! Hydraulic resistivity of root tissue (MPa.s.g root/mmol H2O)
    dbh_heartwood        , &  ! Heartwood diameter (m)
    dbh_sapwood          , &  ! Sapwood diameter (m)
    k_sap(p)             , &  ! Stem hydraulic conductivity (kg H2O.m/m2 sapwood/s/MPa)
    tree_h               , &  ! Tree height (m)
    bole_h               , &  ! Bole height (m)
    tot_psi              , &  ! Total soil water potential (MPa)
    minlp(p)             , &  ! Minimum leaf water potential (MPa)
                              !
                              ! *** Input type ***
    universal            , &  ! Universal variables
                              !
                              ! *** Input/Output ***
    flux                   &  ! Flux variables
    )
    r_whole_root_increment = r_whole                                                      ! <= to be returned.

    ! Case 2: Increment of stem diameter (m)

    call plant_hydraulics (&
                              ! *** Input ***
    root_radius(p)       , &  ! Fine root radius (m)
    root_biomass         , &  ! Stand-level fine root biomass (g root/tree)
    root_area            , &  ! Stand-level root coverage area (m2 ground/tree)
    root_depth           , &  ! Rooting depth (m)
    root_density(p)      , &  ! Specific root density (fine root) (g root/m3 root)
    hk                   , &  ! Soil hydraulic conductance (mmol H2O/m/s/MPa)
    root_resist(p)       , &  ! Hydraulic resistivity of root tissue (MPa.s.g root/mmol H2O)
    new_dbh_heartwood    , &  ! New heartwood diameter after biomass increment (m)
    new_dbh_sapwood      , &  ! New sapwood diameter after biomass increment (m)
    k_sap(p)             , &  ! Stem hydraulic conductivity (kg H2O.m/m2 sapwood/s/MPa)
    tree_h               , &  ! Tree height (m)
    bole_h               , &  ! Bole height (m)
    tot_psi              , &  ! Total soil water potential (MPa)
    minlp(p)             , &  ! Minimum leaf water potential (MPa)
                              !
                              ! *** Input type ***
    universal            , &  ! Universal variables
                              !
                              ! *** Input/Output ***
    flux                   &  ! Flux variables
    )
    r_whole_d_increment = r_whole                                                         ! <= to be returned.

    ! Save leaf tempearture profile

    tleaf_out = tleaf

    end associate
  end subroutine spac_photosynthesis

end module mod_spac_photosynthesis
