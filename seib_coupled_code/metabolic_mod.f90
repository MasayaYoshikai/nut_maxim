  !-----------------------------------------------------------------------
  subroutine daily_production ( &
    )
    !
    ! !Description:
    ! Compute daily production for all trees
    !
    ! !Uses:
    use data_structure
    use time_counter
    use grid_status_current2
    use vegi_status_current1
    use vegi_status_current2
    use mod_spac_photosynthesis
    use mod_metabolism
    use mod_crown_morphology
    use mod_growth
    use mod_tree_allometry
    !
    ! !Argumetns:
    implicit none
    !
    ! !Local variables:
    real(8), parameter :: crit_rad_dir = 0.0d0        ! Critical direct radiation at canopy top in a day (W/m2)
    real(8), parameter :: crown_expand_rate = 1.05d0  ! Rate of crown diameter expansion (/year)
    integer :: no                                     ! Tree index
    integer :: hour                                   ! 0 - 23
    integer :: hour_of_year                           ! Hour of year
    real(8) :: seib_lai                               ! Leaf area index (m2 leaf/m2 ground)
    real(8) :: tree_h_limit                           ! Tree height limitation based on proximate trees (m)
    real(8) :: rad_dir                                ! Hourly direct radiation at canopy top (W/m2)
    real(8) :: max_rad_dir                            ! Maximum direct radiation at canopy top in a day (W/m2)
    real(8) :: available_c                            ! Daily available carbon for tree growth (mol C/tree/day)
    real(8) :: available_n                            ! Daily available nitrogen for tree growth (mol N/tree/day)
    real(8) :: remaining_c                            ! Remaining carbon for tree growth (mol C/tree/day)
    real(8) :: remaining_n                            ! Remaining nitrogen for tree growth (mol N/tree/day)
    real(8) :: d_trunk                                ! dTB/dt: Daily trunk biomass change (g trunk/tree/day)
    real(8) :: d_leaf                                 ! dLB/dt: Daily leaf biomass change (g leaf/tree/day)
    real(8) :: d_root                                 ! dRB/dt: Daily root biomass change (g root/tree/day)
    real(8) :: d_coarse_root                          ! dCRB/dt: Daily coarse root biomass change (g coarse root/tree/day)
    real(8) :: d_above_root                           ! dAGR/dt: Daily above-ground root biomass change (g/tree/day)
    real(8) :: d_stock                                ! dSTOCK/dt: Daily stock biomass change (g stock/tree/day)
    integer :: purge_flag                             ! Flag of perge crown bottom layer (0: no, 1: yes)
    integer :: n_limit_flag                           ! Flag of nitrogen limitation for tree growth (1: N-limited, 0: C-limited)
    integer :: expand_flag                            ! Flag of stem expansion (1: expand, 0: extend)
    real(8) :: new_dbh_heartwood                      ! New heartwood diameter after biomass increment (m)
    real(8) :: new_dbh_sapwood                        ! New sapwood diameter after biomass increment (m)
    real(8) :: new_tree_h                             ! New tree height after biomass increment (m)
    real(8) :: crown_d_max                            ! Potential crown diameter based on allometric relation (m)
    !
    ! !Local outputs from stand-level SPAC-Photosynthesis model:
    integer :: n_leaf_layer                       ! Number of leaf layers
    real(8) :: max_water_uptake_day               ! Maximum stand-level water uptake rate (mol H2O/tree/s)
    real(8) :: max_stand_et_day                   ! Maximum stand-level transpiration rate of the day (mol H2O/tree/s)
    real(8) :: sapflow_day                        ! Daily sapflow rate at stand-level (m3 H2O/tree/day)
    real(8) :: an_mean_day_max                    ! Daily maximum canopy mean leaf net photosynthesis (umol CO2/m2 leaf/s)
    real(8) :: an_top_day_max                     ! Daily maximum canopy top leaf net photosynthesis (umol CO2/m2 leaf/s)
    real(8) :: c_uptake_day                       ! Daily carbon uptake rate by photosynthesis (mol C/tree/day)
    real(8) :: c_uptake_bottom_day                ! Daily bottom layer carbon uptake rate by photosynthesis (mol CO2/m2 leaf/day)
    real(8) :: n_uptake_bottom_day                ! Daily bottom layer nitrogen uptake rate (mol N/m2 leaf/day)
    real(8) :: n_uptake_day                       ! Daily nitrogen uptake rate (mol N/tree/day)
    real(8) :: r_whole_root_increment             ! Whole plant resistance after root biomass increment (MPa.s.tree/mmol H2O)
    real(8) :: r_whole_d_increment                ! Whole plant resistance after diameter increment (MPa.s.tree/mmol H2O)
    real(8) :: r_root_out                         ! Root-stem resistance (MPa.s.tree/mmol H2O)
    real(8) :: r_sap_out                          ! Whole-plant stem resistance (MPa.s.tree/mmol H2O)
    real(8), dimension(Max_hgt) :: tleaf_out      ! Leaf temperature profile (K)
    !---------------------------------------------------------------------

    ! Maximum direct radiation at canopy top in a day (W/m2)

    max_rad_dir = 0.0d0
    do hour = 0, 23
       hour_of_year = (doy -1)*24 + hour + 1
       rad_dir = time_series_rad_dir(hour_of_year)
       max_rad_dir = max(max_rad_dir, rad_dir)
    end do

    ! Computation of all trees

    do no = 1, Max_no

       if ( .not. tree_exist(no) ) cycle
       if ( mass_leaf(no) <= 0.0 ) cycle

       ! Leaf area index of each tree (m2 leaf/m2 ground)

       seib_lai = mass_leaf(no) * SLA(pft(no)) / crown_area(no)

       ! Tree height limitation based on proximate trees (m)
       ! <= came from subroutine spatial_limitation

       tree_h_limit = real(height_limit(no)) * STEP + 1.3d0

       !---------------------------------------------------------------------
       ! SPAC-Photosynthesis model
       !---------------------------------------------------------------------

!write(*,*) 'tree_index', no
!write(*,*) 'doy', doy

       call spac_photosynthesis ( &
                                            ! *** Input ***
       o2air                           , &  ! Atmospheric O2 (mmol/mol)
       co2air                          , &  ! Atmospheric CO2 (umol/mol)
       hksat                           , &  ! Soil hydraulic conductivity at saturation (mm H2O/s)
       moist                           , &  ! Relative soil moisture (-)
       bsw                             , &  ! Soil layer Clapp and Hornberger "b" parameter (-)
       psisat                          , &  ! Soil matric potential at saturation (MPa)
       soil_t                          , &  ! Soil water temperature (K)
       sal                             , &  ! Pore-water salinity (mol/m3)
       din                             , &  ! DIN concentration in pore-water (mol N/m3) 1.e-03 for (umol/L) -> (mol/m3)
       root_filter                     , &  ! Filtering rate of salt at root surface (-)
       root_resist                     , &  ! Hydraulic resistivity of root tissue (MPa.s.g/mmol H2O)
       root_density                    , &  ! Specific root density (fine root) (g biomass/m3 root)
       root_radius                     , &  ! Fine root radius (m)
       root_depth                      , &  ! Rooting depth (m)
       k_sap                           , &  ! Stem hydraulic conductivity (kg H2O.m/m2 sapwood/s/MPa)
       wood_rho                        , &  ! Wood density (g/cm3)
       c_n_stem                        , &  ! C/N ratio in mol in stem
       c_n_root                        , &  ! C/N ratio in mol in root
       t_growth                        , &  ! Growth temperature (degree), mean temperature for 30 days before
       t_home                          , &  ! Home temperature (degree), mean maximum temperature of the warmest month
       minlp                           , &  ! Minimum leaf water potential (MPa)
       dleaf                           , &  ! Leaf dimension (m)
       gs_max                          , &  ! Maximum stomatal conductance (mol H2O/m2/s)
       vcmaxpft                        , &  ! Maximum carboxylation rate at 25C (umol/m2/s)
       iota                            , &  ! Stomatal water-use efficiency (umol CO2/ mol H2O)
       dble(C_in_drymass)              , &  ! Proportion of carbon in biomass (g C/g DW)
       dble(ALM5)                      , &  ! Sapwood diameter proportion (m sapwood/m dbh)
                                            !
                                            ! *** From SEIB-DGVM ***
       doy                             , &  ! Day of the year (1 - 365)
       no                              , &  ! Tree index
       pft(no)                         , &  ! Species index
       dble(par_direct_rel(no,:))      , &  ! Profile of relative intensity of direct PAR within canopy compared to canopy top (fraction)
       dble(par_diffuse_rel)           , &  ! Profile of relative intensity of diffused PAR within canopy compared to canopy top (fraction)
       dble(nir_direct_rel(no,:))      , &  ! Profile of relative intensity of direct NIR within canopy compared to canopy top (fraction)
       dble(nir_diffuse_rel)           , &  ! Profile of relative intensity of diffused NIR within canopy compared to canopy top (fraction)
       irleaf                          , &  ! Leaf absorbed longwave radiation for canopy layer (W/m2 leaf)
       dble(kbm_vis)                   , &  ! Scattering adjusted light extinction coefficient for VIS (-)
       dble(kbm_nir)                   , &  ! Scattering adjusted light extinction coefficient for NIR (-)
       rwind_profile                   , &  ! Relative wind speed profile to the canopy top (-)
       tree_h(no)                      , &  ! Tree height (m)
       dble(mass_root(no))             , &  ! Stand-level fine root biomass (g/tree)
       dble(seib_lai)                  , &  ! Leaf area index (m2 leaf/m2 ground)
       dble(crown_area(no))            , &  ! Stand-level leaf area (m2/tree)
       dble(dbh_heartwood(no))         , &  ! Heartwood diameter (m)
       dble(dbh_sapwood(no))           , &  ! Sapwood diameter (m)
       dble(mass_trunk(no))            , &  ! Stand-level trunk biomass (g/tree)
       height(no)                      , &  ! SEIB-DGVM specific tree height (the unit is STEP!!)
       bole(no)                        , &  ! SEIB-DGVM specific bole height (the unit is STEP!!)
       time_series_tair                , &  ! Time-series air temperature (K)
       time_series_pair                , &  ! Time-series air pressure (Pa)
       time_series_wind                , &  ! Time-series wind speed (m/s)
       time_series_eair                , &  ! Time-series vapor pressure in air (Pa)
       time_series_rad_dir             , &  ! Time-series direct radiation at canopy top (W/m2)
       time_series_rad_dif             , &  ! Time-series diffused radiation at canopy top (W/m2)
                                            !
                                            ! *** Output ***
       n_leaf_layer                    , &  ! Number of leaf layers
       max_water_uptake_day            , &  ! Maximum stand-level water uptake rate (mol H2O/tree/s)
       max_stand_et_day                , &  ! Maximum stand-level transpiration rate of the day (mol H2O/tree/s)
       sapflow_day                     , &  ! Daily sapflow rate at stand-level (m3 H2O/tree/day)
       an_mean_day_max                 , &  ! Daily maximum canopy mean leaf net photosynthesis (umol CO2/m2 leaf/s)
       an_top_day_max                  , &  ! Daily maximum canopy top leaf net photosynthesis (umol CO2/m2 leaf/s)
       c_uptake_day                    , &  ! Daily carbon uptake rate by photosynthesis (mol C/tree/day)
       c_uptake_bottom_day             , &  ! Daily bottom layer carbon uptake rate by photosynthesis (mol CO2/m2 leaf/day)
       n_uptake_bottom_day             , &  ! Daily bottom layer nitrogen uptake rate (mol N/m2 leaf/day)
       n_uptake_day                    , &  ! Daily nitrogen uptake rate (mol N/tree/day)
       r_whole_root_increment          , &  ! Whole plant resistance after root biomass increment (MPa.s.tree/mmol H2O)
       r_whole_d_increment             , &  ! Whole plant resistance after diameter increment (MPa.s.tree/mmol H2O)
       r_root_out                      , &  ! Root-stem resistance (MPa.s.tree/mmol H2O)
       r_sap_out                       , &  ! Whole-plant stem resistance (MPa.s.tree/mmol H2O)
       tleaf_out                         &  ! Leaf temperature profile (K)
       )

       ! Save outputs to be returned to module.

!       tleaf_all(no,:) = tleaf_out           ! Leaf temperature profile
       plant_water_uptake(no) = sapflow_day   ! Daily plant water uptake rate (m3 H2O/tree/day)

       ! Sapflow and transpiration potential
       ! Use value of the previous day if the weather condition is not good for transpiration.

       if (max_rad_dir > crit_rad_dir) then
          potential_sapflow(no) = max_water_uptake_day
          potential_stand_et(no) = max_stand_et_day
       end if

!write(*,*) 'sal (psu)', sal * 58.44d0 / 1000.0d0
!write(*,*) 'plant_water_uptake, m3/tree/day', water_uptake_day


       !---------------------------------------------------------------------
       ! Woody respiration
       !---------------------------------------------------------------------

       ! Daily available carbon for tree growth (mol C/tree/day)

       available_c = c_uptake_day

       call woody_respiration ( &
                                            ! *** Input ***
                                            !
                                            ! *** Input parameters ***
       main_resp_stem                  , &  ! Maintenance stem respiration rate (day-1)
       main_resp_root                  , &  ! Maintenance stem respiration rate (day-1)
       dble(C_in_drymass)              , &  ! Proportion of carbon in biomass (g C/g DW)
                                            !
                                            ! *** Local input ***
       available_c                     , &  ! Daily available carbon for tree growth (mol C/tree/day)
                                            !
                                            ! *** From SEIB-DGVM ***
       pft(no)                         , &  ! Species index (1: Rh, 2: Br)
       dble(mass_trunk(no))            , &  ! Stand-level trunk biomass (g trunk/tree)
       dble(mass_leaf(no))             , &  ! Stand-level leaf biomass (g leaf/tree)
       dble(mass_root(no))             , &  ! Stand-level fine root biomass (g root/tree)
       dble(mass_coarse_root(no))      , &  ! Stand-level coarse root biomass (g/tree)
       dble(mass_stock(no))            , &  ! Stand-level stock biomass (g stock/tree)
                                            !
                                            ! *** Output ***
       d_leaf                          , &  ! dLB/dt: Daily leaf biomass change (g leaf/tree/day)
       d_root                          , &  ! dRB/dt: Daily root biomass change (g root/tree/day)
       d_stock                         , &  ! dSTOCK/dt: Daily stock biomass change (g stock/tree/day)
       remaining_c                       &  ! Remaining carbon for tree growth (mol C/tree/day)
       )

       ! Update results

       available_c = remaining_c
       if (growth_calc) then
          mass_leaf(no) = max((mass_leaf(no) + d_leaf), 1.0d0)
          mass_root(no) = max((mass_root(no) + d_root), 1.0d0)
          mass_stock(no) = max((mass_stock(no) + d_stock), 1.0d0)
       end if
       net_production(no) = net_production(no) + d_leaf + d_root + d_stock

       ! Zero-out

       d_leaf = 0.0d0
       d_root = 0.0d0
       d_stock = 0.0d0

       ! Available carbon and nitrogen for tree growth (mol C or N/tree/day)
       ! Subtract growth respiration for carbon.

       available_c = available_c * (1.0d0 - grow_resp(pft(no)))
       available_n = n_uptake_day

       ! For monitoring tree

       if (no == monitor) then
          c_uptake_bottom_day_monitor = c_uptake_bottom_day
          n_uptake_bottom_day_monitor = n_uptake_bottom_day
          max_water_uptake_day_monitor = max_water_uptake_day
          max_stand_et_day_monitor = max_stand_et_day
          water_uptake_day_monitor = sapflow_day
          available_c_monitor = available_c
          available_n_monitor = available_n
          lai_monitor = seib_lai
       end if

       !---------------------------------------------------------------------
       ! Nutrient resorption from leaves to fall
       !---------------------------------------------------------------------

       if (age(no) > 2) then ! age filter of turnover

          available_n = available_n + &
                        ((mass_leaf(no) * leaf_turn(pft(no)) * C_in_drymass / 12.0d0) &  ! Leaves to fall in a day (mol C/tree/day)
                        / c_n_leaf(pft(no))) * leaf_resorp(pft(no))                      ! Resorbed N (mol N/tree/day)

       end if ! For age filter of turnover

       !---------------------------------------------------------------------
       ! Decide whether to purge the bottom layer or not
       !---------------------------------------------------------------------

       call crown_bottom_purge ( &
                                            ! *** Input ***
                                            !
                                            ! *** Input parameters ***
       leaf_turn                       , &  ! Leaf turnover rate (day-1)
       grow_resp                       , &  ! Growth respiration rate (g DW/g DW)
       c_n_leaf                        , &  ! C/N ratio in mol in leaf
       dble(LA_max)                    , &  ! Maximum leaf area index (m2 leaf/m2 ground)
       dble(SLA)                       , &  ! Specific leaf area (m2 leaf one-sided /g leaf)
       dble(C_in_drymass)              , &  ! Proportion of carbon in biomass (g C/g DW)
                                            !
                                            ! *** Local input ***
       c_uptake_bottom_day             , &  ! Daily bottom layer carbon uptake rate by photosynthesis (mol CO2/m2 leaf/day)
       n_uptake_bottom_day             , &  ! Daily bottom layer nitrogen uptake rate (mol N/m2 leaf/day)
       seib_lai                        , &  ! Leaf area index (m2 leaf/m2 ground)
                                            !
                                            ! *** From SEIB-DGVM ***
       counter                         , &  ! Days since the simulation start (days)
       no                              , &  ! Tree index
       pft(no)                         , &  ! Species index (1: Rh, 2: Br)
       dble(STEP)                      , &  ! Canopy layer thickness (m)
       height(no)                      , &  ! SEIB-DGVM specific tree height (unit is STEP!!)
       bole(no)                        , &  ! SEIB-DGVM specific bole height (unit is STEP!!)
                                            !
                                            ! *** In/Output ***
       gpp_bottom                      , &  ! Gross production of 1 m2 of leaves at crown bottom layer (g/m2 leaf)
                                            !
                                            ! *** Output ***
       purge_flag                        &  ! Flag of perge crown bottom layer (0: no, 1: yes)
       )

       ! Case to perge crown bottom layer
       if (purge_flag == 1) then

          if (growth_calc) then
             bole(no) = bole(no) + 1
             mass_leaf(no) = mass_leaf(no) - mass_leaf(no) / n_leaf_layer
             n_leaf_layer = n_leaf_layer - 1
          end if
          net_production(no) = net_production(no) - mass_leaf(no) / n_leaf_layer

          if (n_leaf_layer < 1) then
             write(*,*) 'ERROR: no leaf layer after crown purge'
          end if
       end if

       !---------------------------------------------------------------------
       ! Coarse root, leaf and fine root turnover
       !---------------------------------------------------------------------

       if (age(no) > 2) then ! age filter of turnover

       call turnover ( &
                                            ! *** Input ***
                                            !
                                            ! *** Input parameters ***
       coarse_root_turn                , &  ! Coarse root turnover rate (day-1)
       leaf_turn                       , &  ! Leaf turnover rate (day-1)
       root_turn                       , &  ! Root turnover rate (day-1)
       c_n_stem                        , &  ! C/N ratio in mol in stem
       c_n_leaf                        , &  ! C/N ratio in mol in leaf
       c_n_root                        , &  ! C/N ratio in mol in root
       dble(C_in_drymass)              , &  ! Proportion of carbon in biomass (g C/g DW)
                                            !
                                            ! *** Local input ***
       available_c                     , &  ! Daily available carbon for tree growth (mol C/tree/day)
       available_n                     , &  ! Daily available nitrogen for tree growth (mol N/tree/day)
                                            !
                                            ! *** From SEIB-DGVM ***
       pft(no)                         , &  ! Species index (1: Rh, 2: Br)
       dble(mass_coarse_root(no))      , &  ! Stand-level coarse root biomass (g/tree)
       dble(mass_leaf(no))             , &  ! Stand-level leaf biomass (g leaf/tree)
       dble(mass_root(no))             , &  ! Stand-level fine root biomass (g root/tree)
       dble(mass_stock(no))            , &  ! Stand-level stock biomass (g stock/tree)
       potential_sapflow(no)           , &  ! Potential of daily sapflow rate under normal weather condition (mol H2O/tree/day)
       potential_stand_et(no)          , &  ! Potential of daily transpiration rate under normal weather condition (mol H2O/tree/day)
                                            !
                                            ! *** Output ***
       d_coarse_root                   , &  ! dCRB/dt: Daily coarse root biomass change (g/tree/day)
       d_leaf                          , &  ! dLB/dt: Daily leaf biomass change (g leaf/tree/day)
       d_root                          , &  ! dRB/dt: Daily root biomass change (g root/tree/day)
       d_stock                         , &  ! dSTOCK/dt: Daily stock biomass change (g stock/tree/day)
       remaining_c                     , &  ! Remaining carbon for tree growth (mol C/tree/day)
       remaining_n                       &  ! Remaining nitrogen for tree growth (mol N/tree/day)
       )

       ! Update results

       available_c = remaining_c
       available_n = remaining_n
       if (growth_calc) then
          mass_coarse_root(no) = max((mass_coarse_root(no) + d_coarse_root), 1.0d0)
          mass_leaf(no) = max((mass_leaf(no) + d_leaf), 1.0d0)
          mass_root(no) = max((mass_root(no) + d_root), 1.0d0)
          mass_stock(no) = max((mass_stock(no) + d_stock), 1.0d0)
       end if
       net_production(no) = net_production(no) + d_coarse_root + d_leaf + d_root + d_stock

     end if ! For age filter of turnover

       ! Zero-out

       d_coarse_root = 0.0d0
       d_leaf = 0.0d0
       d_root = 0.0d0
       d_stock = 0.0d0

       ! When there is still available resource.

       if (min(remaining_c, remaining_n) > 0.0d0) then

          !---------------------------------------------------------------------
          ! Reload to stock biomass
          ! If stock_biomass is insufficient, a portion of available resource is
          ! invested to stock_biomass.
          !---------------------------------------------------------------------

          call stock_reload ( &
                                               ! *** Input ***
                                               !
                                               ! *** Input parameters ***
          stock_trunk_ratio               , &  ! Desirable stock / trunk ratio (g stock/g trunk)
          c_n_leaf                        , &  ! C/N ratio in mol in leaf
          dble(C_in_drymass)              , &  ! Proportion of carbon in biomass (g C/g DW)
                                               !
                                               ! *** Local input ***
          available_c                     , &  ! Daily available carbon for tree growth (mol C/tree/day)
          available_n                     , &  ! Daily available nitrogen for tree growth (mol N/tree/day)
                                               !
                                               ! *** From SEIB-DGVM ***
          pft(no)                         , &  ! Species index (1: Rh, 2: Br)
          dble(mass_trunk(no))            , &  ! Stand-level trunk biomass (g trunk/tree)
          dble(mass_coarse_root(no))      , &  ! Stand-level coarse root biomass (g/tree)
          dble(mass_stock(no))            , &  ! Stand-level stock biomass (g stock/tree)
                                               !
                                               ! *** Output ***
          d_stock                         , &  ! dSTOCK/dt: Daily stock biomass change (g stock/tree/day)
          remaining_c                     , &  ! Remaining carbon for tree growth (mol C/tree/day)
          remaining_n                       &  ! Remaining nitrogen for tree growth (mol N/tree/day)
          )

          ! Update results

          available_c = remaining_c
          available_n = remaining_n
          if (growth_calc) then
             mass_stock(no) = mass_stock(no) + d_stock
          end if
          net_production(no) = net_production(no) + d_stock

          ! Zero-out

          d_stock = 0.0d0

          ! Flag of nitrogen limitation for tree growth (1: N-limited, 0: C-limited)

          if (available_c / available_n > c_n_stem(pft(no))) then
             n_limit_flag = 1
          else
             n_limit_flag = 0
          end if

          !---------------------------------------------------------------------
          ! Below-ground limitation (N-limited) case
          !---------------------------------------------------------------------

          if (n_limit_flag == 1) then

             call belowground_limit ( &
                                                  ! *** Input ***
                                                  !
                                                  ! *** Input parameters ***
             optimum_dpai                    , &  ! Optimum layer leaf area index (m2 leaf in a crown disk/m2 ground)
             c_n_stem                        , &  ! C/N ratio in mol in stem
             c_n_leaf                        , &  ! C/N ratio in mol in leaf
             c_n_root                        , &  ! C/N ratio in mol in root
             fine_root_ratio                 , &  ! Fraction of fine root in below-ground biomass (-)
             stock_trunk_ratio               , &  ! Desirable stock / trunk ratio (g stock/g trunk)
             dble(SLA)                       , &  ! Specific leaf area (m2 leaf one-sided /g leaf)
             dble(DBH_limit)                 , &  ! Limitation of trunk diameter (m)
             dble(C_in_drymass)              , &  ! Proportion of carbon in biomass (g C/g DW)
                                                  !
                                                  ! *** Local input ***
             n_leaf_layer                    , &  ! Number of leaf layers
             available_c                     , &  ! Daily available carbon for tree growth (mol C/tree/day)
             available_n                     , &  ! Daily available nitrogen for tree growth (mol N/tree/day)
!             r_whole_root_increment          , &  ! Whole plant resistance after root biomass increment (MPa.s.tree/mmol H2O)
!             r_whole_d_increment             , &  ! Whole plant resistance after diameter increment (MPa.s.tree/mmol H2O)
             r_sap_out                       , &  ! Root-stem resistance (MPa.s.tree/mmol H2O)
             r_root_out                      , &  ! Whole-plant stem resistance (MPa.s.tree/mmol H2O)
             tree_h_limit                    , &  ! Tree height limitation based on proximate trees (m)
                                                  !
                                                  ! *** From SEIB-DGVM ***
             pft(no)                         , &  ! Species index (1: Rh, 2: Br)
             tree_h(no)                      , &  ! Tree height (m)
             dble(mass_trunk(no))            , &  ! Stand-level trunk biomass (g trunk/tree)
             dble(mass_coarse_root(no))      , &  ! Stand-level coarse root biomass (g/tree)
             dble(mass_leaf(no))             , &  ! Stand-level leaf biomass (g leaf/tree)
             dble(mass_root(no))             , &  ! Stand-level fine root biomass (g root/tree)
             dble(mass_stock(no))            , &  ! Stand-level stock biomass (g stock/tree)
             dble(crown_area(no))            , &  ! Stand-level leaf crown area (m2 ground/tree)
             dble(dbh_heartwood(no))         , &  ! Heartwood diameter (m)
             dble(dbh_sapwood(no))           , &  ! Sapwood diameter (m)
             potential_sapflow(no)           , &  ! Potential of daily sapflow rate under normal weather condition (mol H2O/tree/day)
             potential_stand_et(no)          , &  ! Potential of daily transpiration rate under normal weather condition (mol H2O/tree/day)
                                                  !
                                                  ! *** Output ***
             d_leaf                          , &  ! dLB/dt: Daily leaf biomass change (g leaf/tree/day)
             d_root                          , &  ! dRB/dt: Daily root biomass change (g root/tree/day)
             d_trunk                         , &  ! dTB/dt: Daily trunk biomass change (g trunk/tree/day)
             d_coarse_root                   , &  ! dCRB/dt: Daily coarse root biomass change (g/tree/day)
             d_stock                         , &  ! dSTOCK/dt: Daily stock biomass change (g stock/tree/day)
             expand_flag                       &  ! Flag of stem expansion (1: expand, 0: extend)
             )

             ! Calculate needed prop root biomass corresponding to the DBH
             ! based on allometric relationship.

             if (d_trunk > 0.0d0) then

                call proot_allometry ( &
                                                     ! *** Input ***
                                                     !
                                                     ! *** Input parameters ***
                wood_rho                        , &  ! Wood density (g/cm3)
                pr_s_a                          , &  ! Slope for the scaling factor of prop root system
                pr_s_b                          , &  ! Intercept for the scaling factor of prop root system
                pr_h_a                          , &  ! Slope for the maximum root height of prop root system
                pr_h_b                          , &  ! Intercept for the maximum root height of prop root system
                pr_d                            , &  ! Mean prop root diameter (m)
                                                     !
                                                     ! *** From SEIB-DGVM ***
                pft(no)                         , &  ! Species index (1: Rh, 2: Br)
                dble(dbh_heartwood(no))         , &  ! Heartwood diameter (m)
                dble(dbh_sapwood(no))           , &  ! Sapwood diameter (m)
                dble(mass_above_root(no))       , &  ! Stand-level above-ground root biomass (g/tree)
                                                     !
                                                     ! *** In/Output ***
                d_trunk                         , &  ! dTB/dt: Daily trunk biomass change (g trunk/tree/day)
                                                     !
                                                     ! *** Output ***
                d_above_root                      &  ! dAGR/dt: Daily above-ground root biomass change (g/tree/day)
                )

             end if

             ! Morphological change based on trunk biomass change

             if (d_trunk > 0.0d0) then

                call tree_allometry (  &
                                                    ! *** Input ***
                dble(dbh_heartwood(no))        , &  ! Heartwood diameter (m)
                dble(dbh_sapwood(no))          , &  ! Sapwood diameter (m)
                tree_h(no)                     , &  ! Tree height (m)
                dble(mass_trunk(no))           , &  ! Stand-level trunk biomass (g/tree)
                wood_rho(pft(no))              , &  ! Wood density (g/cm3)
                dble(ALM5(pft(no)))            , &  ! Sapwood diameter proportion (m sapwood/m dbh)
                d_trunk                        , &  ! Increment of trunk biomass (g/tree)
                                                    !
                                                    ! *** Output ***
                new_dbh_heartwood              , &  ! New heartwood diameter after biomass increment (m)
                new_dbh_sapwood                , &  ! New sapwood diameter after biomass increment (m)
                new_tree_h                       &  ! New tree height after biomass increment (m)
                )

                ! Update results
                ! Case plant transpiration ability is limiting sapflow rate.

                if (potential_stand_et(no) < potential_sapflow(no)) then

                   ! When tree can extend
                   ! => Tree height grows.

                   if (expand_flag == 0) then

                      ! Convert tree_h (m) to height (STEP)

                      if (growth_calc) then
                         tree_h(no) = new_tree_h
                         height(no) = int((tree_h(no) - 1.3) / STEP)
                         height(no) = max(2, height(no))
                      end if

                   ! When tree cannot extend due to allometric limitation
                   ! => Stem diameter expands.

                   else

                      if (growth_calc) then
                         dbh_heartwood(no) = new_dbh_heartwood
                         dbh_sapwood(no) = new_dbh_sapwood
                      end if

                   end if

                ! Case plant hydraulics is limiting sapflow rate.
                ! => Stem diameter expands.

                else

                   if (growth_calc) then
                      dbh_heartwood(no) = new_dbh_heartwood
                      dbh_sapwood(no) = new_dbh_sapwood
                   end if

                end if
             end if

             ! Update results

             if (growth_calc) then
                mass_leaf(no) = mass_leaf(no) + d_leaf
                mass_root(no) = mass_root(no) + d_root
                mass_trunk(no) = mass_trunk(no) + d_trunk
                mass_coarse_root(no) = mass_coarse_root(no) + d_coarse_root
                mass_above_root(no) = mass_above_root(no) + d_above_root
                mass_stock(no) = mass_stock(no) + d_stock
             end if
             net_production(no) = net_production(no) + d_coarse_root + d_above_root &
                                  + d_leaf + d_root + d_stock + d_trunk

          !---------------------------------------------------------------------
          ! Above-ground limitation (C-limited) case
          !---------------------------------------------------------------------

          else

             call aboveground_limit ( &
                                                  ! *** Input ***
                                                  !
                                                  ! *** Input parameters ***
             optimum_dpai                    , &  ! Optimum layer leaf area index (m2 leaf in a crown disk/m2 ground)
             c_n_stem                        , &  ! C/N ratio in mol in stem
             c_n_leaf                        , &  ! C/N ratio in mol in leaf
             stock_trunk_ratio               , &  ! Desirable stock / trunk ratio (g stock/g trunk)
             dble(SLA)                       , &  ! Specific leaf area (m2 leaf one-sided /g leaf)
             dble(DBH_limit)                 , &  ! Limitation of trunk diameter (m)
             dble(C_in_drymass)              , &  ! Proportion of carbon in biomass (g C/g DW)
                                                  !
                                                  ! *** Local input ***
             an_mean_day_max                 , &  ! Daily maximum canopy mean leaf net photosynthesis (umol CO2/m2 leaf/s)
             an_top_day_max                  , &  ! Daily maximum canopy top leaf net photosynthesis (umol CO2/m2 leaf/s)
             n_leaf_layer                    , &  ! Number of leaf layers
             available_c                     , &  ! Daily available carbon for tree growth (mol C/tree/day)
             available_n                     , &  ! Daily available nitrogen for tree growth (mol N/tree/day)
             tree_h_limit                    , &  ! Tree height limitation based on proximate trees (m)
                                                  !
                                                  ! *** From SEIB-DGVM ***
             pft(no)                         , &  ! Species index (1: Rh, 2: Br)
             tree_h(no)                      , &  ! Tree height (m)
             dble(mass_trunk(no))            , &  ! Stand-level trunk biomass (g trunk/tree)
             dble(mass_coarse_root(no))      , &  ! Stand-level coarse root biomass (g/tree)
             dble(mass_leaf(no))             , &  ! Stand-level leaf biomass (g leaf/tree)
             dble(mass_stock(no))            , &  ! Stand-level stock biomass (g stock/tree)
             dble(crown_area(no))            , &  ! Stand-level leaf crown area (m2 ground/tree)
             dble(dbh_heartwood(no))         , &  ! Heartwood diameter (m)
             dble(dbh_sapwood(no))           , &  ! Sapwood diameter (m)
                                                  !
                                                  ! *** Output ***
             d_leaf                          , &  ! dLB/dt: Daily leaf biomass change (g leaf/tree/day)
             d_trunk                         , &  ! dTB/dt: Daily trunk biomass change (g trunk/tree/day)
             d_stock                         , &  ! dSTOCK/dt: Daily stock biomass change (g stock/tree/day)
             expand_flag                       &  ! Flag of stem expansion (1: expand, 0: extend)
             )

             ! Calculate needed prop root biomass corresponding to the DBH
             ! based on allometric relationship.

             if (d_trunk > 0.0d0) then

                call proot_allometry ( &
                                                     ! *** Input ***
                                                     !
                                                     ! *** Input parameters ***
                wood_rho                        , &  ! Wood density (g/cm3)
                pr_s_a                          , &  ! Slope for the scaling factor of prop root system
                pr_s_b                          , &  ! Intercept for the scaling factor of prop root system
                pr_h_a                          , &  ! Slope for the maximum root height of prop root system
                pr_h_b                          , &  ! Intercept for the maximum root height of prop root system
                pr_d                            , &  ! Mean prop root diameter (m)
                                                     !
                                                     ! *** From SEIB-DGVM ***
                pft(no)                         , &  ! Species index (1: Rh, 2: Br)
                dble(dbh_heartwood(no))         , &  ! Heartwood diameter (m)
                dble(dbh_sapwood(no))           , &  ! Sapwood diameter (m)
                dble(mass_above_root(no))       , &  ! Stand-level above-ground root biomass (g/tree)
                                                     !
                                                     ! *** In/Output ***
                d_trunk                         , &  ! dTB/dt: Daily trunk biomass change (g trunk/tree/day)
                                                     !
                                                     ! *** Output ***
                d_above_root                      &  ! dAGR/dt: Daily above-ground root biomass change (g/tree/day)
                )

             end if

             ! Morphological change based on trunk biomass change

             if (d_trunk > 0.0d0) then

                call tree_allometry (  &
                                                    ! *** Input ***
                dble(dbh_heartwood(no))        , &  ! Heartwood diameter (m)
                dble(dbh_sapwood(no))          , &  ! Sapwood diameter (m)
                tree_h(no)                     , &  ! Tree height (m)
                dble(mass_trunk(no))           , &  ! Stand-level trunk biomass (g/tree)
                wood_rho(pft(no))              , &  ! Wood density (g/cm3)
                dble(ALM5(pft(no)))            , &  ! Sapwood diameter proportion (m sapwood/m dbh)
                d_trunk                        , &  ! Increment of trunk biomass (g/tree)
                                                    !
                                                    ! *** Output ***
                new_dbh_heartwood              , &  ! New heartwood diameter after biomass increment (m)
                new_dbh_sapwood                , &  ! New sapwood diameter after biomass increment (m)
                new_tree_h                       &  ! New tree height after biomass increment (m)
                )

                ! When tree can extend
                ! => Tree height grows.

                if (expand_flag == 0) then

                   ! Convert tree_h (m) to height (STEP)

                   if (growth_calc) then
                      tree_h(no) = new_tree_h
                      height(no) = int((tree_h(no) - 1.3) / STEP)
                      height(no) = max(2, height(no))
                   end if

                ! When tree cannot extend due to allometric limitation
                ! => Stem diameter expands.

                else

                   if (growth_calc) then
                      dbh_heartwood(no) = new_dbh_heartwood
                      dbh_sapwood(no) = new_dbh_sapwood
                   end if

                end if
             end if

             ! Update results

             if (growth_calc) then
                mass_leaf(no) = mass_leaf(no) + d_leaf
                mass_trunk(no) = mass_trunk(no) + d_trunk
                mass_above_root(no) = mass_above_root(no) + d_above_root
                mass_stock(no) = mass_stock(no) + d_stock
             end if
             net_production(no) = net_production(no) + d_leaf + d_stock + d_trunk + d_above_root

          end if
       end if

       !---------------------------------------------------------------------
       ! New crown diameter and area based on allometric relation
       ! (only when the crown is not suppressed by the neighbor trees)
       !---------------------------------------------------------------------

       ! Potential crown diameter based on allometric relation (m)

       crown_d_max = crown_allometry ( &
                                            ! *** Input ***
                                            !
                                            ! *** Input parameters ***
       crown_a                         , &  ! Tree crown allometric parameter
       crown_b                         , &  ! Tree crown allometric parameter
                                            !
                                            ! *** From SEIB-DGVM ***
       pft(no)                         , &  ! Species index (1: Rh, 2: Br)
       dble(dbh_heartwood(no))         , &  ! Heartwood diameter (m)
       dble(dbh_sapwood(no))             &  ! Sapwood diameter (m)
       )

       ! See whether the crown is suppressed or not.
       ! Only trees larger than a certain size could be suppressed.

       if (crown_d_max >= radius_limit(no) * 2.0d0 .and. &
          (dbh_heartwood(no)+dbh_sapwood(no)) >= 0.03d0) then

          if (.not. crown_limit_flag(no)) then
             if (growth_calc) then
                crown_diameter(no) = min(crown_diameter(no), radius_limit(no) * 2.0d0)
                crown_area(no) = (crown_diameter(no) * 0.5d0) * (crown_diameter(no) * 0.5d0) * PI
             end if
          end if

          crown_limit_flag(no) = .true.  ! The crown is suppressed by the neighbor trees.

       end if

       ! New crown diameter and area based on allometric relationship
       ! when the crown is not suppressed by the neighbor trees.

       if (.not. crown_limit_flag(no)) then  ! = .false.: The crown is not suppressed.

          if (growth_calc) then
             crown_diameter(no) = crown_d_max
             crown_area(no) = (crown_diameter(no) * 0.5d0) * (crown_diameter(no) * 0.5d0) * PI
          end if

       end if

       !---------------------------------------------------------------------
       ! Yearly crown expansion.
       ! !!! Only for trees whose crowns are suppressed. !!!
       !---------------------------------------------------------------------

       if (doy == Day_in_Year) then

          if (crown_limit_flag(no)) then  ! = .true.: The crown is suppressed.

              if (growth_calc) then
                 crown_diameter(no) = min(crown_d_max, crown_expand_rate * crown_diameter(no))  ! Potential crown diameter (m)
                 crown_diameter(no) = min(crown_diameter(no), radius_limit(no) * 2.0d0)         ! New crown diameter based on expansion rate and space availability
                 crown_area(no) = (crown_diameter(no) * 0.5d0) * (crown_diameter(no) * 0.5d0) * PI
              end if

          end if

       end if

       ! New leaf area per tree (m2 leaf/tree)

       la(no) = mass_leaf(no) * SLA(pft(no))

!write(*,*) 'trunk', mass_trunk(no), 'leaf', mass_leaf(no), 'root', mass_root(no), 'stock', mass_stock(no)

!write(*,*) 'year', year, 'p', PFT(no), 'dbh', dbh_heartwood(no)+dbh_sapwood(no), 'tree_h', tree_h(no), &
!           'crown_diameter', &
!           crown_diameter(no), 'crown_depth', (height(no) - bole(no))*STEP

    end do

  end subroutine daily_production
