module mod_growth

  !-----------------------------------------------------------------------
  ! !DESCRIPTION:
  ! Biomass investment based on tree growth optimization
  !
  ! !USES:
  !
  !
  ! !PUBLIC TYPES:
  implicit none
  !
  ! !PUBLIC MEMBER FUNCTIONS
  public  :: stock_reload
  public  :: belowground_limit
  public  :: aboveground_limit
  private :: invest_pattern1
  !
  ! !PRIVATE MEMBER FUNCTIONS:
  !
  ! !LOCAL VARIABLES
  private
  real(8), parameter :: invest_ratio = 0.3d0   ! Ratio of available resource to be invested to stock_biomass (-)
  real(8), parameter :: dpai_max = 1.1d0       ! Maximum ratio of dpai/dpai_optimum (-)
  !
  !-----------------------------------------------------------------------

  contains

  !-----------------------------------------------------------------------
  subroutine stock_reload ( &
                                       ! *** Input ***
                                       !
                                       ! *** Input parameters ***
    stock_trunk_ratio             , &  ! Desirable stock / trunk ratio (g stock/g trunk)
    c_n_leaf                      , &  ! C/N ratio in mol in leaf
    C_in_drymass                  , &  ! Proportion of carbon in biomass (g C/g DW)
                                       !
                                       ! *** Local input ***
    available_c                   , &  ! Daily available carbon for tree growth (mol C/tree/day)
    available_n                   , &  ! Daily available nitrogen for tree growth (mol N/tree/day)
                                       !
                                       ! *** From SEIB-DGVM ***
    p                             , &  ! Species index (1: Rh, 2: Br)
    trunk_biomass                 , &  ! Stand-level trunk biomass (g trunk/tree)
    coarse_root_biomass           , &  ! Stand-level coarse root biomass (g/tree)
    stock_biomass                 , &  ! Stand-level stock biomass (g stock/tree)
                                       !
                                       ! *** Output ***
    d_stock                       , &  ! dSTOCK/dt: Daily stock biomass change (g stock/tree/day)
    remaining_c                   , &  ! Remaining carbon for tree growth (mol C/tree/day)
    remaining_n                     &  ! Remaining nitrogen for tree growth (mol N/tree/day)
    )
    !
    ! !Description:
    ! Reload available resource to stock biomass to realize stock_trunk_ratio
    ! If stock biomass is lacking for stock_trunk_ratio, a portion of available
    ! resource (= invest_ratio) is invested to stock_biomass.
    !
    ! !Uses:
    !
    !
    ! !Argumetns:
    implicit none
    real(8), intent(in)  :: stock_trunk_ratio(:), c_n_leaf(:), C_in_drymass
    real(8), intent(in)  :: available_c, available_n
    integer, intent(in)  :: p
    real(8), intent(in)  :: trunk_biomass, coarse_root_biomass, stock_biomass
    real(8), intent(out) :: d_stock, remaining_c, remaining_n
    !
    ! !Local variables:
    real(8) :: optimum_stock                     ! Optimum stock biomass in terms of stock_trunk_ratio (g/tree)
    !---------------------------------------------------------------------

    ! Optimum stock biomass in terms of stock_trunk_ratio (g/tree)
    ! Note: It is considered that coarse root is part of trunk.

    optimum_stock = (trunk_biomass + coarse_root_biomass) * stock_trunk_ratio(p)

    ! Case 1: stock_biomass is sufficient
    if (stock_biomass >= optimum_stock) then

       d_stock = 0.0d0
       remaining_c = available_c
       remaining_n = available_n

    ! Case 2: stock_biomass is insufficient
    else

       d_stock = min(available_c, available_n * c_n_leaf(p)) * invest_ratio &
                 * 12.0d0 / C_in_drymass
       remaining_c = available_c - min(available_c, available_n * c_n_leaf(p)) * invest_ratio
       remaining_n = available_n - min(available_c / c_n_leaf(p), available_n) * invest_ratio

    end if

  end subroutine stock_reload

  !-----------------------------------------------------------------------
  subroutine belowground_limit ( &
                                       ! *** Input ***
                                       !
                                       ! *** Input parameters ***
    optimum_dpai                  , &  ! Optimum layer leaf area index (m2 leaf in a crown disk/m2 ground)
    c_n_stem                      , &  ! C/N ratio in mol in stem
    c_n_leaf                      , &  ! C/N ratio in mol in leaf
    c_n_root                      , &  ! C/N ratio in mol in root
    fine_root_ratio               , &  ! Fraction of fine root in below-ground biomass (-)
    stock_trunk_ratio             , &  ! Desirable stock / trunk ratio (g stock/g trunk)
    SLA                           , &  ! Specific leaf area (m2 leaf one-sided /g leaf)
    DBH_limit                     , &  ! Limitation of trunk diameter (m)
    C_in_drymass                  , &  ! Proportion of carbon in biomass (g C/g DW)
                                       !
                                       ! *** Local input ***
    n_leaf_layer                  , &  ! Number of leaf layers
    available_c                   , &  ! Daily available carbon for tree growth (mol C/tree/day)
    available_n                   , &  ! Daily available nitrogen for tree growth (mol N/tree/day)
    r_whole_root_increment        , &  ! Whole plant resistance after root biomass increment (MPa.s.tree/mmol H2O)
    r_whole_d_increment           , &  ! Whole plant resistance after diameter increment (MPa.s.tree/mmol H2O)
    tree_h_limit                  , &  ! Tree height limitation based on proximate trees (m)
                                       !
                                       ! *** From SEIB-DGVM ***
    p                             , &  ! Species index (1: Rh, 2: Br)
    tree_h                        , &  ! Tree height (m)
    trunk_biomass                 , &  ! Stand-level trunk biomass (g trunk/tree)
    coarse_root_biomass           , &  ! Stand-level coarse root biomass (g/tree)
    leaf_biomass                  , &  ! Stand-level leaf biomass (g leaf/tree)
    root_biomass                  , &  ! Stand-level fine root biomass (g root/tree)
    stock_biomass                 , &  ! Stand-level stock biomass (g stock/tree)
    crown_area                    , &  ! Stand-level leaf crown area (m2 ground/tree)
    dbh_heartwood                 , &  ! Heartwood diameter (m)
    dbh_sapwood                   , &  ! Sapwood diameter (m)
    potential_sapflow             , &  ! Potential of daily sapflow rate under normal weather condition (mol H2O/tree/day)
    potential_stand_et            , &  ! Potential of daily transpiration rate under normal weather condition (mol H2O/tree/day)
                                       !
                                       ! *** Output ***
    d_leaf                        , &  ! dLB/dt: Daily leaf biomass change (g leaf/tree/day)
    d_root                        , &  ! dRB/dt: Daily root biomass change (g root/tree/day)
    d_trunk                       , &  ! dTB/dt: Daily trunk biomass change (g trunk/tree/day)
    d_coarse_root                 , &  ! dCRB/dt: Daily coarse root biomass change (g/tree/day)
    d_stock                       , &  ! dSTOCK/dt: Daily stock biomass change (g stock/tree/day)
    expand_flag                     &  ! Flag of stem expansion (1: expand, 0: extend)
    )
    !
    ! !Description:
    ! Biomass investment in below-ground limitation case
    !
    ! !Uses:
    !
    !
    ! !Argumetns:
    implicit none
    real(8), intent(in)  :: optimum_dpai(:), c_n_stem(:), c_n_leaf(:), c_n_root(:)
    real(8), intent(in)  :: fine_root_ratio(:), stock_trunk_ratio(:)
    real(8), intent(in)  :: SLA(:), DBH_limit(:), C_in_drymass
    integer, intent(in)  :: n_leaf_layer
    real(8), intent(in)  :: available_c, available_n
    real(8), intent(in)  :: r_whole_root_increment, r_whole_d_increment, tree_h_limit
    integer, intent(in)  :: p
    real(8), intent(in)  :: tree_h, trunk_biomass, coarse_root_biomass, leaf_biomass, root_biomass
    real(8), intent(in)  :: stock_biomass, crown_area, dbh_heartwood, dbh_sapwood
    real(8), intent(in)  :: potential_sapflow, potential_stand_et
    real(8), intent(out) :: d_leaf, d_root, d_trunk, d_coarse_root, d_stock
    integer, intent(out) :: expand_flag
    !
    ! !Local variables:
    integer :: trunk_invest                  ! Flag of trunk investment (1: okay, 2: not okay)
    real(8) :: dpai                          ! Layer leaf area index (m2 leaf/m2 ground)
    real(8) :: leaf_required                 ! Required carbon investment to leaf to realize optimum_dpai (mol C/tree)
    real(8) :: optimum_coarse_root           ! Optimum coarse root biomass that realizes fine_root_ratio (g/tree)
    real(8) :: coarse_root_required          ! Required carbon investment to coarse root to realize optimum_coarse_root (mol C/tree)
    real(8) :: remaining_c                   ! Remaining available carbon for tree growth (mol C/tree/day)
    real(8) :: remaining_n                   ! Remaining available nitrogen for tree growth (mol N/tree/day)
    integer :: case_check                    ! Error check
    !---------------------------------------------------------------------

    ! Flag of trunk investment (1: okay, 0: not okay)
    ! - If stock biomass is low, plant does not invest biomass to trunk.
    ! - Note: Also if DBH or tree height exceeds the maximum limitation,
    !         plant does not invest biomass to trunk not to exceed the maximum.

    trunk_invest = 1
    if (stock_biomass < (trunk_biomass + coarse_root_biomass) * stock_trunk_ratio(p) * 0.5d0) then
       trunk_invest = 0
    end if

    ! Current layer leaf area index (m2 leaf/m2 ground)

    dpai = (leaf_biomass * SLA(p) / crown_area) / real(n_leaf_layer)

    ! Required investment of leaf biomass to realize optimum_dpai (mol C/tree)

    if (dpai < optimum_dpai(p)) then

       leaf_required = ((optimum_dpai(p) * crown_area * real(n_leaf_layer) / SLA(p)) &
                       - leaf_biomass) * C_in_drymass / 12.0d0

       if (leaf_required < 0.0d0) then
          write(*,*) 'ERROR: Negative leaf_required in belowground_limit.'
       end if

    else

       leaf_required = 0.0d0

    end if

    ! Optimum coarse root biomass that realizes fine_root_ratio (g/tree)

    optimum_coarse_root = root_biomass * (1.0d0 - fine_root_ratio(p)) / fine_root_ratio(p)

    ! Required biomass investment to coarse root to realize fine_root_ratio (mol C/tree)

    if (coarse_root_biomass < optimum_coarse_root) then

       coarse_root_required = (optimum_coarse_root - coarse_root_biomass) * C_in_drymass / 12.0d0

    else

       coarse_root_required = 0.0d0

    end if

    ! Remaining resource for tree growth (mol C or N/tree/day)

    remaining_c = available_c
    remaining_n = available_n

    ! Initialize variables

    d_trunk = 0.0d0
    d_coarse_root = 0.0d0
    d_leaf = 0.0d0
    d_root = 0.0d0
    d_stock = 0.0d0
    case_check = 0

    ! Initialize expand flag (1: expand, 0: extend)

    expand_flag = 1

    !---------------------------------------------------------------------
    ! When plant transpiration ability is limiting nutrient uptake rate,
    ! plant invests biomass to leaf to realize optimum_dpai.
    ! Then if resource is still remaining, plant invests residual resource
    ! to stem to be higher.
    !---------------------------------------------------------------------

    if (potential_stand_et < potential_sapflow) then

       call invest_pattern1 ( &
                                            ! *** Input ***
                                            !
                                            ! *** Input parameters ***
       optimum_dpai                    , &  ! Optimum layer leaf area index (m2 leaf in a crown disk/m2 ground)
       c_n_stem                        , &  ! C/N ratio in mol in stem
       c_n_leaf                        , &  ! C/N ratio in mol in leaf
       stock_trunk_ratio               , &  ! Desirable stock / trunk ratio (g stock/g trunk)
       DBH_limit                       , &  ! Limitation of trunk diameter (m)
       C_in_drymass                    , &  ! Proportion of carbon in biomass (g C/g DW)
                                            !
                                            ! *** Local input ***
       leaf_required                   , &  ! Required carbon investment to leaf to realize optimum_dpai (mol C/tree)
       dpai                            , &  ! Layer leaf area index (m2 leaf/m2 ground)
       trunk_invest                    , &  ! Flag of trunk investment (1: okay, 2: not okay)
       tree_h_limit                    , &  ! Tree height limitation based on proximate trees (m)
                                            !
                                            ! *** From SEIB-DGVM ***
       p                               , &  ! Species index (1: Rh, 2: Br)
       tree_h                          , &  ! Tree height (m)
       trunk_biomass                   , &  ! Stand-level trunk biomass (g trunk/tree)
       coarse_root_biomass             , &  ! Stand-level coarse root biomass (g/tree)
       stock_biomass                   , &  ! Stand-level stock biomass (g stock/tree)
       dbh_heartwood                   , &  ! Heartwood diameter (m)
       dbh_sapwood                     , &  ! Sapwood diameter (m)
                                            !
                                            ! *** In/Output ***
       remaining_c                     , &  ! Remaining available carbon for tree growth (mol C/tree/day)
       remaining_n                     , &  ! Remaining available nitrogen for tree growth (mol N/tree/day)
       d_leaf                          , &  ! dLB/dt: Daily leaf biomass change (g leaf/tree/day)
       d_trunk                         , &  ! dTB/dt: Daily trunk biomass change (g trunk/tree/day)
       d_stock                         , &  ! dSTOCK/dt: Daily stock biomass change (g stock/tree/day)
       case_check                      , &  ! Error check
       expand_flag                       &  ! Flag of stem expansion (1: expand, 0: extend)
       )

    end if

    !---------------------------------------------------------------------
    ! When plant water take ability is limiting nutrient uptake rate,
    ! plant invests biomass based on plant hydraulics optimization.
    !---------------------------------------------------------------------

    if (min(remaining_c, remaining_n) > 0.0d0) then

       ! Case 1: Investing to root is the most effective.

       if (r_whole_root_increment < r_whole_d_increment) then

          ! Case 1-1: However, coarse_root_biomass is lacking,
          !           and it is okay to invest to trunk.
          ! => All remaining resource goes to coarse root.

          if (coarse_root_required > 0.0d0 .and. trunk_invest == 1) then

             d_coarse_root = d_coarse_root + min(remaining_c, remaining_n * c_n_stem(p)) &
                             * 12.0d0 / C_in_drymass
             remaining_c = 0.0d0
             remaining_n = 0.0d0
             case_check = case_check + 1

          ! Case 1-2: Otherwise, all remaining resource goes to fine root.

          else

             d_root = d_root + min(remaining_c, remaining_n * c_n_root(p)) * 12.0d0 / C_in_drymass
             remaining_c = 0.0d0
             remaining_n = 0.0d0
             case_check = case_check + 1

          end if

       ! Case 2: Investing to trunk is the most effective.

       else

          ! Case 2-1: It is okay to invest to trunk, and DBH is not at limitation.
          ! => All remaining resource goes to trunk.

          if (trunk_invest == 1 .and. (dbh_heartwood + dbh_sapwood) <= DBH_limit(p)) then

             d_trunk = d_trunk + min(remaining_c, remaining_n * c_n_stem(p)) * 12.0d0 / C_in_drymass
             remaining_c = 0.0d0
             remaining_n = 0.0d0
             case_check = case_check + 1

          ! Case 2-2: It is okay to invest to trunk, but DBH is at limitation.
          !           And, coarse_root_biomass is lacking. <= This may rarely happen.
          ! => All remaining resource goes to coarse root.

          elseif (trunk_invest == 1 .and. coarse_root_required > 0.0d0) then

             d_coarse_root = d_coarse_root + min(remaining_c, remaining_n * c_n_stem(p)) &
                             * 12.0d0 / C_in_drymass
             remaining_c = 0.0d0
             remaining_n = 0.0d0
             case_check = case_check + 1

          ! Case 2-3: When it is not okay to invest to trunk based on remaining stock_biomass,
          !           or at DBH limitation, or no need to invest to coarse root,
          !           if stock_biomass is lacking:
          ! => All remaining resource goes to stock.

          elseif (stock_biomass < (trunk_biomass + coarse_root_biomass) &
                                  * stock_trunk_ratio(p)) then

             d_stock = d_stock + min(remaining_c, remaining_n * c_n_leaf(p)) * 12.0d0 / C_in_drymass
             remaining_c = 0.0d0
             remaining_n = 0.0d0
             case_check = case_check + 1

          ! Case 2-4: In other cases:
          ! => All remaining resource goes to root.

          else

             d_root = d_root + min(remaining_c, remaining_n * c_n_root(p)) * 12.0d0 / C_in_drymass
             remaining_c = 0.0d0
             remaining_n = 0.0d0
             case_check = case_check + 1

          end if
       end if
    end if

    ! Error check

    if (min(remaining_c, remaining_n) > 0.0d0) then
       write(*,*) 'ERROR: There is still remaining resource in belowground_limit.'
    end if
    if (case_check > 1) then
       write(*,*) 'ERROR: case_check is more than one in belowground_limit.'
    elseif (case_check < 1) then
       write(*,*) 'ERROR: case_check is still zero in belowground_limit'
    end if

  end subroutine belowground_limit

  !-----------------------------------------------------------------------
  subroutine aboveground_limit ( &
                                       ! *** Input ***
                                       !
                                       ! *** Input parameters ***
    optimum_dpai                  , &  ! Optimum layer leaf area index (m2 leaf in a crown disk/m2 ground)
    c_n_stem                      , &  ! C/N ratio in mol in stem
    c_n_leaf                      , &  ! C/N ratio in mol in leaf
    stock_trunk_ratio             , &  ! Desirable stock / trunk ratio (g stock/g trunk)
    SLA                           , &  ! Specific leaf area (m2 leaf one-sided /g leaf)
    DBH_limit                     , &  ! Limitation of trunk diameter (m)
    C_in_drymass                  , &  ! Proportion of carbon in biomass (g C/g DW)
                                       !
                                       ! *** Local input ***
    an_mean_day_max               , &  ! Daily maximum canopy mean leaf net photosynthesis (umol CO2/m2 leaf/s)
    an_top_day_max                , &  ! Daily maximum canopy top leaf net photosynthesis (umol CO2/m2 leaf/s)
    n_leaf_layer                  , &  ! Number of leaf layers
    available_c                   , &  ! Daily available carbon for tree growth (mol C/tree/day)
    available_n                   , &  ! Daily available nitrogen for tree growth (mol N/tree/day)
    tree_h_limit                  , &  ! Tree height limitation based on proximate trees (m)
                                       !
                                       ! *** From SEIB-DGVM ***
    p                             , &  ! Species index (1: Rh, 2: Br)
    tree_h                        , &  ! Tree height (m)
    trunk_biomass                 , &  ! Stand-level trunk biomass (g trunk/tree)
    coarse_root_biomass           , &  ! Stand-level coarse root biomass (g/tree)
    leaf_biomass                  , &  ! Stand-level leaf biomass (g leaf/tree)
    stock_biomass                 , &  ! Stand-level stock biomass (g stock/tree)
    crown_area                    , &  ! Stand-level leaf crown area (m2 ground/tree)
    dbh_heartwood                 , &  ! Heartwood diameter (m)
    dbh_sapwood                   , &  ! Sapwood diameter (m)
                                       !
                                       ! *** Output ***
    d_leaf                        , &  ! dLB/dt: Daily leaf biomass change (g leaf/tree/day)
    d_trunk                       , &  ! dTB/dt: Daily trunk biomass change (g trunk/tree/day)
    d_stock                       , &  ! dSTOCK/dt: Daily stock biomass change (g stock/tree/day)
    expand_flag                     &  ! Flag of stem expansion (1: expand, 0: extend)
    )
    !
    ! !Description:
    ! Biomass investment in above-ground limitation case
    !
    ! !Uses:
    !
    !
    ! !Argumetns:
    implicit none
    real(8), intent(in)  :: optimum_dpai(:), c_n_stem(:), c_n_leaf(:)
    real(8), intent(in)  :: stock_trunk_ratio(:), SLA(:), DBH_limit(:), C_in_drymass
    real(8), intent(in)  :: an_mean_day_max, an_top_day_max
    integer, intent(in)  :: n_leaf_layer
    real(8), intent(in)  :: available_c, available_n, tree_h_limit
    integer, intent(in)  :: p
    real(8), intent(in)  :: tree_h, trunk_biomass, coarse_root_biomass, leaf_biomass
    real(8), intent(in)  :: stock_biomass, crown_area, dbh_heartwood, dbh_sapwood
    real(8), intent(out) :: d_leaf, d_trunk, d_stock
    integer, intent(out) :: expand_flag
    !
    ! !Local variables:
    real(8), parameter :: an_saturated = 15.0d0    ! Light-saturated net photosynthesis rate (umol CO2/m2 leaf/s)
    real(8), parameter :: critical_ratio = 0.25d0  ! Critical ratio for determining low photosynthesis efficiency (fraction)
    integer :: efficiency_flag                     ! Flag of photosynthesis efficiency (1: Still efficient, 0: Inefficient)
    integer :: trunk_invest                        ! Flag of trunk investment (1: okay, 2: not okay)
    real(8) :: dpai                                ! Layer leaf area index (m2 leaf/m2 ground)
    real(8) :: leaf_required                       ! Required carbon investment to leaf to realize optimum_dpai (mol C/tree)
    real(8) :: remaining_c                         ! Remaining available carbon for tree growth (mol C/tree/day)
    real(8) :: remaining_n                         ! Remaining available nitrogen for tree growth (mol N/tree/day)
    integer :: case_check                          ! Error check

    real(8) :: required_leaf_invest             ! Required investment of leaf biomass to realize optimum_dpai (g/tree/day)
    real(8) :: remaining_resource               ! Remaining resource for tree growth (g/tree/day)
    real(8) :: daily_mean_an                    ! Daily mean CO2 assimilation rate over the canopy (umol CO2/m2 leaf/s)
    real(8), parameter :: critical_an = 2.0d0   ! Critical daily mean CO2 net asssimilation rate over canopy (umol CO2/m2 leaf/s)
    !---------------------------------------------------------------------

    !---------------------------------------------------------------------
    ! Flag of photosynthesis efficiency.
    ! 1: Still efficient => Low photosynthesis is due to small amount of leaves.
    ! 0: Inefficient => Low photosynthesis is due to low light availability (or intensity).
    !---------------------------------------------------------------------

    ! When daily maximum canopy (mean or top) leaf net photosynthesis
    ! is higher than XX % of light-saturated photosynthesis rate:
    ! => Photosynthesis is still efficient

    if (an_mean_day_max > critical_ratio * an_saturated) then

       efficiency_flag = 1

    ! Otherwise, photosynthesis is inefficient.

    else

       efficiency_flag = 0

    end if

    ! Flag of trunk investment (1: okay, 0: not okay)
    ! - If stock biomass is low, plant does not invest biomass to trunk.
    ! - Note: Also if tree height exceeds the maximum limitation,
    !         plant does not invest biomass to trunk not to exceed the maximum.

    trunk_invest = 1
    if (stock_biomass < (trunk_biomass + coarse_root_biomass) * stock_trunk_ratio(p) * 0.5d0) then
       trunk_invest = 0
    end if

    ! Current layer leaf area index (m2 leaf/m2 ground)

    dpai = (leaf_biomass * SLA(p) / crown_area) / real(n_leaf_layer)

    ! Required investment of leaf biomass to realize optimum_dpai (mol C/tree)

    if (dpai < optimum_dpai(p)) then

       leaf_required = ((optimum_dpai(p) * crown_area * real(n_leaf_layer) / SLA(p)) &
                       - leaf_biomass) * C_in_drymass / 12.0d0

       if (leaf_required < 0.0d0) then
          write(*,*) 'ERROR: Negative leaf_required in belowground_limit.'
       end if

    else

       leaf_required = 0.0d0

    end if

    ! Remaining resource for tree growth (mol C or N/tree/day)

    remaining_c = available_c
    remaining_n = available_n

    ! Initialize variables

    d_trunk = 0.0d0
    d_leaf = 0.0d0
    d_stock = 0.0d0
    case_check = 0

    ! Initialize expand flag (1: expand, 0: extend)

    expand_flag = 0

    !---------------------------------------------------------------------
    ! When low productivity is due to small amount of leaves,
    ! plant invests biomass to leaf to realize optimum_dpai.
    ! Then if resource is still remaining, plant invests residual resource
    ! to stem to be higher.
    !---------------------------------------------------------------------

    if (efficiency_flag == 1) then

       call invest_pattern1 ( &
                                            ! *** Input ***
                                            !
                                            ! *** Input parameters ***
       optimum_dpai                    , &  ! Optimum layer leaf area index (m2 leaf in a crown disk/m2 ground)
       c_n_stem                        , &  ! C/N ratio in mol in stem
       c_n_leaf                        , &  ! C/N ratio in mol in leaf
       stock_trunk_ratio               , &  ! Desirable stock / trunk ratio (g stock/g trunk)
       DBH_limit                       , &  ! Limitation of trunk diameter (m)
       C_in_drymass                    , &  ! Proportion of carbon in biomass (g C/g DW)
                                            !
                                            ! *** Local input ***
       leaf_required                   , &  ! Required carbon investment to leaf to realize optimum_dpai (mol C/tree)
       dpai                            , &  ! Layer leaf area index (m2 leaf/m2 ground)
       trunk_invest                    , &  ! Flag of trunk investment (1: okay, 2: not okay)
       tree_h_limit                    , &  ! Tree height limitation based on proximate trees (m)
                                            !
                                            ! *** From SEIB-DGVM ***
       p                               , &  ! Species index (1: Rh, 2: Br)
       tree_h                          , &  ! Tree height (m)
       trunk_biomass                   , &  ! Stand-level trunk biomass (g trunk/tree)
       coarse_root_biomass             , &  ! Stand-level coarse root biomass (g/tree)
       stock_biomass                   , &  ! Stand-level stock biomass (g stock/tree)
       dbh_heartwood                   , &  ! Heartwood diameter (m)
       dbh_sapwood                     , &  ! Sapwood diameter (m)
                                            !
                                            ! *** In/Output ***
       remaining_c                     , &  ! Remaining available carbon for tree growth (mol C/tree/day)
       remaining_n                     , &  ! Remaining available nitrogen for tree growth (mol N/tree/day)
       d_leaf                          , &  ! dLB/dt: Daily leaf biomass change (g leaf/tree/day)
       d_trunk                         , &  ! dTB/dt: Daily trunk biomass change (g trunk/tree/day)
       d_stock                         , &  ! dSTOCK/dt: Daily stock biomass change (g stock/tree/day)
       case_check                      , &  ! Error check
       expand_flag                       &  ! Flag of stem expansion (1: expand, 0: extend)
       )

    !---------------------------------------------------------------------
    ! When low productivity is due to lack of light,
    ! plant invests biomass to stem to be higher.
    !---------------------------------------------------------------------

    else

       ! Case 1: It is okay to invest to trunk, and tree height is not at limitation.
       ! => All remaining resource goes to trunk to be higher.

       if (trunk_invest == 1 .and. tree_h < tree_h_limit) then

          d_trunk = d_trunk + min(remaining_c, remaining_n * c_n_stem(p)) * 12.0d0 / C_in_drymass
          remaining_c = 0.0d0
          remaining_n = 0.0d0
          case_check = case_check + 1
          expand_flag = 0

       ! Case 2: It is okay to invest to trunk, but tree height is at limitation,
       !         and DBH is not at limitation.
       ! => All remaining resource goes to trunk to be thicker.

       elseif (trunk_invest == 1 .and. (dbh_heartwood + dbh_sapwood) <= DBH_limit(p)) then

          d_trunk = d_trunk + min(remaining_c, remaining_n * c_n_stem(p)) * 12.0d0 / C_in_drymass
          remaining_c = 0.0d0
          remaining_n = 0.0d0
          case_check = case_check + 1
          expand_flag = 1

       ! Case 3: It is not okay to invest to trunk due to lack of stock biomass or allometric limitation,
       !         and stock_biomass is lacking.
       ! => All remaining resource goes to stock.

       elseif (stock_biomass < (trunk_biomass + coarse_root_biomass) &
                               * stock_trunk_ratio(p)) then

          d_stock = d_stock + min(remaining_c, remaining_n * c_n_leaf(p)) * 12.0d0 / C_in_drymass
          remaining_c = 0.0d0
          remaining_n = 0.0d0
          case_check = case_check + 1

       ! Case 4: It is not okay to invest to trunk due to lack of stock biomass or allometric limitation,
       !         and stock_biomass is enough,
       !         and dpai is not reached to maximum.
       ! => All remaining resource goes to stock.

       elseif (dpai > optimum_dpai(p) * dpai_max) then

          d_stock = d_stock + min(remaining_c, remaining_n * c_n_leaf(p)) * 12.0d0 / C_in_drymass
          remaining_c = 0.0d0
          remaining_n = 0.0d0
          case_check = case_check + 1

       ! Case 5: It is not okay to invest to trunk due to lack of stock biomass or allometric limitation,
       !         and stock_biomass is enough,
       !         and doai is not reaached to maximum.
       ! => All remaining resource goes to leaf.

       else

          d_leaf = d_leaf + min(remaining_c, remaining_n * c_n_leaf(p)) * 12.0d0 / C_in_drymass
          remaining_c = 0.0d0
          remaining_n = 0.0d0
          case_check = case_check + 1

       end if
    end if

    ! Error check

    if (min(remaining_c, remaining_n) > 0.0d0) then
       write(*,*) 'ERROR: There is still remaining resource in aboveground_limit.'
    end if
    if (case_check > 1) then
       write(*,*) 'ERROR: case_check is more than one in aboveground_limit.'
    elseif (case_check < 1) then
       write(*,*) 'ERROR: case_check is still zero in aboveground_limit'
    end if

  end subroutine aboveground_limit

  !-----------------------------------------------------------------------
  subroutine invest_pattern1 ( &
                                       ! *** Input ***
                                       !
                                       ! *** Input parameters ***
    optimum_dpai                  , &  ! Optimum layer leaf area index (m2 leaf in a crown disk/m2 ground)
    c_n_stem                      , &  ! C/N ratio in mol in stem
    c_n_leaf                      , &  ! C/N ratio in mol in leaf
    stock_trunk_ratio             , &  ! Desirable stock / trunk ratio (g stock/g trunk)
    DBH_limit                     , &  ! Limitation of trunk diameter (m)
    C_in_drymass                  , &  ! Proportion of carbon in biomass (g C/g DW)
                                       !
                                       ! *** Local input ***
    leaf_required                 , &  ! Required carbon investment to leaf to realize optimum_dpai (mol C/tree)
    dpai                          , &  ! Layer leaf area index (m2 leaf/m2 ground)
    trunk_invest                  , &  ! Flag of trunk investment (1: okay, 2: not okay)
    tree_h_limit                  , &  ! Tree height limitation based on proximate trees (m)
                                       !
                                       ! *** From SEIB-DGVM ***
    p                             , &  ! Species index (1: Rh, 2: Br)
    tree_h                        , &  ! Tree height (m)
    trunk_biomass                 , &  ! Stand-level trunk biomass (g trunk/tree)
    coarse_root_biomass           , &  ! Stand-level coarse root biomass (g/tree)
    stock_biomass                 , &  ! Stand-level stock biomass (g stock/tree)
    dbh_heartwood                 , &  ! Heartwood diameter (m)
    dbh_sapwood                   , &  ! Sapwood diameter (m)
                                       !
                                       ! *** In/Output ***
    remaining_c                   , &  ! Remaining available carbon for tree growth (mol C/tree/day)
    remaining_n                   , &  ! Remaining available nitrogen for tree growth (mol N/tree/day)
    d_leaf                        , &  ! dLB/dt: Daily leaf biomass change (g leaf/tree/day)
    d_trunk                       , &  ! dTB/dt: Daily trunk biomass change (g trunk/tree/day)
    d_stock                       , &  ! dSTOCK/dt: Daily stock biomass change (g stock/tree/day)
    case_check                    , &  ! Error check
    expand_flag                     &  ! Flag of stem expansion (1: expand, 0: extend)
    )
    !
    ! !Description:
    ! Biomass investment pattern 1
    !
    ! !Uses:
    !
    !
    ! !Argumetns:
    implicit none
    real(8), intent(in)    :: optimum_dpai(:), c_n_stem(:), c_n_leaf(:)
    real(8), intent(in)    :: stock_trunk_ratio(:), DBH_limit(:), C_in_drymass
    real(8), intent(in)    :: leaf_required, dpai
    integer, intent(in)    :: trunk_invest
    real(8), intent(in)    :: tree_h_limit
    integer, intent(in)    :: p
    real(8), intent(in)    :: tree_h, trunk_biomass, coarse_root_biomass, stock_biomass
    real(8), intent(in)    :: dbh_heartwood, dbh_sapwood
    real(8), intent(inout) :: remaining_c, remaining_n
    real(8), intent(inout) :: d_leaf, d_trunk, d_stock
    integer, intent(inout) :: case_check
    integer, intent(inout) :: expand_flag
    !
    ! !Local variables:
    !
    !
    !---------------------------------------------------------------------

    ! Case 1: available resource can compensate leaf_required.

    if (min(remaining_c, remaining_n * c_n_leaf(p)) > leaf_required) then

       d_leaf = d_leaf + leaf_required * 12.0d0 / C_in_drymass   ! (mol C/tree) -> (g C/tree) -> (g DW/tree)
       remaining_c = remaining_c - leaf_required
       remaining_n = remaining_n - leaf_required / c_n_leaf(p)

       if (min(remaining_c, remaining_n) < 0.0d0) then
          write(*,*) 'ERROR: Negative remaining resource in Case 1 in belowground_limit.'
       end if

       ! Case 1-1: It is okay to invest to trunk, and tree height is not at limitation.
       ! => Residual goes to trunk to be higher.

       if (trunk_invest == 1 .and. tree_h < tree_h_limit) then

          d_trunk = d_trunk + min(remaining_c, remaining_n * c_n_stem(p)) * 12.0d0 / C_in_drymass
          remaining_c = 0.0d0
          remaining_n = 0.0d0
          case_check = case_check + 1
          expand_flag = 0

        ! Case 1-2: It is okay to invest to trunk, but tree height is at limitation,
        !           and DBH is not at limitation.
        ! => Residual goes to trunk to be thicker.

        elseif (trunk_invest == 1 .and. (dbh_heartwood + dbh_sapwood) <= DBH_limit(p)) then

          d_trunk = d_trunk + min(remaining_c, remaining_n * c_n_stem(p)) * 12.0d0 / C_in_drymass
          remaining_c = 0.0d0
          remaining_n = 0.0d0
          case_check = case_check + 1
          expand_flag = 1

       ! Case 1-3: It is not okay to invest to trunk due to lack of stock biomass or allometric limitation,
       !           and stock_biomass is lacking.
       ! => Residual goes to stock.

       elseif (stock_biomass < (trunk_biomass + coarse_root_biomass) &
                               * stock_trunk_ratio(p)) then

          d_stock = d_stock + min(remaining_c, remaining_n * c_n_leaf(p)) * 12.0d0 / C_in_drymass
          remaining_c = 0.0d0
          remaining_n = 0.0d0
          case_check = case_check + 1

       ! Case 1-4: It is not okay to invest to trunk due to lack of stock biomass or allometric limitation,
       !           and stock_biomass is enough,
       !           but dpai already reached maximum.
       ! => Residual goes to stock.

       elseif (dpai > optimum_dpai(p) * dpai_max) then

          d_stock = d_stock + min(remaining_c, remaining_n * c_n_leaf(p)) * 12.0d0 / C_in_drymass
          remaining_c = 0.0d0
          remaining_n = 0.0d0
          case_check = case_check + 1

       ! Case 1-5: It is not okay to invest to trunk due to lack of stock biomass or allometric limitation,
       !           and stock_biomass is enough,
       !           and dpai is not reached to maximum.
       ! => Residual goes to leaf.

       else

          d_leaf = d_leaf + min(remaining_c, remaining_n * c_n_leaf(p)) * 12.0d0 / C_in_drymass
          remaining_c = 0.0d0
          remaining_n = 0.0d0
          case_check = case_check + 1

       end if

    ! Case 2: available resource is lacking for leaf_required.
    ! => All remaining resource goes to leaf.

    else

       d_leaf = d_leaf + min(remaining_c, remaining_n * c_n_leaf(p)) * 12.0d0 / C_in_drymass
       remaining_c = 0.0d0
       remaining_n = 0.0d0
       case_check = case_check + 1

    end if

  end subroutine invest_pattern1

end module mod_growth
