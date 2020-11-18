module mod_metabolism

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
  ! !PUBLIC MEMBER FUNCTIONS
  public  :: woody_respiration
  public  :: turnover
  private :: turnover_x
  !
  ! !PRIVATE MEMBER FUNCTIONS
  !
  !-----------------------------------------------------------------------

  contains

  !-----------------------------------------------------------------------
  subroutine woody_respiration ( &
                                       ! *** Input ***
                                       !
                                       ! *** Input parameters ***
    main_resp_stem                , &  ! Maintenance stem respiration rate (day-1)
    main_resp_root                , &  ! Maintenance stem respiration rate (day-1)
    C_in_drymass                  , &  ! Proportion of carbon in biomass (g C/g DW)
                                       !
                                       ! *** Local input ***
    available_c                   , &  ! Daily available carbon for tree growth (mol C/tree/day)
                                       !
                                       ! *** From SEIB-DGVM ***
    p                             , &  ! Species index (1: Rh, 2: Br)
    trunk_biomass                 , &  ! Stand-level trunk biomass (g trunk/tree)
    leaf_biomass                  , &  ! Stand-level leaf biomass (g leaf/tree)
    root_biomass                  , &  ! Stand-level fine root biomass (g root/tree)
    coarse_root_biomass           , &  ! Stand-level coarse root biomass (g/tree)
    stock_biomass                 , &  ! Stand-level stock biomass (g stock/tree)
                                       !
                                       ! *** Output ***
    d_leaf                        , &  ! dLB/dt: Daily leaf biomass change (g leaf/tree/day)
    d_root                        , &  ! dRB/dt: Daily root biomass change (g root/tree/day)
    d_stock                       , &  ! dSTOCK/dt: Daily stock biomass change (g stock/tree/day)
    remaining_c                     &  ! Remaining carbon for tree growth (mol C/tree/day)
    )
    !
    ! !Description:
    ! Calculate carbon consumption by stem and root respiration
    ! (Leaf respiration rate was alraedy considered in SPAC-Photosynthesis model.)
    ! Step 1: available_c is used.
    ! Step 2: stock_biomass is used.
    ! Step 3: A portion of leaf and root is dropped (same procedure in SEIB-DGVM).
    !
    ! !Uses:
    !
    !
    ! !Argumetns:
    implicit none
    real(8), intent(in)  :: main_resp_stem(:), main_resp_root(:)
    real(8), intent(in)  :: C_in_drymass, available_c
    integer, intent(in)  :: p
    real(8), intent(in)  :: trunk_biomass, leaf_biomass, root_biomass
    real(8), intent(in)  :: coarse_root_biomass, stock_biomass
    real(8), intent(out) :: d_leaf, d_root, d_stock, remaining_c
    !
    ! !Local variables:
    real(8) :: required_c                              ! Required carbon for stem and root respiration (g C/tree/day)
    real(8) :: required_c_remaining                    ! Remaining required carbon for stem and root respiration (g C/tree/day)
    real(8), parameter :: frac_organ_remove = 0.01d0   ! Proportion of biomass to be dropped
    !---------------------------------------------------------------------

    ! Case 1: available_c can compensate woody respiration.
    ! Case 2: available_c is lacking.
    !         Case 2-1: stock_biomass can compensate remaining required respiration.
    !         Case 2-2: stock_biomass is also lacking.

    ! Available carbon for tree growth (g C/tree/day)
    ! Multiplying 12 for (mol C/tree/day) -> (g C/tree/day)
    ! Reversing unit will be done later.

    remaining_c = available_c * 12.0d0

    ! Required carbon for stem and root respiration (g C/tree/day)
    ! Multuplying C_in_drymass (g C/g DW) for (g DW/tree/day) -> (g C/tree/day)

    required_c = main_resp_stem(p) * (trunk_biomass + coarse_root_biomass) * C_in_drymass &
                 + main_resp_root(p) * root_biomass * C_in_drymass

    ! When daily gained resource is negative due to leaf respiration.

    if (remaining_c < 0.0d0) then
       required_c = required_c + abs(remaining_c)
       remaining_c = 0.0d0
    end if

    ! Case 1: remaining_c can compensate respiration.
    if (remaining_c > required_c) then

       ! Step 1: remaining_c is used.

       remaining_c = remaining_c - required_c
       required_c_remaining = 0.0d0
       d_leaf = 0.0d0
       d_root = 0.0d0
       d_stock = 0.0d0

    ! Case 2: remaining_c is lacking.
    else

       ! Step 1: remaining_c is used.

       required_c_remaining = required_c - remaining_c
       remaining_c = 0.0d0

       ! Case 2-1: stock_biomass can compensate remaining required respiration.
       if (stock_biomass * C_in_drymass > required_c_remaining) then

          ! Step 2: stock_biomass is used.

          d_stock = -1.0d0 * required_c_remaining / C_in_drymass
          required_c_remaining = 0.0d0
          d_leaf = 0.0d0
          d_root = 0.0d0

       ! Case 2-2: stock_biomass is also lacking.
       else

          ! Step 2: stock_biomass is used.

          d_stock = -1.0d0 * stock_biomass
          required_c_remaining = required_c_remaining - stock_biomass * C_in_drymass

          ! Step 3: A portion of leaf and root is dropped (same procedure in SEIB-DGVM).

          d_leaf = -1.0d0 * frac_organ_remove * leaf_biomass
          d_root = -1.0d0 * frac_organ_remove * root_biomass
          required_c_remaining = 0.0d0
          ! <= SEIB-DGVMによると、leafとrootの一部を切り落とせば呼吸の不足分は帳消しになるらしい。

       end if
    end if

    ! Unit conversion (g C/tree/day) -> (mol C/tree/day)

    remaining_c = remaining_c / 12.0d0

    if (required_c_remaining > 0.0d0) then
       write(*,*) 'ERROR: required carbon for respiration is still remaining'
    end if

  end subroutine woody_respiration

  !-----------------------------------------------------------------------
  subroutine turnover ( &
                                       ! *** Input ***
                                       !
                                       ! *** Input parameters ***
    coarse_root_turn              , &  ! Coarse root turnover rate (day-1)
    leaf_turn                     , &  ! Leaf turnover rate (day-1)
    root_turn                     , &  ! Root turnover rate (day-1)
    c_n_stem                      , &  ! C/N ratio in mol in stem
    c_n_leaf                      , &  ! C/N ratio in mol in leaf
    c_n_root                      , &  ! C/N ratio in mol in root
    C_in_drymass                  , &  ! Proportion of carbon in biomass (g C/g DW)
                                       !
                                       ! *** Local input ***
    available_c                   , &  ! Daily available carbon for tree growth (mol C/tree/day)
    available_n                   , &  ! Daily available nitrogen for tree growth (mol N/tree/day)
                                       !
                                       ! *** From SEIB-DGVM ***
    p                             , &  ! Species index (1: Rh, 2: Br)
    coarse_root_biomass           , &  ! Stand-level coarse root biomass (g/tree)
    leaf_biomass                  , &  ! Stand-level leaf biomass (g leaf/tree)
    root_biomass                  , &  ! Stand-level fine root biomass (g root/tree)
    stock_biomass                 , &  ! Stand-level stock biomass (g stock/tree)
    potential_sapflow             , &  ! Potential of daily sapflow rate under normal weather condition (mol H2O/tree/day)
    potential_stand_et            , &  ! Potential of daily transpiration rate under normal weather condition (mol H2O/tree/day)
                                       !
                                       ! *** Output ***
    d_coarse_root                 , &  ! dCRB/dt: Daily coarse root biomass change (g/tree/day)
    d_leaf                        , &  ! dLB/dt: Daily leaf biomass change (g leaf/tree/day)
    d_root                        , &  ! dRB/dt: Daily root biomass change (g root/tree/day)
    d_stock                       , &  ! dSTOCK/dt: Daily stock biomass change (g stock/tree/day)
    remaining_c                   , &  ! Remaining carbon for tree growth (mol C/tree/day)
    remaining_n                     &  ! Remaining nitrogen for tree growth (mol N/tree/day)
    )
    !
    ! !Description:
    ! Coarse root, leaf, and root turnover.
    ! Firstly, turnover of coarse root will be compensated.
    ! Then, lead or root is compensated depending on nutrient limitation and water uptake ability.
    !
    ! !Uses:
    !
    !
    ! !Argumetns:
    implicit none
    real(8), intent(in)  :: coarse_root_turn(:), leaf_turn(:), root_turn(:)
    real(8), intent(in)  :: c_n_stem(:), c_n_leaf(:), c_n_root(:), C_in_drymass
    real(8), intent(in)  :: available_c, available_n
    integer, intent(in)  :: p
    real(8), intent(in)  :: coarse_root_biomass, leaf_biomass, root_biomass, stock_biomass
    real(8), intent(in)  :: potential_sapflow, potential_stand_et
    real(8), intent(out) :: d_coarse_root, d_leaf, d_root, d_stock
    real(8), intent(out) :: remaining_c, remaining_n
    !
    ! !Local variables:
    integer :: n_limit_flag                      ! Flag of nitrogen limitation (1: N-limited, 0: C-limited)
    integer :: root_flag                         ! Flag of priority of turnover compensation (1: root, 0: leaf)
    real(8) :: coarse_root_c_required            ! Required carbon for coarse root turnover (mol C/tree/day)
    real(8) :: root_c_required                   ! Required carbon for root turnover (mol C/tree/day)
    real(8) :: leaf_c_required                   ! Required carbon for leaf turnover (mol C/tree/day)
    real(8) :: stock_remaining                   ! Remaining stock biomass (g stock/tree/day)
    !---------------------------------------------------------------------

    ! Remaining available carbon and nitrogen for tree growth (mol C or N/tree/day)

    remaining_c = available_c
    remaining_n = available_n

    ! Required carbon for turnover (mol C/tree/day)

    coarse_root_c_required = coarse_root_turn(p) * coarse_root_biomass * C_in_drymass / 12.0d0
    root_c_required = root_turn(p) * root_biomass * C_in_drymass / 12.0d0
    leaf_c_required = leaf_turn(p) * leaf_biomass * C_in_drymass / 12.0d0

    ! Prevent negative values.

    coarse_root_c_required = max(coarse_root_c_required, 0.0d0)
    root_c_required = max(root_c_required, 0.0d0)
    leaf_c_required = max(leaf_c_required, 0.0d0)

    ! Initialize variables.

    d_coarse_root = 0.0d0
    d_root = 0.0d0
    d_leaf = 0.0d0
    stock_remaining = stock_biomass

    !---------------------------------------------------------------------
    ! Coarse root turnover
    ! C/N/P ratio in stem is used for coarse root.
    !---------------------------------------------------------------------

    call turnover_x ( &
                                         ! *** Input ***
                                         !
                                         ! *** Input parameters ***
    c_n_stem(p)                     , &  ! C/N ratio in mol in x (x = coarse root, leaf, root)
    C_in_drymass                    , &  ! Proportion of carbon in biomass (g C/g DW)
                                         !
                                         ! *** Local input ***
    coarse_root_c_required          , &  ! Required carbon for turnover (mol C/tree/day)
                                         !
                                         ! *** In/Output ***
    remaining_c                     , &  ! Remaining available carbon for tree growth (mol C/tree/day)
    remaining_n                     , &  ! Remaining available nitrogen for tree growth (mol N/tree/day)
    d_coarse_root                   , &  ! dx/dt: Daily biomass change (x = coarse root, leaf, root) (g/tree/day)
    stock_remaining                   &  ! Remaining stock biomass (g stock/tree/day)
    )

    ! Flag of nitrogen limitation (1: N-limited, 0: C-limited)

    if (remaining_c / remaining_n > c_n_leaf(p)) then
       n_limit_flag = 1
    else
       n_limit_flag = 0
    end if

    ! Flag of priority of turnover compensation (1: root, 0: leaf)

    if (n_limit_flag == 1 .and. potential_sapflow <= potential_stand_et) then
       root_flag = 1
    else
       root_flag = 0
    end if

    ! When compensation for root turnover is priority.
    if (root_flag == 1) then

       !---------------------------------------------------------------------
       ! Fine root turnover
       !---------------------------------------------------------------------

       call turnover_x ( &
                                            ! *** Input ***
                                            !
                                            ! *** Input parameters ***
       c_n_root(p)                     , &  ! C/N ratio in mol in x (x = coarse root, leaf, root)
       C_in_drymass                    , &  ! Proportion of carbon in biomass (g C/g DW)
                                            !
                                            ! *** Local input ***
       root_c_required                 , &  ! Required carbon for turnover (mol C/tree/day)
                                            !
                                            ! *** In/Output ***
       remaining_c                     , &  ! Remaining available carbon for tree growth (mol C/tree/day)
       remaining_n                     , &  ! Remaining available nitrogen for tree growth (mol N/tree/day)
       d_root                          , &  ! dx/dt: Daily biomass change (x = coarse root, leaf, root) (g/tree/day)
       stock_remaining                   &  ! Remaining stock biomass (g stock/tree/day)
       )

       !---------------------------------------------------------------------
       ! Leaf root turnover
       !---------------------------------------------------------------------

       call turnover_x ( &
                                            ! *** Input ***
                                            !
                                            ! *** Input parameters ***
       c_n_leaf(p)                     , &  ! C/N ratio in mol in x (x = coarse root, leaf, root)
       C_in_drymass                    , &  ! Proportion of carbon in biomass (g C/g DW)
                                            !
                                            ! *** Local input ***
       leaf_c_required                 , &  ! Required carbon for turnover (mol C/tree/day)
                                            !
                                            ! *** In/Output ***
       remaining_c                     , &  ! Remaining available carbon for tree growth (mol C/tree/day)
       remaining_n                     , &  ! Remaining available nitrogen for tree growth (mol N/tree/day)
       d_leaf                          , &  ! dx/dt: Daily biomass change (x = coarse root, leaf, root) (g/tree/day)
       stock_remaining                   &  ! Remaining stock biomass (g stock/tree/day)
       )


    ! When compensation for leaf turnover is priority.
    else

       !---------------------------------------------------------------------
       ! Leaf root turnover
       !---------------------------------------------------------------------

       call turnover_x ( &
                                            ! *** Input ***
                                            !
                                            ! *** Input parameters ***
       c_n_leaf(p)                     , &  ! C/N ratio in mol in x (x = coarse root, leaf, root)
       C_in_drymass                    , &  ! Proportion of carbon in biomass (g C/g DW)
                                            !
                                            ! *** Local input ***
       leaf_c_required                 , &  ! Required carbon for turnover (mol C/tree/day)
                                            !
                                            ! *** In/Output ***
       remaining_c                     , &  ! Remaining available carbon for tree growth (mol C/tree/day)
       remaining_n                     , &  ! Remaining available nitrogen for tree growth (mol N/tree/day)
       d_leaf                          , &  ! dx/dt: Daily biomass change (x = coarse root, leaf, root) (g/tree/day)
       stock_remaining                   &  ! Remaining stock biomass (g stock/tree/day)
       )

       !---------------------------------------------------------------------
       ! Fine root turnover
       !---------------------------------------------------------------------

       call turnover_x ( &
                                            ! *** Input ***
                                            !
                                            ! *** Input parameters ***
       c_n_root(p)                     , &  ! C/N ratio in mol in x (x = coarse root, leaf, root)
       C_in_drymass                    , &  ! Proportion of carbon in biomass (g C/g DW)
                                            !
                                            ! *** Local input ***
       root_c_required                 , &  ! Required carbon for turnover (mol C/tree/day)
                                            !
                                            ! *** In/Output ***
       remaining_c                     , &  ! Remaining available carbon for tree growth (mol C/tree/day)
       remaining_n                     , &  ! Remaining available nitrogen for tree growth (mol N/tree/day)
       d_root                          , &  ! dx/dt: Daily biomass change (x = coarse root, leaf, root) (g/tree/day)
       stock_remaining                   &  ! Remaining stock biomass (g stock/tree/day)
       )

    end if

    ! Calculate d_stock (g stock/tree/day)

    d_stock = stock_remaining - stock_biomass

  end subroutine turnover

  !-----------------------------------------------------------------------
  subroutine turnover_x ( &
                                       ! *** Input ***
                                       !
                                       ! *** Input parameters ***
    c_n_ratio                     , &  ! C/N ratio in mol in x (x = coarse root, leaf, root)
    C_in_drymass                  , &  ! Proportion of carbon in biomass (g C/g DW)
                                       !
                                       ! *** Local input ***
    c_required                    , &  ! Required carbon for turnover (mol C/tree/day)
                                       !
                                       ! *** In/Output ***
    remaining_c                   , &  ! Remaining available carbon for tree growth (mol C/tree/day)
    remaining_n                   , &  ! Remaining available nitrogen for tree growth (mol N/tree/day)
    d_x                           , &  ! dx/dt: Daily biomass change (x = coarse root, leaf, root) (g/tree/day)
    stock_remaining                 &  ! Remaining stock biomass (g stock/tree/day)
    )
    !
    ! !Description:
    ! Turnover of component x (x = coarse root, leaf, root)
    ! Step 1: Gained resource is used.
    ! Step 2: stock_biomass is used.
    ! Step 3: Portion of remaining required resource is dropped.
    !
    ! !Uses:
    !
    !
    ! !Argumetns:
    implicit none
    real(8), intent(in)    :: c_n_ratio, C_in_drymass
    real(8), intent(in)    :: c_required
    real(8), intent(inout) :: remaining_c, remaining_n
    real(8), intent(inout) :: d_x, stock_remaining
    !
    ! !Local variables:
    real(8) :: c_required_remaining   ! Remaining required carbon for turnover (mol C/tree/day)
    !---------------------------------------------------------------------

    ! Step 1: Gained resource is used for turnover.

    ! When remaining resource can compensate for turnover.
    if (min(remaining_c, remaining_n * c_n_ratio) > c_required) then

       c_required_remaining = 0.0d0
       remaining_c = remaining_c - c_required
       remaining_n = remaining_n - c_required / c_n_ratio

       if (min(remaining_c, remaining_n) < 0.0d0) then
          write(*,*) 'ERROR: Negative remaining resource in Step 1.'
       end if

    ! When remaining resource is lacking for turnover
    else

       c_required_remaining = c_required - min(remaining_c, remaining_n * c_n_ratio)
       remaining_c = 0.0d0
       remaining_n = 0.0d0

    end if

    ! Step 2: remaining stock biomass is used for turnover.

    if (c_required_remaining > 0.0d0) then

       ! When remaining stock biomass can compensate for remaining turnover.
       if (stock_remaining * C_in_drymass / 12.0d0 > c_required_remaining) then

          c_required_remaining = 0.0d0
          stock_remaining = stock_remaining - c_required_remaining * 12.0d0 / C_in_drymass

          if (stock_remaining < 0.0d0) then
             write(*,*) 'ERROR: Negative remaining stock biomass in Step 2.'
          end if

       ! When remaining stock biomass is lacking for remaining turnover.
       else

          c_required_remaining = c_required_remaining - stock_remaining * C_in_drymass / 12.0d0
          stock_remaining = 0.0d0

       end if
    end if

    ! Step 3: Portion of remaining required resource is dropped.

    if (c_required_remaining > 0.0d0) then

       c_required_remaining = 0.0d0
       d_x = d_x - c_required_remaining * 12.0d0 / C_in_drymass

    end if

  end subroutine turnover_x

end module mod_metabolism
