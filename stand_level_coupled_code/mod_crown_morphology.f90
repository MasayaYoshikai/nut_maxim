module mod_crown_morphology

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
  public :: crown_bottom_purge
  public :: crown_allometry
  !
  ! !PRIVATE MEMBER FUNCTIONS
  !
  !-----------------------------------------------------------------------

  contains

  !-----------------------------------------------------------------------
  subroutine crown_bottom_purge ( &
                                       ! *** Input ***
                                       !
                                       ! *** Input parameters ***
    leaf_turn                     , &  ! Leaf turnover rate (day-1)
    grow_resp                     , &  ! Growth respiration rate (g DW/g DW)
    c_n_leaf                      , &  ! C/N ratio in mol in leaf
    LA_max                        , &  ! Maximum leaf area index (m2 leaf/m2 ground)
    SLA                           , &  ! Specific leaf area (m2 leaf one-sided /g leaf)
    C_in_drymass                  , &  ! Proportion of carbon in biomass (g C/g DW)
                                       !
                                       ! *** Local input ***
    c_uptake_bottom_day           , &  ! Daily bottom layer carbon uptake rate by photosynthesis (mol CO2/m2 leaf/day)
    n_uptake_bottom_day           , &  ! Daily bottom layer nitrogen uptake rate (mol N/m2 leaf/day)
    seib_lai                      , &  ! Leaf area index (m2 leaf/m2 ground)
                                       !
                                       ! *** From SEIB-DGVM ***
    day_count                     , &  ! Days since the simulation start (days)
    no                            , &  ! Tree index
    p                             , &  ! Species index (1: Rh, 2: Br)
    STEP                          , &  ! Canopy layer thickness (m)
    seib_height                   , &  ! SEIB-DGVM specific tree height (unit is STEP!!)
    seib_bole                     , &  ! SEIB-DGVM specific bole height (unit is STEP!!)
                                       !
                                       ! *** In/Output ***
    gpp_bottom                    , &  ! Gross production of 1 m2 of leaves at crown bottom layer (g/m2 leaf)
                                       !
                                       ! *** Output ***
    purge_flag                      &  ! Flag of perge crown bottom layer (0: no, 1: yes)
    )
    !
    ! !Description:
    ! Calculate net production of 1 m2 of leaves at crown bottom layer.
    ! Then decide whether purge the bottom layer or not.
    ! Also, if LAI is already reached at the maximum, the bottom layer will be purged.
    !
    ! !Uses:
    !
    !
    ! !Argumetns:
    implicit none
    real(8), intent(in)    :: leaf_turn(:), grow_resp(:), c_n_leaf(:), LA_max(:), SLA(:), C_in_drymass
    real(8), intent(in)    :: c_uptake_bottom_day, n_uptake_bottom_day, seib_lai
    integer, intent(in)    :: day_count, no, p
    real(8), intent(in)    :: STEP
    integer, intent(in)    :: seib_height, seib_bole
    real(8), intent(inout) :: gpp_bottom(:)
    integer, intent(out)   :: purge_flag
    !
    ! !Local variables:
    integer, parameter :: crown_day = 182               ! Interval of days to purge crown bottom layer
    real(8) :: turn_cost                               ! Turnover cost of 1 m2 of leaves (g/m2 leaf)
    real(8) :: leaf_production_c                       ! Leaf biomass production from C uptake rate per unit leaf area (g leaf/m2 leaf/day)
    real(8) :: leaf_production_n                       ! Leaf biomass production from N uptake rate per unit leaf area (g leaf/m2 leaf/day)
    real(8) :: crown_depth                             ! Current crown depth (m)
    real(8) :: crown_depth_min = 0.5d0                 ! Minimum crown  depth (m)
    real(8), parameter :: crit_ratio = 1.2d0           ! Critical ratio to purge crown bottom layer (gpp_bottom / turn_cost)
    !---------------------------------------------------------------------

    ! Turnover cost of 1 m2 of leaves (g leaf/m2 leaf)

    turn_cost = (1.0d0 / SLA(p)) * leaf_turn(p) * real(crown_day)

    ! Leaf biomass production from C and N uptake rates per unit leaf area (g leaf/m2 leaf/day)
    ! For carbon, growth respiration is subtracted.

    leaf_production_c = c_uptake_bottom_day * (1.0d0 - grow_resp(p)) * 12.0d0 / C_in_drymass
    leaf_production_n = n_uptake_bottom_day * c_n_leaf(p) * 12.0d0 / C_in_drymass

    ! Gross production of 1 m2 of leaves at crown bottom layer (g/m2 leaf)

    gpp_bottom(no) = gpp_bottom(no) + min(leaf_production_c, leaf_production_n)

    ! Current crown depth (m)

    crown_depth = (real(seib_height) - real(seib_bole)) * STEP

    ! Default of flag of purge crown bottom layer

    purge_flag = 0

    ! Day to purge crown bottom layer

    if (mod(day_count-2, crown_day) == (crown_day-1)) then

!write(*,*) 'day_count, crown_day', day_count, crown_day

       ! Case bottom layer production is not effective
       if ((gpp_bottom(no) / turn_cost) < crit_ratio) then

          purge_flag = 1

          if (crown_depth <= crown_depth_min) then

             purge_flag = 0
          end if

       ! Case bottom layer production is effective
       else
          purge_flag = 0
       end if

       ! When LAI is already reached at the maximum

       if (seib_lai > LA_max(p)) then

          purge_flag = 1

          if (crown_depth <= crown_depth_min) then

             purge_flag = 0
          end if
       end if

       ! Initialize gpp_bottom

       gpp_bottom(no) = 0.0d0

    end if

  end subroutine crown_bottom_purge

  !-----------------------------------------------------------------------
  function crown_allometry ( &
                                       ! *** Input ***
                                       !
                                       ! *** Input parameters ***
    crown_a                       , &  ! Tree crown allometric parameter
    crown_b                       , &  ! Tree crown allometric parameter
                                       !
                                       ! *** From SEIB-DGVM ***
    p                             , &  ! Species index (1: Rh, 2: Br)
    dbh_heartwood                 , &  ! Heartwood diameter (m)
    dbh_sapwood                     &  ! Sapwood diameter (m)
    ) result(ans)
    !
    ! !DESCRIPTION:
    ! Calculate potential crown diameter based on allometric relation.
    !
    ! !USES:
    !
    !
    ! !ARGUMENTS:
    implicit none
    real(8), intent(in) :: crown_a(:), crown_b(:)
    integer, intent(in) :: p
    real(8), intent(in) :: dbh_heartwood, dbh_sapwood
    !
    ! !LOCAL VARIABLES:
    real(8) :: dbh              ! DBH (m)
    real(8) :: ans              ! Potential crown diameter (m)
    !---------------------------------------------------------------------

    dbh = dbh_heartwood + dbh_sapwood
    ans = crown_a(p) * log(dbh) + crown_b(p)
    ans = max(ans, 0.4d0)

  end function crown_allometry

end module mod_crown_morphology
