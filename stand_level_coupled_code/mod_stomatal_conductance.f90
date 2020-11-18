module mod_stomatal_conductance

  !-----------------------------------------------------------------------
  ! !DESCRIPTION:
  !
  ! !USES:
  !
  !
  ! !PUBLIC TYPES:
  implicit none
  !
  ! !PUBLIC MEMBER FUNCTIONS:
  public  :: stomatal_conductance1
  public  :: stomatal_conductance2
  !
  ! !PRIVATE MEMBER FUNCTIONS:
  private :: stomata_efficiency        ! Water-use efficiency for maximum gs
  !-----------------------------------------------------------------------

contains

  !-----------------------------------------------------------------------
  subroutine stomatal_conductance1 ( &
                                       ! *** Input ***
    crown_area                    , &  ! Stand-level leaf crown area (m2 ground/tree)
    gs_max                        , &  ! Maximum stomatal conductance (mol H2O/m2/s)
                                       !
                                       ! *** Input type ***
    universal                     , &  ! Universal variables
    atmos                         , &  ! Atmospheric variables
    layer                         , &  ! Layer variables
                                       !
                                       ! *** Input/Output ***
    flux                            &  ! Flux variables
    )
    !
    ! !DESCRIPTION:
    ! Optimizing stomatal conductance for the maximum water uptake rate,
    ! then compute leaf temperature and transpiration profile using
    ! the optimized gs.
    !
    ! !USES:
    use mod_param
    use mod_energy_water_balance
    use mod_math_tools
    !
    ! !ARGUMENTS:
    implicit none
    real(8), intent(in) :: crown_area
    real(8), intent(in) :: gs_max
    type(universal_type), intent(in)    :: universal
    type(atmos_type),     intent(in)    :: atmos
    type(layer_type),     intent(in)    :: layer
    type(flux_type),      intent(inout) :: flux
    !
    ! !LOCAL VARIABLES:
    integer :: flag_energy                 ! Flag (1: for stomatal optimization, 0: other purposes)
    real(8) :: gs_min                      ! Minimum gs (mol H2O/m2/s)
    real(8) :: check1, check2              ! For brent's method
    real(8), parameter :: tol = 0.0001d0   ! Accuracy for water balance (mol H2O/tree/s)
    real(8) :: val                         ! Dummy
    !---------------------------------------------------------------------

    associate ( &
                                 ! *** Output ***
    gs     => flux%gs         &  ! Stand-level stomatal conductance (mol H2O/m2/s)
    )

    ! Minimum gs (mol H2O/m2/s)

    gs_min = 0.002d0

    ! Check if there is a optimal gs value for realizing
    ! the maximum sap flow rate within gs1 and gs2

    flag_energy = 0
    check1 = energy_water_balance (flag_energy, crown_area, gs_min, &
                                   universal, atmos, layer, flux)
    check2 = energy_water_balance (flag_energy, crown_area, gs_max, &
                                   universal, atmos, layer, flux)

    if (check1 * check2 < 0.0d0) then

       ! Calculate gs using the function Energy_water_balance to iterate gs
       ! to an accuracy of tol (mol H2O/tree/s)

       gs = zbrent1 ('Stomatal optimaization1', energy_water_balance, flag_energy, &
                    crown_area, universal, atmos, layer, flux, gs_min, gs_max, tol)

    else if (check2 > 0.0d0) then

       write(*,*) 'Error: Maximum gs in stomatal optimization1'
       gs = gs_max

    else

       ! Stomata need to close

       gs = gs_min

    end if

    ! Compute leaf temperature and transpiration profile using the optimized gs

    flag_energy = 1
    val = energy_water_balance (flag_energy, crown_area, gs, &
                                universal, atmos, layer, flux)

    end associate
  end subroutine stomatal_conductance1

  !-----------------------------------------------------------------------
  subroutine stomatal_conductance2 ( &
                                       ! *** Input ***
    o2air                         , &  ! Atmospheric O2 (mmol/mol)
    co2air                        , &  ! Atmospheric CO2 (umol/mol)
    iota                          , &  ! Stomatal water-use efficiency (umol CO2/ mol H2O)
    crown_area                    , &  ! Stand-level leaf crown area (m2 ground/tree)
    n_leaf_layer                  , &  ! Number of leaf layers
                                       !
                                       ! *** Input type ***
    universal                     , &  ! Universal variables
    atmos                         , &  ! Atmospheric variables
    layer                         , &  ! Layer variables
                                       !
                                       ! *** Input/Output ***
    flux                            &  ! Flux variables
    )
    !
    ! !DESCRIPTION:
    ! Photosynthesis and stomatal conductance with optimization of water use efficiency.
    !
    ! !USES:
    use mod_param
    use mod_math_tools
    use mod_energy_water_balance
    use mod_photosynthesis
    !
    ! !ARGUMENTS:
    implicit none
    real(8), intent(in) :: o2air, co2air
    real(8), intent(in) :: iota
    real(8), intent(in) :: crown_area
    integer, intent(in) :: n_leaf_layer
    type(universal_type), intent(in)    :: universal
    type(atmos_type),     intent(in)    :: atmos
    type(layer_type),     intent(in)    :: layer
    type(flux_type),      intent(inout) :: flux
    !
    ! !LOCAL VARIABLES:
    real(8) :: gs1, gs2                    ! Initial guess for gs (mol H2O/m2/s)
    real(8) :: check1, check2              ! For brent's method
    real(8), parameter :: tol = 0.004d0    ! gs is updated to accuracy tol (mol H2O/m2/s)
    integer :: flag_energy                 ! Flag (1: for stomatal optimization, 0: other purposes)
    real(8) :: val                         ! Dummy
    !---------------------------------------------------------------------

    associate ( &
                                 ! *** Output ***
    gs     => flux%gs         &  ! Stand-level stomatal conductance (mol H2O/m2/s)
    )

    ! Low and high initial estimates for gs (mol H2O/m2/s)

    gs1 = 0.002d0
    gs2 = gs                     ! Maximum gs based on water uptake maximization.

    ! Check for minimum stomatal conductance linked to low light
    ! based on the water-use efficiency for gs1 and gs2

    check1 = stomata_efficiency (    &
                                       ! *** Input ***
    o2air                         , &  ! Atmospheric O2 (mmol/mol)
    co2air                        , &  ! Atmospheric CO2 (umol/mol)
    iota                          , &  ! Stomatal water-use efficiency (umol CO2/ mol H2O)
    crown_area                    , &  ! Stand-level leaf crown area (m2 ground/tree)
    n_leaf_layer                  , &  ! Number of leaf layers
    gs1                           , &  ! Value for gs to use in calculations
                                       !
                                       ! *** Input type ***
    universal                     , &  ! Universal variables
    atmos                         , &  ! Atmospheric variables
    layer                         , &  ! Layer variables
                                       !
                                       ! *** Input/Output ***
    flux                            &  ! Flux variables
    )

    check2 = stomata_efficiency (    &
                                       ! *** Input ***
    o2air                         , &  ! Atmospheric O2 (mmol/mol)
    co2air                        , &  ! Atmospheric CO2 (umol/mol)
    iota                          , &  ! Stomatal water-use efficiency (umol CO2/ mol H2O)
    crown_area                    , &  ! Stand-level leaf crown area (m2 ground/tree)
    n_leaf_layer                  , &  ! Number of leaf layers
    gs2                           , &  ! Value for gs to use in calculations
                                       !
                                       ! *** Input type ***
    universal                     , &  ! Universal variables
    atmos                         , &  ! Atmospheric variables
    layer                         , &  ! Layer variables
                                       !
                                       ! *** Input/Output ***
    flux                            &  ! Flux variables
    )

    if (check1 * check2 < 0.0d0) then

       ! Calculate gs using the function stomata_efficiency to iterate gs
       ! to an accuracy of tol (mol H2O/m2/s)

       gs = zbrent2 ('Stomatal optimaization2', stomata_efficiency, o2air, co2air, iota, &
                    crown_area, n_leaf_layer, universal, atmos, layer, flux, gs1, gs2, tol)

    else

       ! Low light. Set gs to minimum conductance

       gs = 0.002d0

    end if

    ! Compute leaf temperature and transpiration profile,
    ! and photosynthetic rate for this gs

    flag_energy = 0
    val = energy_water_balance (flag_energy, crown_area, gs, &
                                universal, atmos, layer, flux)

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

    end associate
  end subroutine stomatal_conductance2

  !-----------------------------------------------------------------------
  function stomata_efficiency ( &
                                       ! *** Input ***
    o2air                         , &  ! Atmospheric O2 (mmol/mol)
    co2air                        , &  ! Atmospheric CO2 (umol/mol)
    iota                          , &  ! Stomatal water-use efficiency (umol CO2/ mol H2O)
    crown_area                    , &  ! Stand-level leaf crown area (m2 ground/tree)
    n_leaf_layer                  , &  ! Number of leaf layers
    gs_val                        , &  ! Value for gs to use in calculations
                                       !
                                       ! *** Input type ***
    universal                     , &  ! Universal variables
    atmos                         , &  ! Atmospheric variables
    layer                         , &  ! Layer variables
                                       !
                                       ! *** Input/Output ***
    flux                            &  ! Flux variables
    ) result(val)
    !
    ! !DESCRIPTION:
    ! Stomata water-use efficiency check to determine maximum gs.
    ! For the stomatal conductance gs_val, calculate photosynthesis
    ! for an increase in stomatal conductance equal to "delta".
    ! The returned value is positive if this increase produces a change in
    ! photosynthesis > iota*vpd*delta.
    ! The returned value is negative if the increase produces a change in
    ! photosynthesis < iota*vpd*delta.
    !
    ! !USES:
    use mod_param
    use mod_energy_water_balance
    use mod_photosynthesis
    !
    ! !ARGUMENTS:
    implicit none
    real(8), intent(in) :: o2air, co2air
    real(8), intent(in) :: iota
    real(8), intent(in) :: crown_area
    integer, intent(in) :: n_leaf_layer
    real(8), intent(in) :: gs_val
    type(universal_type), intent(in)    :: universal
    type(atmos_type),     intent(in)    :: atmos
    type(layer_type),     intent(in)    :: layer
    type(flux_type),      intent(inout) :: flux
    !
    ! !LOCAL VARIABLES:
    integer :: flag_energy           ! Flag (1: for stomatal optimization, 0: other purposes)
    integer :: i                     ! Layer index
    real(8) :: val_dummy             ! Dummy
    real(8) :: delta                 ! Small difference for gs (mol H2O/m2/s)
    real(8) :: gs2                   ! Lower value for gs (mol H2O/m2/s)
    real(8) :: an_mean2              ! Canopy mean leaf net photosynthesis at gs2 (umol CO2/m2 leaf/s)
    real(8) :: gs1                   ! Higher value for gs (mol H2O/m2/s)
    real(8) :: an_mean1              ! Canopy mean leaf net photosynthesis at gs1 (umol CO2/m2 leaf/s)
    real(8) :: wue                   ! Water-use efficiency check
    real(8) :: val                   ! Returned value
    !---------------------------------------------------------------------

    associate ( &
    n_layer    => universal%n_layer  , &  ! Number of layers (-) (Max_hgt in SEIB-DGVM)
    leaf_layer => layer%leaf_layer   , &  ! 1: Leaf, 0: No leaf
    gs         => flux%gs            , &  ! Stand-level stomatal conductance (mol H2O/m2/s)
    an         => flux%an              &  ! Leaf net photosynthesis (umol CO2/m2 leaf/s)
    )

    ! Specify "delta" as a small difference in gs (mol H2O/m2/s)

    delta = 0.001d0

    ! --- Flux calculation with lower gs (gs_val - delta)

    gs2 = gs_val - delta
    gs  = gs2

    ! Leaf temperature profile

    flag_energy = 0
    val_dummy = energy_water_balance (flag_energy, crown_area, gs2, &
                                      universal, atmos, layer, flux)

    ! Leaf photosynthesis

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

    ! Canopy mean leaf net photosynthesis (umol CO2/m2 leaf/s)

    an_mean2 = 0.0d0
    do i = 1, n_layer
       if (leaf_layer(i) == 1) then ! Leaf layer
          an_mean2 = an_mean2 + an(i)
       end if
    end do
    an_mean2 = an_mean2 / real(n_leaf_layer)

    ! --- Flux calculation with higher gs (gs_val)

    gs1 = gs_val
    gs  = gs1

    ! Leaf temperature profile

    flag_energy = 0
    val_dummy = energy_water_balance (flag_energy, crown_area, gs1, &
                                      universal, atmos, layer, flux)

    ! Leaf photosynthesis

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

    ! Canopy mean leaf net photosynthesis (umol CO2/m2 leaf/s)

    an_mean1 = 0.0d0
    do i = 1, n_layer
       if (leaf_layer(i) == 1) then ! Leaf layer
          an_mean1 = an_mean1 + an(i)
       end if
    end do
    an_mean1 = an_mean1 / real(n_leaf_layer)

    ! Efficiency check: wue < 0 when d(An) / d(gs) < iota

    wue = (an_mean1 - an_mean2) - iota * delta

    ! Return value

    val = wue

    end associate
  end function stomata_efficiency

end module mod_stomatal_conductance
