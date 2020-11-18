module mod_energy_water_balance

  !-----------------------------------------------------------------------
  ! !DESCRIPTION:
  ! Compute leaf temperature and transpiration rate profile
  ! based on given stomatal conductance
  !
  ! !USES:
  !
  !
  ! !PUBLIC TYPES:
  implicit none
  !
  ! !PUBLIC MEMBER FUNCTIONS:
  public :: energy_water_balance
  !-----------------------------------------------------------------------

contains

  !-----------------------------------------------------------------------
  function energy_water_balance ( &
                                       ! *** Input ***
    flag_energy                   , &  ! Flag (1: for stomatal optimization, 0: other purposes)
    crown_area                    , &  ! Stand-level leaf crown area (m2 ground/tree)
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
    !
    ! !USES:
    use mod_param
    use mod_water_vapor
    !
    ! !ARGUMENTS:
    implicit none
    integer, intent(in) :: flag_energy
    real(8), intent(in) :: crown_area
    real(8), intent(in) :: gs_val
    type(universal_type), intent(in)    :: universal
    type(atmos_type),     intent(in)    :: atmos
    type(layer_type),     intent(in)    :: layer
    type(flux_type),      intent(inout) :: flux
    !
    ! !LOCAL VARIABLES:
    real(8) :: lambda                               ! Latent heat of vaporization (J/mol)
    real(8) :: esat                                 ! Saturation vapor pressure (Pa)
    real(8) :: desat                                ! Temperature derivative of saturation vapor pressure (Pa/K)
    real(8) :: term1, term2, term3, term4, term5    ! For computation of leaf temerature
    real(8) :: energy_check                         ! Energy balance check
    real(8) :: stand_et                             ! Total transpiration rate at stand-level (mol H2O/tree/s)
    real(8) :: val                                  ! Returned value of stand-level water balance error
    integer :: i                                    ! Layer index
    !---------------------------------------------------------------------

    associate ( &
                                         ! *** Input ***
    tfrz       => universal%tfrz    , &  ! Freezing point of water (K)
    n_layer    => universal%n_layer , &  ! Number of leaf layers
    tair       => atmos%tair        , &  ! Air temperature (K)
    pair       => atmos%pair        , &  ! Air pressure (Pa)
    cp         => atmos%cp          , &  ! Specific heat of air at constant pressure (J/mol/K)
    eair       => atmos%eair        , &  ! Vapor pressure in air (Pa)
    leaf_layer => layer%leaf_layer  , &  ! 1: Leaf, 0: No leaf
    dpai       => layer%dpai        , &  ! Layer leaf area index (m2 leaf/m2 ground)
    rn         => flux%rn           , &  ! Leaf net radiation profile (W/m2 leaf)
    gbh        => flux%gbh          , &  ! Leaf boundary layer conductance, heat (mol/m2 leaf/s)
    gbv        => flux%gbv          , &  ! Leaf boundary layer conductance, H2O (mol H2O/m2 leaf/s)
    sapmax     => flux%sapmax       , &  ! Stand-level maximum water uptake rate (mol H2O/s/tree)
                                         !
                                         ! *** Output ***
    tleaf      => flux%tleaf        , &  ! Leaf temperature (K)
    etflux     => flux%etflux         &  ! Leaf transpiration rate (mol H2O/m2 leaf/s)
    )

    ! Latent heat of vaporization (J/mol)

    lambda = 56780.3d0 - 42.84d0 * tair

    ! Saturation vapor pressure (Pa)

    call sat_vap (tfrz, tair, esat, desat)

    ! Computation for leaf layer

    do i = 1, n_layer
       if (leaf_layer(i) == 1) then ! Leaf layer

          ! Leaf temperature (K)

!          term1 = (rn(i) + 2*cp*gbh*tair) * (gbv + gs_val)
!          term2 = (lambda * (esat - desat*tair - eair) * gbv * gs_val) / pair
!          term3 = 2*cp*gbh*(gbv + gs_val)
!          term4 = (lambda * desat *gbv * gs_val) / pair
!          tleaf(i) = (term1 - term2) / (term3 + term4)
          term1 = rn(i) * pair * (gs_val + gbv(i))
          term2 = 2.0d0 * cp * gbh(i) * tair * pair * (gs_val + gbv(i))
          term3 = lambda * gs_val * gbv(i) * (esat - desat * tair - eair)
          term4 = 2.0d0 * cp * gbh(i) * pair * (gs_val + gbv(i))
          term5 = desat * lambda * gs_val * gbv(i)
          tleaf(i) = (term1 + term2 - term3) / (term4 + term5)

          ! Leaf transpiration rate (mol H2O/m2 leaf/s)

          etflux(i) = (esat + desat*(tleaf(i) - tair) - eair) * gs_val * gbv(i) &
                      / (pair * (gs_val + gbv(i)))

          ! Energy balance check

          energy_check = rn(i) - 2.0d0*cp*gbh(i)*(tleaf(i) - tair) - lambda*etflux(i)

          if (abs(energy_check) > 0.001d0) then
             write(*,*) 'Error: Energy_water_balance: Energy balance error', energy_check
          end if

       else ! No leaf layer

          tleaf(i) = 0.0d0
          etflux(i) = 0.0d0

       end if

    end do

    ! Total transpiration rate (mol H2O/tree/s)

    stand_et = 0.0d0
    do i = 1, n_layer
       stand_et = stand_et + etflux(i) * dpai * crown_area
    end do

    ! Water balance check

    val = sapmax - stand_et
    if (flag_energy == 1 .and. abs(val) > 0.001d0) then
       write(*,*) 'Error: Water balance error', val
       write(*,*) 'sapmax, stand_et, gs', sapmax, stand_et, gs_val
    end if

    end associate
  end function energy_water_balance

end module mod_energy_water_balance
