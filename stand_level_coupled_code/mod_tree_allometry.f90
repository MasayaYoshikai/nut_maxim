module mod_tree_allometry

  !-----------------------------------------------------------------------
  ! !DESCRIPTION:
  ! Calculate increment of D and H after stem biomass increment
  ! based on allometric equation.
  !
  ! Reference:
  ! Chave et al. (2005) Oecologia 145: 87-99
  ! Tree allometry and improved estimation of carbon stocks and balance
  ! in tropical forests
  !
  ! The allometric model is expresses as:
  !
  ! M = 0.051 * rho * D^2 * H
  ! where
  ! M: above-ground biomass (kg)
  ! rho: wood density (g/cm3)
  ! D: stem diameter (cm)
  ! H: tree height (m)
  !
  ! The values of rho were from Global Wood Density Database (Zanne et al., 2009)
  ! The mean value of wood density for the focal species was applied.
  ! Rhizophora stylosa: 0.840 g/cm3
  ! Bruguiera gymnorhiza: 0.764 g/cm3
  !
  ! !USES:
  !
  !
  ! !PUBLIC TYPES:
  implicit none
  !
  ! !PUBLIC MEMBER FUNCTIONS:
  public :: tree_allometry
  public :: height_allometry
  public :: proot_allometry
  !-----------------------------------------------------------------------

contains

  !-----------------------------------------------------------------------
  subroutine tree_allometry ( &
                                       ! *** Input ***
    dbh_heartwood                 , &  ! Heartwood diameter (m)
    dbh_sapwood                   , &  ! Sapwood diameter (m)
    tree_h                        , &  ! Tree height (m)
    trunk_biomass                 , &  ! Stand-level trunk biomass (g/tree)
    wood_rho                      , &  ! Wood density (g/cm3)
    ALM5                          , &  ! Sapwood diameter proportion (m sapwood/m dbh)
    biomass_increment             , &  ! Increment of trunk biomass (g/tree)
                                       !
                                       ! *** Output ***
    new_dbh_heartwood             , &  ! New heartwood diameter after biomass increment (m)
    new_dbh_sapwood               , &  ! New sapwood diameter after biomass increment (m)
    new_tree_h                      &  ! New tree height after biomass increment (m)
    )
    !
    ! !DESCRIPTION:
    !
    ! !USES:
    !
    !
    ! !ARGUMENTS:
    implicit none
    real(8), intent(in)  :: dbh_heartwood, dbh_sapwood, tree_h
    real(8), intent(in)  :: trunk_biomass, wood_rho, ALM5, biomass_increment
    real(8), intent(out) :: new_dbh_heartwood, new_dbh_sapwood, new_tree_h
    !
    ! !LOCAL VARIABLES:
    real(8) :: dbh                          ! Current DBH (cm)
    real(8) :: new_dbh                      ! New DBH after biomass increment (m)
    real(8) :: new_trunk_biomass            ! New trunk biomass after biomass increment (kg/tree)
    !---------------------------------------------------------------------

    ! Above-ground biomass after biomass increment (kg/tree)
    ! !!! Caution1: the unit is kg/tree !!!
    ! !!! Caution2: trunk_biomass is supposed to include below-ground woody roots.
    !               Currently, it is assumed to be above-ground biomass.

    new_trunk_biomass = (trunk_biomass + biomass_increment) / 1000.0d0

    ! Current DBH (cm)
    ! !!! Caution: the unit is cm !!!

    dbh = (dbh_heartwood + dbh_sapwood) * 100.0d0

    ! New DBH after biomass increment (m)
    ! !!! Caution: the unit is m !!!

    new_dbh = (sqrt(new_trunk_biomass / (0.051d0 * wood_rho * tree_h))) * 0.01d0
    new_dbh_heartwood = (1.0d0 - ALM5) * new_dbh
    new_dbh_sapwood = ALM5 * new_dbh

    ! New tree height after biomass increment (m)

    new_tree_h = new_trunk_biomass / (0.051d0 * wood_rho * dbh**2.0d0)

    ! Error check
!    write(*,*) 'biomass_increment', biomass_increment
!    write(*,*) 'trunk_biomass', trunk_biomass
!    write(*,*) 'new_trunk_biomass', new_trunk_biomass
!    write(*,*) 'dbh', dbh
!    write(*,*) 'tree_h', tree_h
!    write(*,*) 'wood_rho', wood_rho
!    write(*,*) 'new_dbh', new_dbh
!    write(*,*) 'new_tree_h', new_tree_h
    if (new_dbh < dbh/100.0d0) then
       write(*,*) 'Error: Negative D increment'
       write(*,*) 'biomass_inc, dbh, new_dbh', biomass_increment, dbh, new_dbh
    end if
    if (new_tree_h < tree_h) then
       write(*,*) 'Error: Negative H increment'
       write(*,*) 'biomass_inc, tree_h, new_tree_h', biomass_increment, tree_h, new_tree_h
    end if

  end subroutine tree_allometry

  !-----------------------------------------------------------------------
  function height_allometry ( &
                                       ! *** Input ***
                                       !
                                       ! *** Input parameters ***
    tree_h_a                      , &  ! Tree height allometric parameter
    tree_h_b                      , &  ! Tree height allometric parameter
    HGT_max                       , &  ! Maximum tree height (m)
                                       !
                                       ! *** From SEIB-DGVM ***
    p                             , &  ! Species index (1: Rh, 2: Br)
    dbh_heartwood                 , &  ! Heartwood diameter (m)
    dbh_sapwood                     &  ! Sapwood diameter (m)
    ) result(ans)
    !
    ! !DESCRIPTION:
    ! Calculate potential tree height based on allometric relation.
    !
    ! !USES:
    !
    !
    ! !ARGUMENTS:
    implicit none
    real(8), intent(in) :: tree_h_a(:), tree_h_b(:), HGT_max(:)
    integer, intent(in) :: p
    real(8), intent(in) :: dbh_heartwood, dbh_sapwood
    !
    ! !LOCAL VARIABLES:
    real(8) :: dbh              ! DBH (m)
    real(8) :: ans              ! Potential tree height (m)
    !---------------------------------------------------------------------

    dbh = dbh_heartwood + dbh_sapwood
    ans = tree_h_a(p) * (dbh ** tree_h_b(p))
    ans = max(ans, 1.80d0)
    ans = min(ans, HGT_max(p))

  end function height_allometry

  !-----------------------------------------------------------------------
  subroutine proot_allometry ( &
                                       ! *** Input ***
                                       !
                                       ! *** Input parameters ***
    wood_rho                      , &  ! Wood density (g/cm3)
    pr_s_a                        , &  ! Slope for the scaling factor of prop root system
    pr_s_b                        , &  ! Intercept for the scaling factor of prop root system
    pr_h_a                        , &  ! Slope for the maximum root height of prop root system
    pr_h_b                        , &  ! Intercept for the maximum root height of prop root system
    pr_d                          , &  ! Mean prop root diameter (m)
                                       !
                                       ! *** From SEIB-DGVM ***
    p                             , &  ! Species index (1: Rh, 2: Br)
    dbh_heartwood                 , &  ! Heartwood diameter (m)
    dbh_sapwood                   , &  ! Sapwood diameter (m)
    mass_above_root               , &  ! Stand-level above-ground root biomass (g/tree)
                                       !
                                       ! *** In/Output ***
    d_trunk                       , &  ! dTB/dt: Daily trunk biomass change (g trunk/tree/day)
                                       !
                                       ! *** Output ***
    d_above_root                    &  ! dAGR/dt: Daily above-ground root biomass change (g/tree/day)
    )
    !
    ! !DESCRIPTION:
    ! Calculate needed prop root biomass corresponding to the DBH
    ! based on allometric relationship.
    !
    ! !USES:
    use data_structure, only : PI
    !
    !
    ! !ARGUMENTS:
    implicit none
    real(8), intent(in)    :: wood_rho(:), pr_s_a(:), pr_s_b(:)
    real(8), intent(in)    :: pr_h_a(:), pr_h_b(:), pr_d(:)
    integer, intent(in)    :: p
    real(8), intent(in)    :: dbh_heartwood, dbh_sapwood, mass_above_root
    real(8), intent(inout) :: d_trunk
    real(8), intent(out)   :: d_above_root
    !
    ! !LOCAL VARIABLES:
    real(8) :: dbh                          ! Current DBH (m)
    real(8) :: s_factor                     ! Scaling factor for prop root system
    real(8) :: hr_max                       ! Maximum prop root height (cm)
    real(8) :: hr                           ! Prop root height (cm)
    integer :: z                            ! Layer index
    integer :: r                            ! Root index
    real(8), parameter :: hr_min = 5.0d0    ! Minimum prop root height (cm)
    integer, parameter :: n_layer = 20      ! Number of vertical layers
    real(8), parameter :: dz = 10.0d0       ! Layer thickness (cm)
    real(8), parameter :: l_factor = 14.5d0 ! Factor for converting prop root number to total length in cm
    real(8), dimension(n_layer) :: h_layer  ! Height of upper boundary of each vertical layer (cm)
    real(8), dimension(n_layer) :: n_root   ! Number of prop roots in each vertical layer
    real(8), dimension(n_layer) :: v_root   ! Total volume of prop roots in each vertical layer (cm3)
    real(8) :: total_v_root                 ! Total volume of prop roots (cm3)
    real(8) :: proot_biomass                ! Prop root biomass needed based on allometric relationship (g/tree)
    !---------------------------------------------------------------------

    if (pr_s_a(p) < 0.0d0) then ! For Rhizophora genus

       ! Current DBH (m)

       dbh = dbh_heartwood + dbh_sapwood

       ! Scaling factor for prop root system

       s_factor = 1.0d0 - 10.0d0 **(pr_s_a(p)*log10(dbh)+pr_s_b(p))
       s_factor = max(s_factor, 0.3d0)  ! Prevent very low values

       ! Maximum prop root height (cm)

       hr_max = (pr_h_a(p) * dbh + pr_h_b(p)) * 100.0d0

       ! Height of upper boundary of each vertical layer (cm)

       h_layer(:) = 0.0d0 ! Zero out

       do z = 1, n_layer

          h_layer(z) = real(n_layer - z + 1) * dz

       end do

       !---------------------------------------------------------------------
       ! Computation of number of prop roots in each vertical layer
       !---------------------------------------------------------------------

       ! Initialize prop root height (cm)

       hr = hr_max

       ! Zero out

       n_root(:) = 0.0d0
       v_root(:) = 0.0d0

       ! Loop from the highest to the lowest prop root

       do r = 1, 1000

          ! Loop for the vertical layer

          do z = 1, n_layer

             if (hr > h_layer(z)) then
                n_root(z) = n_root(z) + 1.0d0
             else
                if (z == n_layer) then
                   n_root(z) = n_root(z) + hr/dz
                elseif (hr > h_layer(z+1)) then
                   n_root(z) = n_root(z) + (hr - h_layer(z+1))/dz
                end if
             end if

          end do

          ! Next prop root height

          hr = hr * s_factor

          ! Exit the loop when the prop root height is lower than minimum value.

          if (hr < hr_min) exit

       end do

       !---------------------------------------------------------------------
       ! Total volume (cm3) and biomass (g) of prop roots
       !---------------------------------------------------------------------

       ! Total volume of prop roots in each vertical layer (cm3)

       do z = 1, n_layer

          v_root(z) = l_factor * n_root(z) * PI * (pr_d(p)*100.0d0 * 0.5d0)**2.0d0
!                     --------------------   -------------------------------------
!                       Total length (cm)         Cross-sectional area (cm2)
       end do

       ! Integrate for the vertical layers

       total_v_root = sum(v_root)  ! (cm3)

       ! Prop root biomass needed based on allometric relationship (g/tree)

       proot_biomass = total_v_root * wood_rho(p)

       ! When the current prop root biomass is lacking with respect to the current DBH

       if (mass_above_root < proot_biomass) then

          if (d_trunk >= (proot_biomass - mass_above_root)) then

             d_above_root = proot_biomass - mass_above_root
             d_trunk = d_trunk - d_above_root

          else

             d_above_root = 0.50d0 * d_trunk
             d_trunk = d_trunk - d_above_root

          end if

       ! When the current prop root biomass is sufficient

       else

          d_above_root = 0.0d0

       end if

    else ! For other genus

       d_above_root = 0.0d0

    end if

  end subroutine proot_allometry

end module mod_tree_allometry
