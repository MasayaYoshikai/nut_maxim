module mod_radiation

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
  public :: extinction_coefficient
  public :: longwave_rad
  !-----------------------------------------------------------------------

contains

  !-----------------------------------------------------------------------
  subroutine extinction_coefficient ( &
                                       ! *** Input ***
    n_spe                         , &  ! Number of mangrove species
    rhol_vis                      , &  ! Leaf relflectance to VIS (-)
    rhol_nir                      , &  ! Leaf relflectance to NIR (-)
    taul_vis                      , &  ! Leaf transmittance to VIS (-)
    taul_nir                      , &  ! Leaf transmittance to NIR (-)
    xl                            , &  ! Departure of leaf angle from spherical orientation (-)
                                       !
                                       ! *** From SEIB-DGVM ***
    Max_loc                       , &  ! Dimension of the virtual forest (m)
    sl_hgt                        , &  ! Solar elevation angle at midday (degree)
    tree_exist                    , &  ! Flag of tree presence
    la                            , &  ! Leaf area per tree (m2 leaf/tree)
    height                        , &  ! SEIB-DGVM specific tree height (the unit is STEP!!)
    bole                          , &  ! SEIB-DGVM specific bole height (the unit is STEP!!)
    pft                           , &  ! Species index (1: Rh, 2: Br)
                                       !
                                       ! *** Output ***
    dpai_layer                    , &  ! Layer leaf area index for each PFT for each forest layer (m2 leaf/me ground)
    dpai_layer_sum                , &  ! Sum of layer leaf area index for each forest layer (m2 leaf/me ground)
    kbm_vis                       , &  ! Scattering adjusted light extinction coefficient for VIS (-)
    kbm_nir                       , &  ! Scattering adjusted light extinction coefficient for NIR (-)
    albvegb_vis                   , &  ! (Direct beam) vegetation albedo for VIS, non-horizontal leaves
    albvegb_nir                   , &  ! (Direct beam) vegetation albedo for NIR, non-horizontal leaves
    nbot                          , &  ! Index for bottom leaf layer
    ntop                          , &  ! Index for top leaf layer
    td                              &  ! Exponential transmittance of diffuse radiation through a single leaf layer
    )
    !
    ! !DESCRIPTION:
    ! Scattering adjusted light extinction coefficient for VIS/NIR, and
    ! exponential transmittance of diffuse radiation through a single leaf layer
    ! which is needed for calculating longwave radiation distribution.
    !
    ! !USES:
    use data_structure, only : Max_no, Max_hgt, PFT_no
    !
    ! !ARGUMENTS:
    implicit none
    integer, intent(in)  :: n_spe
    real, intent(in)     :: rhol_vis(:), rhol_nir(:), taul_vis(:), taul_nir(:), xl(:)
    integer, intent(in)  :: Max_loc
    real, intent(in)     :: sl_hgt
    logical, intent(in)  :: tree_exist(:)
    real(8), intent(in)  :: la(:)
    integer, intent(in)  :: height(:), bole(:), pft(:)
    real(8), intent(out) :: dpai_layer(:,:), dpai_layer_sum(:)
    real, intent(out)    :: kbm_vis(:), kbm_nir(:), albvegb_vis, albvegb_nir
    integer, intent(out) :: nbot, ntop
    real(8), intent(out) :: td(:)
    !
    ! !LOCAL VARIABLES:
    real    :: solar_zen                     ! Solar zenith angle (radians)
    integer :: no                            ! Tree index
    real(8) :: dpai_tree                     ! Layer leaf area index for a tree (m2 leaf/m2 ground)
    integer :: p                             ! Species index
    integer :: i                             ! Layer index
    real, dimension(PFT_no) :: omega_vis     ! Leaf scattering coefficient for VIS
    real, dimension(PFT_no) :: omega_nir     ! Leaf scattering coefficient for NIR
    real, dimension(PFT_no) :: chil          ! Departure of leaf angle from spherical orientation (-0.4 <= xl <= 0.6)
    real, dimension(PFT_no) :: phi1          ! Term in Ross-Goudriaan function for gdir
    real, dimension(PFT_no) :: phi2          ! Term in Ross-Goudriaan function for gdir
    real    :: gdir                          ! Relative projected area of leaf elements in the direction of solar beam
    real    :: kb                            ! Direct beam extinction coefficient
    real    :: albvegh_vis                   ! Vegetation albedo for VIS, horizontal leaves
    real    :: albvegh_nir                   ! Vegetation albedo for NIR, horizontal leaves
    integer :: canopy_index                  ! Canopy index
    integer :: j                             ! Sky angle index
!    real(8) :: term1                         ! Term in exponential function
!    real(8) :: angle                         ! Sky angles (5, 15, 25, 35, 45, 55, 65, 75, 85 degrees)
!    real(8), dimension(PFT_no) :: gdirj      ! Relative projected area of leaf elements in the direction of sky angle
!    real(8), parameter :: pi = 3.14159265d0
    !---------------------------------------------------------------------

    !---------------------------------------------------------------------
    ! Solar zenith angle (radians)
    !---------------------------------------------------------------------

    ! Note: The effect of solar zenith angle on light attenuation was
    !       removed as in SEIB-DGVM.

!    solar_zen = (90.0d0 - sl_hgt) * (pi / 180.0d0)
    solar_zen = 0.0d0

    !---------------------------------------------------------------------
    ! Layer leaf area index for each PFT for each forest layer (m2 leaf/me ground)
    !---------------------------------------------------------------------

    ! Zero out for all layers.

    dpai_layer(:,:) = 0.0d0

    ! Summing every tree per pft type

    do no = 1, Max_no

       if ( .not. tree_exist(no) ) cycle

       ! Layer leaf area index for a tree (m2 leaf/m2 ground)

       dpai_tree = la(no) / real(height(no) - bole(no))

       do i = bole(no)+1, height(no)

          dpai_layer(pft(no),i) = dpai_layer(pft(no),i) + dpai_tree

       end do
    end do

    ! Deviding by plot area

    dpai_layer(:,:) = dpai_layer(:,:) / real(Max_loc) / real(Max_loc)

    !------------------------------------------------------------------
    ! Sum of layer leaf area index for each forest layer (m2 leaf/me ground)
    !------------------------------------------------------------------

    do i = 1, Max_hgt

       dpai_layer_sum(i) = 0.0d0
       do p = 1, n_spe
          dpai_layer_sum(i) = dpai_layer_sum(i) + dpai_layer(p,i)
       end do
    end do

    !---------------------------------------------------------------------
    ! Leaf scattering coefficient (-)
    !---------------------------------------------------------------------

    do p = 1, n_spe

       omega_vis(p) = rhol_vis(p) + taul_vis(p)
       omega_nir(p) = rhol_nir(p) + taul_nir(p)

    end do

    !---------------------------------------------------------------------
    ! Light extinction coefficient
    !---------------------------------------------------------------------

    do p = 1, n_spe

       chil(p) = min(max(xl(p), -0.4d0), 0.6d0)
       if (abs(chil(p)) <= 0.01d0) chil(p) = 0.01d0

       phi1(p) = 0.5d0 - 0.633d0*chil(p) - 0.330d0*chil(p)*chil(p)
       phi2(p) = 0.877d0 * (1.0d0 - 2.0d0*phi1(p))

       gdir = phi1(p) + phi2(p) * cos(solar_zen)
       kb = gdir / cos(solar_zen)
       kb = min(kb, 40.0d0)

       ! Adjust for scattering

       kbm_vis(p) = kb * sqrt(1.0d0 - omega_vis(p))
       kbm_nir(p) = kb * sqrt(1.0d0 - omega_nir(p))

    end do

    !---------------------------------------------------------------------
    ! Vegetation albedo
    !---------------------------------------------------------------------

    ! Vegetation albedo, horizontal leaves
    ! Note: This equation is valid only for canopy with the kb (same leaf angle xl).

    albvegh_vis = (1.0d0 - sqrt(1.0d0 - omega_vis(1))) / (1.0d0 + sqrt(1.0d0 - omega_vis(1)))
    albvegh_nir = (1.0d0 - sqrt(1.0d0 - omega_nir(1))) / (1.0d0 + sqrt(1.0d0 - omega_nir(1)))

    ! (Direct beam) vegetation albedo, non-horizontal leaves
    ! Note: This equation is valid only for canopy with the kb (same leaf angle xl).

!   albvegb = 2.0d0 * kb / (kb + kd) * albvegh
    albvegb_vis = albvegh_vis
    albvegb_nir = albvegh_nir

    !---------------------------------------------------------------------
    ! Identifying top and bottom layer
    !---------------------------------------------------------------------

    ! Identifying canopy bottom layer

    canopy_index = 0
    do i = 1, Max_hgt
       if (dpai_layer_sum(i) > 0.0d0 .and. canopy_index == 0) then
          nbot = i
          canopy_index = 1
       end if
    end do

    ! Identifying canopy top layer

    canopy_index = 0
    do i = Max_hgt, 1, -1
       if (dpai_layer_sum(i) > 0.0d0 .and. canopy_index == 0) then
          ntop = i
          canopy_index = 1
       end if
    end do

    !---------------------------------------------------------------------
    ! Diffuse transmittance for a single layer estimated for nine sky
    ! angles in increments of 10 degrees (also needed for longwave radiation)
    !---------------------------------------------------------------------

    ! Zero out for all layers

    do i = 1, Max_hgt
       td(i) = 0.0d0
    end do

    ! Leaf layers

!!! ------------------- Based on multiple angle for diffused radiation ------------------- !!!
!    if (nbot > 0) then
!       do i = nbot, ntop
!          do j = 1,9
!             angle = (5.0d0 + (j - 1) * 10.0d0) * pi / 180.0d0
!
!             term1 = 0.0d0
!             do p = 1, n_spe
!                gdirj(p) = phi1(p) + phi2(p) * cos(angle)
!                ! term1 = -gdirj(p) / cos(angle) * dpai_layer(p,i)
!                term1 = term1 + (-gdirj(p) / cos(angle) * dpai_layer(p,i))
!             end do
!             td(i) = td(i) + exp(term1) * sin(angle) * cos(angle)
!          end do
!          td(i) = td(i) * 2.0d0 * (10.0d0 * pi / 180.0d0)
!
!          ! Prevent td value higher than 1 which will cause error in
!          ! radiative transfer model
!          if (td(i) > 0.999d0) then
!             td(i) = 0.999d0
!             write(*,*) 'Caution: dpai is too low. td value set to 0.999.'
!          end if
!       end do
!    end if
!!! ------------------- Based on direct beam extinction coefficient ------------------- !!!

    ! Note: This may be valid because SEIB-DGVM assumed the diffused radiation only from
    !       the vertical direction.

    if (nbot > 0) then
       do i = nbot, ntop

          ! Note: This equation is valid only for canopy with the kb (same leaf angle xl).
          td(i) = exp(-kb * dpai_layer_sum(i))

       end do
    end if

  end subroutine extinction_coefficient

    !-----------------------------------------------------------------------
    subroutine longwave_rad ( &
                                         ! *** Input ***
      n_spe                         , &  ! Number of mangrove species
      emleaf                        , &  ! Leaf emissivity (-)
      dpai_layer_sum                , &  ! Sum of layer leaf area index for each forest layer (m2 leaf/me ground)
      nbot                          , &  ! Index for bottom leaf layer
      ntop                          , &  ! Index for top leaf layer
      tg                            , &  ! Soil surface temperature (K)
      irsky                         , &  ! Downward atmospheric longwave radiation at time step of computation (W/m2)
      tleaf_plot                    , &  ! Plot-scale leaf temperature (K)
      td                            , &  ! Exponential transmittance of diffuse radiation through a single leaf layer
                                         !
                                         ! *** Output ***
      irleaf                        , &  ! Leaf absorbed longwave radiation for canopy layer (W/m2 leaf)
      irveg                         , &  ! Absorbed longwave radiation, vegetation (W/m2)
      ircan                         , &  ! Upward longwave radiation above canopy (W/m2)
      irsoi                           &  ! Absorbed longwave radiation, ground (W/m2)
      )
      !
      ! !DESCRIPTION:
      ! Longwave radiation transfer through canopy using Norman (1979)
      !
      ! !USES:
      use data_structure, only : Max_hgt
      use mod_math_tools, only : tridiag
      !
      ! !ARGUMENTS:
      implicit none
      integer, intent(in)  :: n_spe
      real(8), intent(in)  :: emleaf, dpai_layer_sum(:)
      integer, intent(in)  :: nbot, ntop
      real(8), intent(in)  :: tg, irsky, tleaf_plot, td(:)
      real(8), intent(out) :: irleaf(:), irveg, ircan, irsoi
      !
      ! !LOCAL VARIABLES:
!      integer :: p                              ! Species index
      integer :: ic                             ! Layer index
      integer :: icm1                           ! Layer below ic (ic-1)
      real(8), parameter :: sb = 5.67e-8        ! Stefan-Boltzmann constant (W/m2/K4)
      real(8) :: emg                            ! Ground (soil) emissivity
      real(8) :: sumabs                         ! Absorbed radiation for energy conservation check
      real(8) :: error                          ! Error check
      real(8) :: omega                          ! Leaf scattering coefficient
      real(8) :: rho                            ! Leaf reflectance
      real(8) :: tau                            ! Leaf transmittance
      real(8) :: trand                          ! Term for longwave radiation transmitted by layer
      real(8) :: refld                          ! Term for longwave radiation reflected by layer
      real(8) :: ir_source_sun                  ! Longwave radiation emitted by sunlit leaf (W/m2)
!      real(8) :: ir_source_sha                 ! Longwave radiation emitted by shaded leaf (W/m2)
      real(8) :: ir_source(Max_hgt)             ! Longwave radiation emitted by leaf layer (W/m2)
      integer :: m                              ! Index to the tridiagonal matrix
      real(8) :: aic, bic                       ! Intermediate terms for tridiagonal matrix
      real(8) :: eic, fic                       ! Intermediate terms for tridiagonal matrix
      integer, parameter :: neq = (Max_hgt+1)*2 ! Number of tridiagonal equations to solve
      real(8) :: atri(neq), btri(neq)           ! Entries in tridiagonal matrix
      real(8) :: ctri(neq), dtri(neq)           ! Entries in tridiagonal matrix
      real(8) :: utri(neq)                      ! Tridiagonal solution
      real(8) :: irabs                          ! Absorbed longwave flux (W/m2 ground)
      real(8) :: irup(0:Max_hgt)                ! Upward longwave flux above canopy layer (W/m2 ground)
      real(8) :: irdn(0:Max_hgt)                ! Downward longwave flux onto canopy layer (W/m2 ground)
      !---------------------------------------------------------------------

      ! Zero out radiative fluxes for all layers

      irup(0) = 0.0d0
      irdn(0) = 0.0d0

      do ic = 1, Max_hgt
         irup(ic) = 0.0;
         irdn(ic) = 0.0;
         irleaf(ic) = 0.0;
      end do

      ! Ground (soil) emissivity

      emg = 0.96d0

      !------------------------------------------------------------------
      ! Leaf scattering coefficient and terms for longwave radiation reflected
      ! and transmitted by a layer
      !------------------------------------------------------------------

      omega = 1.0d0 - emleaf

      ! Intercepted radiation is reflected

      rho = omega
      tau = 0.0d0

      ! Intercepted radiation is both reflected and transmitted

!      rho = omega * 0.5_r8
!      tau = omega * 0.5_r8

      !------------------------------------------------------------------
      ! Emitted longwave radiation is weighted average of sunlit and shaded leaves
      !------------------------------------------------------------------

      do ic = nbot, ntop
         ir_source_sun = emleaf * sb * tleaf_plot**4
!         ir_source_sha = emleaf(patch%itype(p)) * sb * tleaf(p,ic,isha)**4
!         ir_source(ic) = (ir_source_sun * fracsun(p,ic) + ir_source_sha * fracsha(p,ic)) * (1._r8 - td(p,ic))
         ir_source(ic) = ir_source_sun * (1.0d0 - td(ic))
      end do

      !------------------------------------------------------------------
      ! Set up and solve tridiagonal system of equations for upward and downward fluxes
      !------------------------------------------------------------------

      ! There are two equations for each leaf layer and the soil. The first
      ! equation is the upward flux and the second equation is the downward flux.

      m = 0

      ! Soil: upward flux

      m = m + 1
      atri(m) = 0.0d0
      btri(m) = 1.0d0
      ctri(m) = -(1.0d0 - emg)
      dtri(m) = emg * sb * tg**4

      ! Soil: downward flux

      refld = (1.0d0 - td(nbot)) * rho
      trand = (1.0d0 - td(nbot)) * tau + td(nbot)
      aic = refld - trand * trand / refld
      bic = trand / refld

      m = m + 1
      atri(m) = -aic
      btri(m) = 1.0d0
      ctri(m) = -bic
      dtri(m) = (1.0d0 - bic) * ir_source(nbot)

      ! Leaf layers, excluding top layer

      do ic = nbot, ntop-1

         ! Upward flux

         refld = (1.0d0 - td(ic)) * rho
         trand = (1.0d0 - td(ic)) * tau + td(ic)
         fic = refld - trand * trand / refld
         eic = trand / refld

         m = m + 1
         atri(m) = -eic
         btri(m) = 1.0d0
         ctri(m) = -fic
         dtri(m) = (1.0d0 - eic) * ir_source(ic)

         ! Downward flux

         refld = (1.0d0 - td(ic+1)) * rho
         trand = (1.0d0 - td(ic+1)) * tau + td(ic+1)
         aic = refld - trand * trand / refld
         bic = trand / refld

         m = m + 1
         atri(m) = -aic
         btri(m) = 1.0d0
         ctri(m) = -bic
         dtri(m) = (1.0d0 - bic) * ir_source(ic+1)

      end do

      ! Top canopy layer: upward flux

      ic = ntop
      refld = (1.0d0 - td(ic)) * rho
      trand = (1.0d0 - td(ic)) * tau + td(ic)
      fic = refld - trand * trand / refld
      eic = trand / refld

      m = m + 1
      atri(m) = -eic
      btri(m) = 1.0d0
      ctri(m) = -fic
      dtri(m) = (1.0d0 - eic) * ir_source(ic)

      ! Top canopy layer: downward flux

      m = m + 1
      atri(m) = 0.0d0
      btri(m) = 1.0d0
      ctri(m) = 0.0d0
      dtri(m) = irsky

      ! Solve tridiagonal system of equations for upward and downward fluxes

      call tridiag (atri, btri, ctri, dtri, utri, m)

      ! Now copy the solution (utri) to the upward (irup) and downward (irdn)
      ! fluxes for each layer
      ! irup =  Upward longwave flux above layer
      ! irdn =  Downward longwave flux onto layer

      m = 0

      ! Soil fluxes

      m = m + 1
      irup(0) = utri(m)
      m = m + 1
      irdn(0) = utri(m)

      ! Leaf layer fluxes

      do ic = nbot, ntop
         m = m + 1
         irup(ic) = utri(m)
         m = m + 1
         irdn(ic) = utri(m)
      end do

      !------------------------------------------------------------------
      ! Compute fluxes
      !------------------------------------------------------------------

      ! Absorbed longwave radiation for ground (soil)

      irsoi = irdn(0) - irup(0)

      ! Leaf layer fluxes

      irveg = 0.0d0

      do ic = nbot, ntop

         ! Absorbed longwave radiation for layer. Note special case for first
         ! leaf layer, where the upward flux from below is from the ground.
         ! The ground is ic=0, but nbot-1 will not equal 0 if there are lower
         ! canopy layers without leaves.

         if (ic == nbot) then
            icm1 = 0
         else
            icm1 = ic - 1
         end if
         irabs = emleaf * (irdn(ic)+irup(icm1)) * (1.0d0 - td(ic)) - 2.0d0 * ir_source(ic)
         irleaf(ic) = irabs / dpai_layer_sum(ic)

         ! Sum longwave radiation absorbed by vegetation

         irveg = irveg + irabs

      end do

      ! Canopy emitted longwave radiation

      ircan = irup(ntop)

      !------------------------------------------------------------------
      ! Conservation check
      !------------------------------------------------------------------

      ! Total radiation balance: absorbed = incoming - outgoing

      sumabs = irsky - ircan
      error = sumabs - (irveg + irsoi)
      if (abs(error) > 1.e-03) then
         write(*,*) 'Error: total longwave radiation conservation error'
      end if

    end subroutine longwave_rad

end module mod_radiation
