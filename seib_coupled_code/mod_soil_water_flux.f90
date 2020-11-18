module mod_soil_water_flux

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
  public :: soil_evap
  public :: soil_water_flux
  !-----------------------------------------------------------------------

contains
  !-----------------------------------------------------------------------
  subroutine soil_evap ( &
                                       ! *** Input ***
    ground_elev                   , &  ! Ground elevation (cm)
    soilresis                     , &  ! Soil evaporative resistance (s/m)
    albsoi_vis                    , &  ! Direct beam and diffuse albedo of ground for VIS (soil)
    albsoi_nir                    , &  ! Direct beam and diffuse albedo of ground for NIR (soil)
    moist                         , &  ! Relative soil moisture (-)
    bsw                           , &  ! Soil layer Clapp and Hornberger "b" parameter (-)
    psisat                        , &  ! Soil matric potential at saturation (MPa)
    time_series_tide              , &  ! Dimension of the virtual forest (m)
    time_series_tair              , &  ! Time-series air temperature (K)
    time_series_pair              , &  ! Time-series air pressure (Pa)
    time_series_wind              , &  ! Time-series wind speed (m/s)
    time_series_eair              , &  ! Time-series vapor pressure in air (Pa)
    time_series_rad_dir           , &  ! Time-series direct radiation at canopy top (W/m2)
    time_series_rad_dif           , &  ! Time-series diffused radiation at canopy top (W/m2)
                                       !
                                       ! *** From SEIB-DGVM ***
    STEP                          , &  ! Canopy layer thickness (m)
    day_of_year                   , &  ! Day of the year (1 - 365)
    ntop                          , &  ! Index for top leaf layer
    albvegb_vis                   , &  ! (Direct beam) vegetation albedo for VIS, non-horizontal leaves
    albvegb_nir                   , &  ! (Direct beam) vegetation albedo for NIR, non-horizontal leaves
    ir_tree                       , &  ! VIS-interruption-coefficient by tree corwn (0.0 - 1.0). lower => more absorptance.
    ir_tree_nir                   , &  ! NIR-interruption-coefficient by tree corwn (0.0 - 1.0). lower => more absorptance.
    irsoi                         , &  ! Absorbed longwave radiation, ground (W/m2)
    irveg                         , &  ! Absorbed longwave radiation, vegetation (W/m2)
                                       !
                                       ! *** Output ***
    rnsoi_day                     , &  ! Daily mean net radiation, ground (W/m2)
    rnveg_day                     , &  ! Daily mean net radiation, vegetation (W/m2)
    etsoi_day                       &  ! Daily soil evaporation rate (mm/day)
    )
    !
    ! !DESCRIPTION:
    ! Soil evaporation rate.
    !
    ! !USES:
    use mod_water_vapor
    !
    ! !ARGUMENTS:
    implicit none
    real(8), intent(in)  :: ground_elev, soilresis, albsoi_vis, albsoi_nir, moist, bsw, psisat
    real(8), intent(in)  :: time_series_tide(:), time_series_tair(:), time_series_pair(:)
    real(8), intent(in)  :: time_series_wind(:), time_series_eair(:)
    real(8), intent(in)  :: time_series_rad_dir(:), time_series_rad_dif(:)
    real(8), intent(in)  :: STEP
    integer, intent(in)  :: day_of_year, ntop
    real(8), intent(in)  :: albvegb_vis, albvegb_nir, ir_tree, ir_tree_nir, irsoi, irveg
    real(8), intent(out) :: rnsoi_day, rnveg_day, etsoi_day
    !
    ! !LOCAL VARIABLES:
    real(8), parameter :: denh2o = 1000.0d0            ! Density of liquid water (kg/m3)
    real(8), parameter :: grav = 9.80665d0             ! Gravitational acceleration (m/s2)
    real(8), parameter :: mmh2o = 18.02d0 / 1000.0d0   ! Molecular mass of water (kg/mol)is:starred
    real(8), parameter :: mmdry = 28.97d0 / 1000.0d0   ! Molecular mass of dry air (kg/mol)
    real(8), parameter :: rgas = 8.31446d0             ! Universal gas constant (J/K/mol)
    real(8), parameter :: cpd = 1005.0d0               ! Specific heat of dry air at constant pressure (J/kg/K)
    real(8), parameter :: cpw = 1846.0d0               ! Specific heat of water vapor at constant pressure (J/kg/K)
    real(8), parameter :: von_karman = 0.4d0           ! Von Karman constant
    real(8), parameter :: tfrz = 273.15d0              ! Freezing point of water (K)
    real(8) :: albcanb_vis                             ! Effective canopy albedo including soil for (direct beam) VIS
    real(8) :: albcanb_nir                             ! Effective canopy albedo including soil for (direct beam) NIR
    real(8) :: canopy_h                                ! Canopy height (m)
    real(8) :: dis_h                                   ! Displacement height (m)
    real(8) :: mrough_h                                ! Momentum roughness height (m)
    real(8) :: wrough_h                                ! Water vapor roughness height (m)
    real(8) :: head                                    ! Head of pressure (MPa/m)
    real(8) :: mat_psi                                 ! Soil matrix potential (mm)
    integer :: hour                                    ! 0 - 23
    integer :: hour_of_year
    real(8) :: tair                                    ! Air temperature (K)
    real(8) :: pair                                    ! Air pressure (Pa)
    real(8) :: wind                                    ! Wind speed (m/s)
    real(8) :: eair                                    ! Vapor pressure in air (Pa)
    real(8) :: rad_dir                                 ! Direct radiation at canopy top (W/m2)
    real(8) :: rad_dif                                 ! Diffused radiation at canopy top (W/m2)
    real(8) :: swsky_vis                               ! Atmospheric VIS solar radiation (W/m2)
    real(8) :: swsky_nir                               ! Atmospheric NIR solar radiation (W/m2)
    real(8) :: tide                                    ! Tide level (cm)
    real(8) :: rhomol                                  ! Molar density (mol/m3)
    real(8) :: cp                                      ! Specific heat of air at constant pressure (J/mol/K)
    real(8) :: qair                                    ! Specific humidity (kg/kg)
    real(8) :: rhoair                                  ! Air density (kg/m3)
    real(8) :: mmair                                   ! Molecular mass of air (kg/mol)
    real(8) :: tg                                      ! Soil surface temperature (K)
    real(8) :: tg_old                                  ! Soil surface temperature for previous timestep (K)
    real(8) :: t_soisno                                ! Soil temperature (K)
    real(8) :: ga                                      ! Canopy layer aerodynamic conductance for scalars (mol/m2/s)
    real(8) :: gws                                     ! Soil conductance for water vapor (mol H2O/m2/s)
    real(8) :: gw                                      ! Total conductance for water vapor (mol H2O/m2/s)
    real(8) :: swsoi                                   ! Solar radiation absorbed by ground (W/m2)
    real(8) :: rnsoi                                   ! Net radiation, ground (W/m2)
    real(8) :: swveg                                   ! Solar radiation absorbed by vegetaion (W/m2)
    real(8) :: rnveg                                   ! Net radiation, vegetation (W/m2)
    real(8) :: lambda                                  ! Latent heat of vaporization (J/mol)
    real(8) :: rhg                                     ! Relative humidity of airspace at soil surface (fraction)
    real(8) :: gamma                                   ! Psychrometric constant (Pa/K)
    real(8) :: esat                                    ! Saturation vapor pressure (Pa)
    real(8) :: desat                                   ! Derivative of saturation vapor pressure (Pa/K)
    real(8) :: qsat                                    ! Saturation vapor pressure of air (mol/mol)
    real(8) :: dqsat                                   ! Temperature derivative of saturation vapor pressure (mol/mol/K)
    real(8) :: num1, num2, num3, num4, den             ! Intermediate terms
    real(8) :: shsoi                                   ! Sensible heat flux, ground (W/m2)
    real(8) :: eg                                      ! Soil surface vapor pressure (Pa)
    real(8) :: lhsoi                                   ! Latent heat flux, ground (W/m2)
    real(8) :: gsoi                                    ! Soil heat flux (W/m2)
    real(8) :: err                                     ! Surface energy imbalance (W/m2)
    real(8) :: etsoi                                   ! Water vapor flux, ground (m3 H2O/m2/s)
    real(8), parameter :: thk = 1.58d0                 ! Soil layer thermal conductivity (W/m/K) .Clay soil (100% water content) in Bonan. Ecological Climatology: Concepts and Applications.
    real(8), parameter :: dz = 0.05d0                  ! Thickness of soil layer (m)
    !---------------------------------------------------------------------

    !---------------------------------------------------------------------
    !
    !    ! **************** Soil absorbed radiation **************** !
    !
    ! Based on Gpudriaan radiative transfer model,
    ! effective canopy albedo including soil is calculated as:
    !
    ! albcan = albveg + (albsoi - albveg) * exp(-2kbm*LAI)
    !
    ! where:
    !
    ! albcan = effective canopy albedo
    ! albveg = vegetation albedo
    ! albsoi = ground albedo
    ! kbm = Direct beam extinction coefficient adjusted for scattering
    ! LAI = Leaf area index
    !
    ! ir_tree provided by SEIB-DGVM corresponds to:
    !
    ! ir_tree = exp(-kbm*LAI)
    !
    ! Then, the effective canopy albedo is given by:
    !
    ! albcan = albveg + (albsoi - albveg) * ir_tree^2
    !
    ! Solar radiation absorbed by ground is:
    !
    ! swsoi = swsky * (1 - albcan) * exp(-kbm*LAI)
    !       = swsky * (1 - albcan) * ir_tree
    !
    ! where:
    !
    ! swsoi = Absorbed solar radiation, ground (W/m2)
    ! swsky = Atmospheric solar radiation (W/m2)
    !
    ! On the other hand, based on SEIB-DGVM, solar radiation absorbed
    ! by ground is given by:
    !
    ! swsoi = swsky * (1 - albsoi) * ir_tree
    !
    !    ! **************** Vegetation absorbed radiation **************** !
    !
    ! Solar radiation absorbed by vegetation is:
    !
    ! swveg = swsky * (1 - albcan) * (1 - exp(-kbm*LAI))
    !       = swsky * (1 - albcan) * (1 - ir_tree)
    !
    ! where:
    !
    ! swveg = Absorbed solar radiation, vegetation (W/m2)
    !
    !    ! **************** Soil evaporation **************** !
    !
    ! Soil evaporation rate is calculated by Penman-Monteith approach as
    ! in SEIB-DGVM.
    ! Note that in SEIB-DGVM, this rate is regarded as maximum evaporation
    ! which might be constrained depending on limited water availability.
    ! In mangroves, constrain of evaporation may be not likely occur
    ! due to frequent tidal water inundation.
    !
    !       s*Rn + cp * rhoair * vpd * ga
    ! E = ----------------------------------
    !        lambda(s + rw (1 + ga/gs))
    !
    ! where:
    !
    ! E = Soil evaporation rate (mol H2O/m2/s)
    ! s = slope of saturation vapor pressure (Pa/K)
    ! Rn = Net radiation, ground (W/m2)
    ! cp = Specific heat of air at constant pressure (J/kg/K)
    ! rhoair = Air density (kg/m3)
    ! vpd = Vapor pressure deficit (potential saturation deficit) of air (Pa)
    ! ga = Canopy layer aerodynamic conductance for scalars (m/s)
    ! lambda = Latent heat of vaporization (J/mol)
    ! rw = Psychrometric constant (Pa/K)
    ! gs = Soil conductance for water vapor (m/s)
    !---------------------------------------------------------------------

    ! Current ground temperature

    tg = 298.15d0
    tg_old = tg

    ! Effective canopy albedo including soil

    albcanb_vis = albvegb_vis + (albsoi_vis - albvegb_vis) * ir_tree**2
    albcanb_nir = albvegb_nir + (albsoi_nir - albvegb_nir) * ir_tree_nir**2

    ! Parameters for aerodynamic conductance (Perri et al. 2017)

    canopy_h = real(ntop) * STEP + 1.3d0
    dis_h = 0.75d0 * canopy_h
!   dis_h = 0.67d0 * canopy_h   ! Bonan et al. (2014)
    mrough_h = 0.10d0 * dis_h
!   mrough_h = 0.055 * canopy_h ! Bonan et al. (2014)
    wrough_h = 0.20d0 * mrough_h

    ! Head of pressure (MPa/m)

    head = denh2o * grav * 1.e-06

    ! Soil matrix potential (mm)
    ! Note: It is assumed that soil moisture is always saturated due to tidal water intrusion.

    mat_psi = psisat * moist ** (-1.0d0 * bsw)   ! Mpa
    mat_psi = mat_psi * head * 1000.d0           ! MPa -> m -> mm

    ! Initialization

    rnsoi_day = 0.0d0
    rnveg_day = 0.0d0
    etsoi_day = 0.0d0

    !---------------------------------------------------------------------
    ! Hourly computation
    !---------------------------------------------------------------------

    do hour = 0, 23

       ! Soil temperature (set to the same as soil surface temperature.)

       t_soisno = tg

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

       swsky_vis = rad_dir * 0.43d0 + rad_dif * 0.57d0
       swsky_nir = rad_dir * (1.0d0 - 0.43d0 - 0.04d0) &
                   + rad_dif * (1.0d0 - 0.57d0 - 0.04d0)

       !---------------------------------------------------------------------
       ! Set atmospheric related parameters
       !---------------------------------------------------------------------

       qair = mmh2o / mmdry * eair / (pair - (1.0d0 - mmh2o/mmdry) * eair)      ! Specific humidity (kg/kg)
       rhomol = pair / (rgas * tair)                                            ! Molar density (mol/m3)
       rhoair = rhomol * mmdry * (1.0d0 - (1.0d0 - mmh2o/mmdry) * eair / pair)  ! Air density (kg/m3)
       mmair = rhoair / rhomol                                                  ! Molecular mass of air (kg/mol)
       cp = cpd * (1.0d0 + (cpw/cpd - 1.0d0) * qair) * mmair;                   ! Specific heat of air at constant pressure (J/mol/K)
       lambda = 56780.3d0 - 42.84d0 * tair                                      ! Latent heat of vaporization
       gamma = cp * pair / lambda                                               ! Psychrometric constant (Pa/K)

       ! Tide level (cm)

       tide = time_series_tide(hour_of_year)

       !---------------------------------------------------------------------
       ! When the ground is exposed
       !---------------------------------------------------------------------
       if (ground_elev >= tide) then

          !---------------------------------------------------------------------
          ! Canopy layer aerodynamic conductance for scalars (mol/m2/s)
          !---------------------------------------------------------------------

          ! Daly et al. (2004)

          ga = 20.0d0 * 1.e-03     ! m/s
          ga = ga * rhomol         ! m/s -> mol H2O/m2/s
!          write(*,*) 'Daly et al (2004) ga', ga

          ! SEIB-DGVM waterbudget (c_aero)

          ga = (1.0d0 + 0.537d0 *wind) / 250.08d0  ! m/s
          ga = ga * rhomol                         ! m/s -> mol H2O/m2/s
!          write(*,*) 'SEIB-DGVM waterbudget (c_aero)', ga

          ! Sato et al. (2007) forest biome

          ga = (wind * von_karman**2.0d0) / (log(17.4d0))**2.0d0  ! m/s
          ga = ga * rhomol                                        ! m/s -> mol H2O/m2/s
!          write(*,*) 'Sato et al. (2007) forest biome', ga

          ! Sato et al. (2007) other biome

          ga = (wind * von_karman**2.0d0) / (log(146.0d0))**2.0d0  ! m/s
          ga = ga * rhomol                                         ! m/s -> mol H2O/m2/s
!          write(*,*) 'Sato et al. (2007) other biome', ga

          ! Perri et al. (2017)

          ga = (wind * von_karman**2.0d0) / &
               (log((canopy_h - dis_h) / mrough_h) * log((canopy_h - dis_h) / wrough_h)) ! m/s
          ga = ga * rhomol             ! m/s -> mol H2O/m2/s
!          write(*,*) 'Perri et al. (2017) ga', ga

          ! Soil conductance to water vapour diffusion

          gws = 1.0d0 / soilresis    ! s/m -> m/s
          gws = gws * rhomol         ! m/s -> mol H2O/m2/s
          gw = ga * gws / (ga + gws) ! Total conductance

          !---------------------------------------------------------------------
          ! Solar radiation absorbed by ground (W/m2)
          !---------------------------------------------------------------------

          ! Method of Bonan et al. (2014) Goudriaan radiative transfer model

          swsoi = swsky_vis * (1.0d0 - albcanb_vis) * ir_tree &
                  + swsky_nir * (1.0d0 - albcanb_nir) * ir_tree_nir
!          write(*,*) 'Bonan et al. (2014) albcanb_vis, albcanb_nir, swsoi', albcanb_vis, albcanb_nir, swsoi

          ! Method of SEIB-DGVM

!          swsoi = swsky_vis * (1.0d0 - albsoi_vis) * ir_tree &
!                  + swsky_nir * (1.0d0 - albsoi_nir) * ir_tree_nir
!          write(*,*) 'SEIB-DGVM albsoi_vis, albsoi_nir, swsoi', albsoi_vis, albsoi_nir, swsoi

          ! Net radiation (sum of soil absorbed short and longwave radiation)

          rnsoi = swsoi + irsoi

          !---------------------------------------------------------------------
          ! Soil evaporation rate based on energy balance
          ! <= Currently, the energy is not balanced.
          !---------------------------------------------------------------------

          ! Relative humidity in soil airspace

          rhg = exp(grav * mmh2o * mat_psi*1.e-03 / (rgas * t_soisno))

          ! Saturation vapor pressure at ground temperature (Pa -> mol/mol)

          call sat_vap (tfrz, tg, esat, desat)
          qsat = esat / pair
          dqsat = desat / pair

          ! Calculate soil surface temperature
          ! Based on Rn = Hg + lambda*E

          num1 = (cp * pair / lambda) * rnsoi
          num2 = cp * rhg * (esat - eair) * gw
          num3 = (cp * pair / lambda) * ga * cp
          tg = (num1 - num2)/num3 + tair
!          write(*,*) 'tg, tair, method1', tg, tair

          ! Calculate soil surface temperature
          ! Based on Rn = Hg + lambda*E + G

          num1 = cp * ga
          num2 = lambda * gw
          num3 = thk / dz
          num4 = rnsoi - num2 * rhg * (qsat - dqsat * tg) + num3 * t_soisno
          den = num1 + num2 * dqsat * rhg + num3
          tg = (num1*tair + num2*eair/pair + num4) / den
!          write(*,*) 'tg, method2', tg

          ! Sensible heat flux

          shsoi = cp * (tg - tair) * ga

          ! Latent heat flux

          eg = rhg * (esat + desat * (tg - tg_old))
          lhsoi = lambda / pair * (eg - eair) * gw

          ! Soil heat flux

          gsoi = thk * (tg - t_soisno) / dz

          ! Error check

          err = rnsoi - shsoi - lhsoi - gsoi
          if (abs(err) > 0.001d0) then
!             write(*,*) 'ERROR: Soil water flux energy balance error', err
          end if

          ! Water vapor flux: W/m2 -> mol H2O/m2/s

          etsoi = lhsoi / lambda               ! (mol H2O/m2/s)
          etsoi = etsoi * mmh2o                ! (kg H2O/m2/s)
          etsoi = etsoi / denh2o               ! (m3 H2O/m2/s)
          etsoi = max(etsoi, 0.0d0)            ! Prevent negative value. (needed?)

!write(*,*) 'etsoi based on energy balance, mm/day', etsoi*1000.0d0*60.0d0*60.0d0*24.0d0

          ! Save current ground temperature

          tg_old = tg

          !---------------------------------------------------------------------
          ! Soil evaporation rate based on Penman-Monteith approach
          !---------------------------------------------------------------------

          ! ga, gws. mol H2O/m2/s -> m/s

          ga = ga / rhomol
          gws = gws / rhomol

          ! cp. J/mol/K -> J/kg/K

          cp = cp / mmair

          ! Saturation vapor pressure at air temperature (Pa)

          call sat_vap (tfrz, tg, esat, desat)

          ! Calculate soil evaporation rate

          num1 = desat * rnsoi
          num2 = cp * rhoair * (esat - eair) * ga
          num3 = lambda * (desat + gamma * (1.0d0 + ga/gws))
          etsoi = (num1 + num2) / num3
          etsoi = etsoi * mmh2o                ! (kg H2O/m2/s)
          etsoi = etsoi / denh2o               ! (m3 H2O/m2/s)
          etsoi = max(etsoi, 0.0d0)            ! Prevent negative value. (needed?)

!write(*,*) 'etsoi mm/day, Penman-Monteith', etsoi*1000.0d0*60.0d0*60.0d0*24.0d0

       !---------------------------------------------------------------------
       ! When the ground is submerged.
       ! Assumed Rn_soil = 0
       !---------------------------------------------------------------------
       else
          rnsoi = 0.0d0
          tg = 298.15d0
          tg_old = tg
          etsoi = 0.0d0
       end if

       !---------------------------------------------------------------------
       ! Solar radiation absorbed by vegetation (W/m2)
       ! Not affected by tide, but the changes in soil albedo by tidal
       ! water inundation is neglected.
       !---------------------------------------------------------------------

       swveg = swsky_vis * (1.0d0 - albcanb_vis) * (1.0d0 - ir_tree) &
               + swsky_nir * (1.0d0 - albcanb_nir) * (1.0d0 - ir_tree_nir)

       ! Net radiation (sum of vegetation absorbed short and longwave radiation)

       rnveg = swveg + irveg

       ! Calculating daily mean net radiation, ground and vegetation (W/m2)

       rnsoi_day = rnsoi_day + rnsoi / 24.0d0
       rnveg_day = rnveg_day + rnveg / 24.0d0

       ! Calculating daily soil evaporation rate (mm/day)
       ! Multiplying 1000 is for (m3 H2O/m2/s) -> (mm H2O/s)
       ! Multiplying 60*60 is for (mm H2O/s) -> (mm H2O/hr)
       ! Summing every hour for (mm H2O/hr) -> (mm H2O/day)

       etsoi_day = etsoi_day + etsoi * 1000.0d0 * 60.0d0 * 60.0d0

!write(*,*) 'etsoi_day, hour, mm/day', etsoi_day, hour

    end do

  end subroutine soil_evap

  !-----------------------------------------------------------------------
  subroutine soil_water_flux ( &
                                       ! *** Input ***
    sal_ini                       , &  ! Initial value of pore-water salinity (psu -> mol/m3)
    din_ini                       , &  ! Initial value of DIN concentration in pore-water (mol N/m3)
    dip_ini                       , &  ! Initial value of DIP concentration in pore-water (mol P/m3)
    flux_fw                       , &  ! Freshwater inputs to soil (mm/hr)
    din_fw                        , &  ! DIN concentration in freshwater (mol N/m3)
    dip_fw                        , &  ! DIP concentration in freshwater (mol P/m3)
    sal_sw                        , &  ! Salinity in seawater (mol/m3)
    din_sw                        , &  ! DIN concentration in seawater (mol N/m3)
    dip_sw                        , &  ! DIP concentration in seawater (mol P/m3)
    root_filter                   , &  ! Filtering rate of salt at root surface (-)
    root_depth                    , &  ! Rooting depth (m)
                                       !
                                       ! *** From SEIB-DGVM ***
    flux_calc_switch              , &  ! Switch of soil water flux calculation (1: ON, 0: OFF)
    Max_loc                       , &  ! Dimension of the virtual forest (m)
    Max_no                        , &  ! Maximum number of individual stands
    tree_exist                    , &  ! Flag of tree presence
    pft                           , &  ! Species index (1: Rh, 2: Br)
    mass_leaf                     , &  ! Stand-level leaf biomass (g leaf/tree)
    plant_water_uptake            , &  ! Daily plant water uptake rate at stand-level (m3 H2O/tree/day)
    etsoi_day                     , &  ! Daily soil evaporation rate (mm/day)
                                       !
                                       ! *** In/Output ***
    etveg_day                     , &  ! Daily transpiration rate (mm/day) <= output only
    sal                           , &  ! Pore-water salinity (mol/m3)
    din                           , &  ! DIN concentration in pore-water (mol N/m3)
    dip                             &  ! DIP concentration in pore-water (mol P/m3)
    )
    !
    ! !DESCRIPTION:
    ! Box model for soil water flux.
    ! Evaporation from soil is not considered yet. Needs to be added in the future.
    !
    ! !USES:
    use time_counter, only : counter
    !
    ! !ARGUMENTS:
    implicit none
    real(8), intent(in)  :: sal_ini, din_ini, dip_ini
    real(8), intent(in)  :: flux_fw, din_fw, dip_fw, sal_sw, din_sw, dip_sw
    real(8), intent(in)  :: root_filter(:), root_depth
    integer, intent(in)  :: flux_calc_switch, Max_loc, Max_no
    logical, intent(in)  :: tree_exist(:)
    integer, intent(in)  :: pft(:)
    real(8), intent(in)  :: mass_leaf(:), plant_water_uptake(:), etsoi_day
    real(8), intent(out) :: etveg_day
    real(8), intent(inout) :: sal, din, dip
    !
    ! !LOCAL VARIABLES:
    real(8) :: box_vol                      ! Volume of the box (m3)
    real(8) :: sal_box                      ! Amounts of salt in the box (mol)
    real(8) :: din_box                      ! Amounts of N in the box (mol)
    real(8) :: dip_box                      ! Amounts of P in the box (mol)
    real(8) :: wat_remove                   ! Daily water uptake by plants (m3 H2O/day)
    real(8) :: sal_remove                   ! Daily salt uptake by plants (mol/day)
    real(8) :: din_remove                   ! Daily N uptake by plants (mol/day)
    real(8) :: dip_remove                   ! Daily P uptake by plants (mol/day)
    integer :: no                           ! Tree index
    real(8) :: wat_input                    ! Freshwater inputs (m3 H2O/day)
    real(8) :: din_input                    ! N inputs (mol/day)
    real(8) :: dip_input                    ! P inputs (mol/day)
    real(8) :: wat_vol                      ! Water volume (m3 H2O)
    real(8) :: wat_recharge                 ! Recharge of water by sea water intrusion (m3 H2O/day)
    real(8) :: sal_recharge                 ! Recharge of salt by sea water intrusion (mol/day)
    real(8) :: din_recharge                 ! Recharge of N by sea water intrusion (mol/day)
    real(8) :: dip_recharge                 ! Recharge of P by sea water intrusion (mol/day)
    !---------------------------------------------------------------------

    !---------------------------------------------------------------------
    ! Amounts of salt/N/P in the box (mol)
    !---------------------------------------------------------------------

    ! Volume of the box (m3)

    box_vol = real(Max_loc) * real(Max_loc) * root_depth

    ! Amounts of salt/N/P in the box (mol)

    sal_box = sal * box_vol
    din_box = din * box_vol
    dip_box = dip * box_vol

    !---------------------------------------------------------------------
    ! Daily water/salt/N/P removal by plants
    !---------------------------------------------------------------------

    ! Computation of all trees

    wat_remove = 0.0d0
    sal_remove = 0.0d0
    din_remove = 0.0d0
    dip_remove = 0.0d0

    do no = 1, Max_no

       if ( .not. tree_exist(no) ) cycle
       if ( mass_leaf(no) <= 0.0d0 ) cycle

       ! Daily water uptake by plants (m3 H2O/day)

       wat_remove = wat_remove + plant_water_uptake(no)

       ! Daily salt/N/P uptake by plants (mol/day)

       sal_remove = sal_remove + plant_water_uptake(no) * sal * (1.0d0 - root_filter(pft(no)))
       din_remove = din_remove + plant_water_uptake(no) * din
       dip_remove = dip_remove + plant_water_uptake(no) * dip

    end do

    ! Daily transpiration rate (mm/day)
    ! Deviding by Max_loc^2 is for m3/day -> m/day
    ! Multiplying 1000 is for m/day -> mm/day

    etveg_day = wat_remove * 1000.0d0 / (real(Max_loc)**2.0d0)

write(*,*) 'Transpiration rate (mm/day)', etveg_day

    !---------------------------------------------------------------------
    ! Add daily soil evaporation rate (m3 H2O/day)
    !---------------------------------------------------------------------

    ! Multiplying 1.e-03 for mm/day -> m/day
    ! Multiplying Max_loc^2 is for m/day -> m3/day

    wat_remove = wat_remove + etsoi_day * 1.e-03 * real(Max_loc) * real(Max_loc)

write(*,*) 'Soil evaporation rate (mm/day)', etsoi_day

    !---------------------------------------------------------------------
    ! Inputs of freshwater/N/P to the box
    !---------------------------------------------------------------------

    ! Freshwater inputs (m3 H2O/day)

    wat_input = flux_fw * real(Max_loc) * real(Max_loc) * 24.0d0 * 1.e-03

    ! N and P inputs (mol/day)

    din_input = din_fw * wat_input
    dip_input = dip_fw * wat_input

    !---------------------------------------------------------------------
    ! Case the water overflows by excess freshwater inputs
    ! 淡水が土の中の塩分・栄養塩濃度を薄めて、そのあと余剰分がボックスから出て行く。
    !---------------------------------------------------------------------

    if (wat_input > wat_remove .and. flux_calc_switch == 1) then

write(*,*) 'Freshwater overflow'

       ! No sea water recharge

       wat_recharge = 0.0d0

!!! ------------------------------- Option 1 ------------------------------ !!!
!
!       ! Water volume just before the overflow (m3 H2O)
!
!       wat_vol = box_vol - wat_remove + wat_input
!
!       ! Amounts of salt/N/P in the box just before the overflow (mol)
!
!       sal_box = sal_box - sal_remove
!       din_box = din_box - din_remove + din_input
!       dip_box = dip_box - dip_remove + dip_input
!
!       ! Recalculation of pore-water salinity and nutrient concentrations (mol/m3)
!
!       sal = sal_box / wat_vol
!       din = din_box / wat_vol
!       dip = dip_box / wat_vol
!
!!! ------------------------------- Option 2 ------------------------------ !!!

       ! Water volume just before the overflow (m3 H2O)

       wat_vol = box_vol - wat_remove

       ! Amounts of salt/N/P in the box just before the overflow (mol)

       sal_box = sal_box - sal_remove
       din_box = din_box - din_remove
       dip_box = dip_box - dip_remove

       ! Amounts of salt/N/P in the box after the overflow (mol)

       sal_box = sal_box - sal_box * (wat_input - wat_remove)/box_vol
       din_box = din_box - din_box * (wat_input - wat_remove)/box_vol + din_input
       dip_box = dip_box - dip_box * (wat_input - wat_remove)/box_vol + dip_input

       ! Recalculation of pore-water salinity and nutrient concentrations (mol/m3)

       sal = sal_box / box_vol
       din = din_box / box_vol
       dip = dip_box / box_vol
!!! ----------------------------------------------------------------------- !!!

    !---------------------------------------------------------------------
    ! Case the water is recharged by sea water intrusion
    !---------------------------------------------------------------------
    else if (flux_calc_switch == 1) then

write(*,*) 'Seawater recharge'

       ! Recharge of water by sea water intrusion (m3 H2O/day)

       wat_recharge = wat_remove - wat_input

       ! Recharge of salt/N/P by sea water intrusion (mol/day)

       sal_recharge = sal_sw * wat_recharge
       din_recharge = din_sw * wat_recharge
       dip_recharge = dip_sw * wat_recharge

       ! Water volume after water recharge

       wat_vol = box_vol - wat_remove + wat_input + wat_recharge

       ! Amounts of salt/N/P in the box after water recharge (mol)

       sal_box = sal_box - sal_remove + sal_recharge
       din_box = din_box - din_remove + din_input + din_recharge
       dip_box = dip_box - dip_remove + dip_input + dip_recharge

       ! Recalculation of pore-water salinity and nutrient concentrations (mol/m3)

       sal = sal_box / wat_vol
       din = din_box / wat_vol
       dip = dip_box / wat_vol

    !---------------------------------------------------------------------
    ! When the flux calculation is off
    !---------------------------------------------------------------------
    else
       sal = sal_ini
       din = din_ini
       dip = dip_ini
    end if

  end subroutine soil_water_flux

end module mod_soil_water_flux
