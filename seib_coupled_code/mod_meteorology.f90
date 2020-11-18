module mod_meteorology

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
  public  :: meteorology_read
  private :: jra55_interpolate
  public  :: t_growth_read
  public  :: wind_profile
  !
  !
  ! !LOCAL VARIABLES
  private
  real(8), parameter :: tfrz = 273.15d0  ! Freezing point of water (K)
  !-----------------------------------------------------------------------

contains

  !-----------------------------------------------------------------------
  subroutine meteorology_read ( &
                                       ! *** From SEIB-DGVM ***
    jra55_use_switch              , &  ! Switch of JRA-55 data (1: use, 0: no use)
    jra55_timezone                , &  ! Time zone of JRA-55 data (correction for local time)
    jra55_lat                     , &  ! Latitude of JRA-55 data
    Fn_jra55                      , &  ! File name of JRA-55 data
    Fn                              &  ! File number for hourly solar radiation
    )
    !
    ! !DESCRIPTION:
    ! Read JMA observation data or JRA-55 data.
    !
    ! !USES:
    use data_structure, only : time_series_tair, time_series_pair, time_series_wind, &
                               time_series_eair, time_series_rad_dir, time_series_rad_dif, &
                               time_series_long_rad, daily_cloud
    !
    ! !ARGUMENTS:
    implicit none
    integer, intent(in) :: jra55_use_switch
    integer, intent(in) :: jra55_timezone
    real(8), intent(in) :: jra55_lat
    character*128       :: Fn_jra55
    integer, intent(in) :: Fn
    !
    ! !LOCAL VARIABLES:
    !
    !---------------------------------------------------------------------

    ! When the case of use of JRA-55 data

    if (jra55_use_switch == 1) then

       call jra55_interpolate ( &
                                         ! *** From SEIB-DGVM ***
       jra55_timezone               , &  ! Time zone of JRA-55 data (correction for local time)
       jra55_lat                    , &  ! Latitude of JRA-55 data
       Fn_jra55                     , &  ! File name of JRA-55 data
       Fn                             &  ! File number for hourly solar radiation
       )

    ! When the case of use of JMA observation data
    else

      ! Read hourly data: Jan1 2014 00:00 - Dec31 2014 23:00

      ! Air temperature (K)
      open (82, file = './input/Meteorological_data/Air_temperature.txt'   , status = 'old')
      read (82,*) time_series_tair
      close (82)

      ! Air pressure (Pa)
      open (83, file = './input/Meteorological_data/Air_pressure.txt'      , status = 'old')
      read (83,*) time_series_pair
      close (83)

      ! Wind speed (m/s)
      open (84, file = './input/Meteorological_data/Wind.txt'              , status = 'old')
      read (84,*) time_series_wind
      close (84)

      ! Vapor pressure in air (Pa)
      open (85, file = './input/Meteorological_data/Vapor_pressure.txt'    , status = 'old')
      read (85,*) time_series_eair
      close (85)

      ! Direct radiation at canopy top (W/m2)
      open (86, file = './input/Meteorological_data/Direct_radiation.txt'  , status = 'old')
      read (86,*) time_series_rad_dir
      close (86)

      ! Diffused radiation at canopy top (W/m2)
      open (87, file = './input/Meteorological_data/Diffused_radiation.txt', status = 'old')
      read (87,*) time_series_rad_dif
      close (87)

      ! Atmospheric downward longwave radiation (W/m2)
      open (88, file = './input/Meteorological_data/Atmospheric_longwave_radiation.txt', status = 'old')
      read (88,*) time_series_long_rad
      close (88)

      ! Daily mean cloud (0 - 10)
      open (89, file = './input/Meteorological_data/Daily_cloud.txt', status = 'old')
      read (89,*) daily_cloud
      close (89)

   end if

  end subroutine meteorology_read

  !-----------------------------------------------------------------------
  subroutine jra55_interpolate ( &
                                       ! *** From SEIB-DGVM ***
    jra55_timezone                , &  ! Time zone of JRA-55 data (correction for local time)
    jra55_lat                     , &  ! Latitude of JRA-55 data
    Fn_jra55                      , &  ! File name of JRA-55 data
    Fn                              &  ! File number for hourly solar radiation
    )
    !
    ! !DESCRIPTION:
    ! Interpolate JRA-55 3-hourly data to hourly.
    ! Compute vapor pressure and solar radiation.
    ! The computation of solar radiation is based on Zillman equation
    ! which is used in ROMS.
    ! The computation of direct and diffused radiation is the same as in SEIB-DGVM.
    !
    ! !USES:
    use data_structure, only : time_series_tair, time_series_pair, time_series_wind, &
                               time_series_eair, time_series_rad_dir, time_series_rad_dif, &
                               daily_cloud, PI, DtoR
    use mod_water_vapor
    !
    ! !ARGUMENTS:
    implicit none
    integer, intent(in) :: jra55_timezone
    real(8), intent(in) :: jra55_lat
    character*128       :: Fn_jra55
    integer, intent(in) :: Fn
    !
    ! !LOCAL VARIABLES:
    real(8), dimension (2920, 5) :: jra55_data          ! 3-hourly cloud, wind, tair, rh, pair
    real(8), dimension (2920*3)  :: cloud_hourly      ! Interpolated hourly cloud (fraction)
    real(8), dimension (2920*3)  :: wind_hourly       ! Interpolated hourly wind (m/s)
    real(8), dimension (2920*3)  :: tair_hourly       ! Interpolated hourly tair (degree)
    real(8), dimension (2920*3)  :: rh_hourly         ! Interpolated hourly rh (fraction)
    real(8), dimension (2920*3)  :: pair_hourly       ! Interpolated hourly pair (Pa)
    real(8), dimension (2920*3)  :: cloud_hourly_adj  ! Interpolated hourly cloud adjusted for local time (fraction)
    real(8), dimension (2920*3)  :: rh_hourly_adj     ! Interpolated hourly rh adjusted for local time (fraction)
    integer :: t, t0                                  ! Time index
    real(8) :: esat                                   ! Saturation vapor pressure (Pa)
    real(8) :: desat                                  ! Temperature derivative of saturation vapor pressure (Pa/K)
    integer :: day_of_year                            ! Day of year
    integer :: hour                                   ! Hour
    real(8) :: sl_dec                                 ! Solar declination (deg)
    real(8) :: h_angle                                ! Hour angle (radian)
    real(8) :: cos_zen                                ! Cosine of zenith angle
    real(8) :: rad_stratosphere                       ! Incoming shortwave radiation under cloudless skies (W/m2)
    real(8) :: rad                                    ! Incoming shortwave radiation corrected for cloudiness (W/m2)
    real(8), parameter :: sl_con = 1353               ! Solar constant (W/m2)
    !---------------------------------------------------------------------

    ! Read hourly data: GMT(+0) Jan1 2014 00:00 - Dec31 2014 21:00

    open (82, file = Fn_jra55, status = 'old')
    do t = 1, 2920
       read (82,*) jra55_data(t,1), jra55_data(t,2), jra55_data(t,3), jra55_data(t,4), jra55_data(t,5)
    end do
    close (82)

    !---------------------------------------------------------------------
    ! Interpolate 3 hourly data to hourly data
    !---------------------------------------------------------------------

    do t = 1, 2920

       cloud_hourly(3*(t-1)+1) = jra55_data(t, 1)
       wind_hourly (3*(t-1)+1) = jra55_data(t, 2)
       tair_hourly (3*(t-1)+1) = jra55_data(t, 3)
       rh_hourly   (3*(t-1)+1) = jra55_data(t, 4)
       pair_hourly (3*(t-1)+1) = jra55_data(t, 5)

       if (t < 2920) then

          cloud_hourly(3*(t-1)+2) = jra55_data(t, 1) + (jra55_data(t+1, 1) - jra55_data(t, 1))/3.0d0
          wind_hourly (3*(t-1)+2) = jra55_data(t, 2) + (jra55_data(t+1, 2) - jra55_data(t, 2))/3.0d0
          tair_hourly (3*(t-1)+2) = jra55_data(t, 3) + (jra55_data(t+1, 3) - jra55_data(t, 3))/3.0d0
          rh_hourly   (3*(t-1)+2) = jra55_data(t, 4) + (jra55_data(t+1, 4) - jra55_data(t, 4))/3.0d0
          pair_hourly (3*(t-1)+2) = jra55_data(t, 5) + (jra55_data(t+1, 5) - jra55_data(t, 5))/3.0d0

          cloud_hourly(3*(t-1)+3) = jra55_data(t, 1) + (jra55_data(t+1, 1) - jra55_data(t, 1))*2.0d0/3.0d0
          wind_hourly (3*(t-1)+3) = jra55_data(t, 2) + (jra55_data(t+1, 2) - jra55_data(t, 2))*2.0d0/3.0d0
          tair_hourly (3*(t-1)+3) = jra55_data(t, 3) + (jra55_data(t+1, 3) - jra55_data(t, 3))*2.0d0/3.0d0
          rh_hourly   (3*(t-1)+3) = jra55_data(t, 4) + (jra55_data(t+1, 4) - jra55_data(t, 4))*2.0d0/3.0d0
          pair_hourly (3*(t-1)+3) = jra55_data(t, 5) + (jra55_data(t+1, 5) - jra55_data(t, 5))*2.0d0/3.0d0

       else

          cloud_hourly(3*(t-1)+2) = jra55_data(t, 1) + (jra55_data(1, 1) - jra55_data(t, 1))/3.0d0
          wind_hourly (3*(t-1)+2) = jra55_data(t, 2) + (jra55_data(1, 2) - jra55_data(t, 2))/3.0d0
          tair_hourly (3*(t-1)+2) = jra55_data(t, 3) + (jra55_data(1, 3) - jra55_data(t, 3))/3.0d0
          rh_hourly   (3*(t-1)+2) = jra55_data(t, 4) + (jra55_data(1, 4) - jra55_data(t, 4))/3.0d0
          pair_hourly (3*(t-1)+2) = jra55_data(t, 5) + (jra55_data(1, 5) - jra55_data(t, 5))/3.0d0

          cloud_hourly(3*(t-1)+3) = jra55_data(t, 1) + (jra55_data(1, 1) - jra55_data(t, 1))*2.0d0/3.0d0
          wind_hourly (3*(t-1)+3) = jra55_data(t, 2) + (jra55_data(1, 2) - jra55_data(t, 2))*2.0d0/3.0d0
          tair_hourly (3*(t-1)+3) = jra55_data(t, 3) + (jra55_data(1, 3) - jra55_data(t, 3))*2.0d0/3.0d0
          rh_hourly   (3*(t-1)+3) = jra55_data(t, 4) + (jra55_data(1, 4) - jra55_data(t, 4))*2.0d0/3.0d0
          pair_hourly (3*(t-1)+3) = jra55_data(t, 5) + (jra55_data(1, 5) - jra55_data(t, 5))*2.0d0/3.0d0

       end if

    end do

    !---------------------------------------------------------------------
    ! Adjusting for local time
    !---------------------------------------------------------------------

    do t = 1, 2920*3

       if (t <= jra55_timezone) then

          do t0 = 1,jra55_timezone

             cloud_hourly_adj(t0) = cloud_hourly(2920*3 - jra55_timezone+t0)
             time_series_wind(t0) = wind_hourly (2920*3 - jra55_timezone+t0)
             time_series_tair(t0) = tair_hourly (2920*3 - jra55_timezone+t0) + tfrz
             rh_hourly_adj   (t0) = rh_hourly   (2920*3 - jra55_timezone+t0)
             time_series_pair(t0) = pair_hourly (2920*3 - jra55_timezone+t0)

          end do
       else

          cloud_hourly_adj(t) = cloud_hourly(t - jra55_timezone)
          time_series_wind(t) = wind_hourly (t - jra55_timezone)
          time_series_tair(t) = tair_hourly (t - jra55_timezone) + tfrz
          rh_hourly_adj   (t) = rh_hourly   (t - jra55_timezone)
          time_series_pair(t) = pair_hourly (t - jra55_timezone)
       end if

    end do

    !---------------------------------------------------------------------
    ! Compute vapor pressure and solar radiation
    !---------------------------------------------------------------------

    do t = 1, 2920*3

       ! Saturation vapor pressure in air (Pa)

       call sat_vap (tfrz, time_series_tair(t), esat, desat)

       ! Vapor pressure in air (Pa)

       time_series_eair(t) = rh_hourly_adj(t) * esat

       ! Day of year and hour

       day_of_year = int((dble(t)-1.0d0)/24.0d0) + 1
       hour = mod(t-1, 24)

       ! Solar declination (deg)

       sl_dec = 23.44d0 * cos((172.0d0 - dble(day_of_year)) * 2.0d0 * PI / 365.0d0)

       ! Hour angle (radian)

       h_angle = dble(12 - hour) * PI / 12.0d0

       ! Cosine of zenith angle

       cos_zen = sin(jra55_lat * DtoR) * sin(sl_dec * DtoR) + &
                 cos(jra55_lat * DtoR) * cos(sl_dec * DtoR) * cos(h_angle)

       ! Incoming radiation under cloudless skies (W/m2)

       rad_stratosphere = (sl_con * cos_zen**2.0d0) / &
                          ((cos_zen + 2.70d0)*time_series_eair(t)*1.e-05 + 1.085d0*cos_zen + 0.10d0)
       if (abs(acos(cos_zen)) > PI/2.0d0) then
          rad_stratosphere = 0.0d0
       end if

       ! Correction of radiation for cloudiness (W/m2)

       rad = rad_stratosphere * (1.0d0 - 0.6d0 * cloud_hourly_adj(t)**3.0d0)

       ! Direct and diffused radiation (W/m2)

       time_series_rad_dif(t) = max(0.0d0 , rad * (0.958d0 - 0.982d0*(rad/rad_stratosphere)))
       time_series_rad_dir(t) = max(0.0d0 , rad - time_series_rad_dif(t))

       ! Error check purpose

!       write(*,*) 't, day_of_year', t, day_of_year, hour
!       write(*,*) 'cloud, eair', cloud_hourly_adj(t), time_series_eair(t)
!       write(*,*) 'solar radiation', rad, sl_dec
       write(*,*) 'doy_of_year, hour, rad', day_of_year, hour, rad
       write (Fn, '(f12.5,a)') rad
!       write (Fn, *) rad

    end do

    !---------------------------------------------------------------------
    ! Daily mean cloud (0 - 10) : cloudness at midday
    !---------------------------------------------------------------------

    do t = 1, 365

       daily_cloud(t) = cloud_hourly_adj((t-1)*24 + 13)

!       write(*,*) 'daily_cloud', daily_cloud(t)

    end do

  end subroutine jra55_interpolate

  !-----------------------------------------------------------------------
  subroutine t_growth_read ( &
                                       ! *** Input ***
    n_spe                         , &  ! Number of mangrove species
    t_acclim                      , &  ! Flag of temperature acclimation (1: Use, 0: No use)
    time_series_tair              , &  ! Time-series air temperature (K)
                                       !
                                       ! *** From SEIB-DGVM ***
    doy                           , &  ! Day of the year (1 - 365)
                                       !
                                       ! *** In/Output
    t_growth                        &  ! Growth temperature (degree), mean temperature for 30 days before
    )
    !
    ! !DESCRIPTION:
    ! Read growth temperature when temperature acclimation is active.
    !
    ! !USES:
    !
    !
    ! !ARGUMENTS:
    implicit none
    integer, intent(in)    :: n_spe
    integer, intent(in)    :: t_acclim(:)
    real(8), intent(in)    :: time_series_tair(:)
    integer, intent(in)    :: doy
    real(8), intent(inout) :: t_growth(:)
    !
    ! !LOCAL VARIABLES:
    integer :: p                             ! Species index
    integer :: hour_of_year                  ! Hour of year
    integer :: t, t0                         ! Time index
    !---------------------------------------------------------------------

    do p = 1, n_spe

       ! When temperature acclimation is active.

       if (t_acclim(p) == 1) then

          ! Initialize t_growth

          t_growth(p) = 0.0d0

          ! Hour of year

          hour_of_year = (doy - 1)*24 + 1

          ! Hour of year of the start of the averaging (30 days before the day)

          if (hour_of_year <= 720) then

             t0 = 8760 + (hour_of_year - 24*30)
          else
             t0  = hour_of_year - 24*30
          end if

          ! Averaging of the temperature

          do t = 1, 720

             if (t0 + t - 1 <= 8760) then

                t_growth(p) = t_growth(p) + (time_series_tair(t0 + t - 1)-tfrz)/720;

              else

                t_growth(p) = t_growth(p) + (time_series_tair(t0 + t - 1 - 8760)-tfrz)/720;

              end if

          end do

! write(*,*) 'doy, tair, t_growth', doy, time_series_tair(hour_of_year)-tfrz, t_growth(1)

       end if
    end do

  end subroutine t_growth_read

  !-----------------------------------------------------------------------
  subroutine wind_profile ( &
                                       ! *** Input ***
    n_spe                         , &  ! Number of mangrove species
    dleaf                         , &  ! Leaf dimension (m)
                                       !
                                       ! *** From SEIB-DGVM ***
    STEP                          , &  ! Canopy layer thickness (m)
    dpai_layer_sum                , &  ! Sum of layer leaf area index for each forest layer (m2 leaf/me ground)
    nbot                          , &  ! Index for bottom leaf layer
    ntop                          , &  ! Index for top leaf layer
                                       !
                                       ! *** In/Output ***
    rwind_profile                   &  ! Relative wind speed profile to the canopy top (-)
    )
    !
    ! !DESCRIPTION:
    ! Calculate wind profile within canopy.
    ! This is based on the Goudriaan (1977) model presented in Barnard et al. (2016)
    ! Barnard and Bauerle. Seasonal variation in canopy aerodynamics and the sensitivity of transpiration
    ! estimates to wind velocity in broadleaved deciduous species. Journal of Hydrometeorology 2016.
    ! Note that the vertically uniform LAI distribution is assumed here.
    !
    ! !USES:
    use data_structure, only : Max_hgt, PI
    !
    ! !ARGUMENTS:
    implicit none
    integer, intent(in)    :: n_spe
    real(8), intent(in)    :: dleaf(:)
    real(8), intent(in)    :: STEP
    real(8), intent(in)    :: dpai_layer_sum(:)
    integer, intent(in)    :: nbot
    integer, intent(in)    :: ntop
    real(8), intent(inout) :: rwind_profile(:)
    !
    ! !LOCAL VARIABLES:
    integer :: p                             ! Species index
    integer :: i                             ! Layer index
    real(8) :: total_lai                     ! Plot-scale LAI (m2 leaf/m2 ground)
    real(8) :: mean_dleaf                    ! Mean leaf dimension (m)
    real(8) :: canopy_h                      ! Canopy height (= thickness) (m)
    real(8) :: canopy_z                      ! Height of the layer (m) <- from the bottom of the canopy
    real(8) :: lm                            ! Mixing length (m)
    real(8) :: alpha_g                       ! Goudriaan exponential coefficient
    !---------------------------------------------------------------------

    ! Plot-scale LAI (m2 leaf/m2 ground)

    total_lai = sum(dpai_layer_sum)

    ! Canopy height (= thickness) (m)

    canopy_h = real(ntop - nbot + 1) * STEP

    ! Mean leaf dimension (m)

    mean_dleaf = 0.0d0
    do p = 1, n_spe
       mean_dleaf = mean_dleaf + dleaf(p)
    end do
    mean_dleaf = mean_dleaf / real(n_spe)

    ! Mixing length (m)

    lm = (6.0d0 * (mean_dleaf**2.0d0) * canopy_h / (PI * total_lai)) ** (1.0d0 / 3.0d0)

    ! Goudriaan exponential coefficient

    alpha_g = (0.2d0 * total_lai * canopy_h / lm) ** 0.5d0

    ! Wind profile

    rwind_profile(:) = 1.0d0

    if (total_lai > 0.001d0) then
       do i = Max_hgt, 1, -1

          if (i <= ntop .and. i >= nbot) then

             canopy_z = real(i - nbot + 1) * STEP
             rwind_profile(i) = exp(alpha_g * (canopy_z / canopy_h - 1.0d0))

          elseif (i < nbot) then

             rwind_profile(i) = rwind_profile(nbot)

          end if

       end do
    end if

  end subroutine wind_profile

end module mod_meteorology
