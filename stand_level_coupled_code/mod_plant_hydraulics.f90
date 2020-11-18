module mod_plant_hydraulics

  !-----------------------------------------------------------------------
  ! !DESCRIPTION:
  ! Calculate plant hydraulic resistance and maximum water uptake rate
  !
  ! !USES:
  !
  !
  ! !PUBLIC TYPES:
  implicit none
  !
  ! !PUBLIC MEMBER FUNCTIONS:
  public :: plant_hydraulics
  !-----------------------------------------------------------------------

contains

  !-----------------------------------------------------------------------
  subroutine plant_hydraulics ( &
                                       ! *** Input ***
    root_radius                   , &  ! Fine root radius (m)
    root_biomass                  , &  ! Stand-level fine root biomass (g root/tree)
    root_area                     , &  ! Stand-level root coverage area (m2 ground/tree)
    root_depth                    , &  ! Rooting depth (m)
    root_density                  , &  ! Specific root density (fine root) (g root/m3 root)
    hk                            , &  ! Soil hydraulic conductance (mmol H2O/m/s/MPa)
    root_resist                   , &  ! Hydraulic resistivity of root tissue (MPa.s.g root/mmol H2O)
    dbh_heartwood                 , &  ! Heartwood diameter (m)
    dbh_sapwood                   , &  ! Sapwood diameter (m)
    k_sap                         , &  ! Stem hydraulic conductivity (kg H2O.m/m2 sapwood/s/MPa)
    tree_h                        , &  ! Tree height (m)
    bole_h                        , &  ! Bole height (m)
    tot_psi                       , &  ! Total soil water potential (MPa)
    minlp                         , &  ! Minimum leaf water potential (MPa)
                                       !
                                       ! *** Input type ***
    universal                     , &  ! Universal variables
                                       !
                                       ! *** Input/Output ***
    flux                            &  ! Flux variables
    )
    !
    ! !DESCRIPTION:
    !
    ! !USES:
    use mod_param
    !
    ! !ARGUMENTS:
    implicit none
    real(8), intent(in) :: root_radius, root_biomass, root_area, root_depth
    real(8), intent(in) :: root_density, hk, root_resist
    real(8), intent(in) :: dbh_heartwood, dbh_sapwood, k_sap
    real(8), intent(in) :: tree_h, bole_h, tot_psi, minlp
    type(universal_type), intent(in) :: universal
    type(flux_type), intent(inout)   :: flux
    !
    ! !LOCAL VARIABLES:
    real(8) :: root_cross_sec_area          ! Cross-sectional area of fine root (m2 root)
    real(8) :: root_length_density          ! Root length density per unit soil volume (m root/m3 soil)
    real(8) :: root_dist                    ! One-half distance between roots (m)
    real(8) :: dbh                          ! Stem diameter (= DBH, m)
    real(8) :: sap_area                     ! Sapwood area (m2/tree)
    real(8) :: k_sap_molar                  ! Stem hydraulic conductivity (mmol H2O.m/m2 sapwood/s/MPa)
    real(8) :: height_psi                   ! Tree height potential (MPa)
    real(8), parameter :: pi = 3.14159265d0
    !---------------------------------------------------------------------

    associate ( &
                                      ! *** Input ***
    denh2o  => universal%denh2o  , &  ! Water density (kg/m3)
    mmh2o   => universal%mmh2o   , &  ! Molecular mass of water (kg/mol)
    grav    => universal%grav    , &  ! Gravitational acceleration (m/s2)
                                      !
                                      ! *** Output ***
    r_soil  => flux%r_soil       , &  ! Soil-to-root resistance (MPa.s.tree/mmol H2O)
    r_root  => flux%r_root       , &  ! Root-stem resistance (MPa.s.tree/mmol H2O)
    dr_sap  => flux%dr_sap       , &  ! Stem resistance per unit path length (MPa.s.tree/mmol H2O/m)
    r_sap   => flux%r_sap        , &  ! Whole-plant stem resistance (MPa.s.tree/mmol H2O)
    r_whole => flux%r_whole      , &  ! Whole plant resistance (MPa.s.tree/mmol H2O)
    sapmax  => flux%sapmax         &  ! Maximum water uptake rate (mol H2O/tree/s)
    )

    !---------------------------------------------------------------------
    ! Soil-root resistance (MPa.s.tree/mmol H2O)
    !---------------------------------------------------------------------

    ! Root cross-sectional area (m2 root)

    root_cross_sec_area = pi * root_radius**2.0d0

    ! Root length density per unit volume of soil (m root/m3 ground)

    root_length_density = root_biomass / (root_area * root_depth * root_density * root_cross_sec_area)

    ! One-half distance between roots (m)

    root_dist = sqrt(1.0d0 / (pi * root_length_density))

    ! Soil-to-root resistance (MPa.s.tree/mmol H2O)

    r_soil = log(root_dist/root_radius) / (2.0d0 * pi * root_length_density * root_depth * hk * root_area)

    !---------------------------------------------------------------------
    ! Root-stem resistance (MPa.s.tree/mmol H2O)
    !---------------------------------------------------------------------

    r_root = root_resist / root_biomass

    !---------------------------------------------------------------------
    ! Stem resistance (MPa.s.tree/mmol H2O)
    !---------------------------------------------------------------------

    ! Stem diameter (= DBH, m)

    dbh = dbh_heartwood + dbh_sapwood

    ! Sapwood area (m2/tree)

    sap_area = pi*(dbh/2.0d0)**2.0d0 - pi*(dbh_heartwood/2.0d0)**2.0d0

    ! Stem hydraulic conductivity
    ! (kg H2O.m/m2 sapwood/s/MPa) -> (mmol H2O.m/m2 sapwood/s/MPa)

    k_sap_molar = (k_sap / mmh2o) * 1000.0d0

    ! Stem resistance per unit path length (MPa.s.tree/mmol H2O/m)

    dr_sap = 1.0d0 / (k_sap_molar * sap_area)

    ! Whole-plant stem resistance (MPa.s.tree/mmol H2O)

    r_sap = dr_sap * ((tree_h + bole_h) * 0.5d0)

    !---------------------------------------------------------------------
    ! Whole plant resistance (MPa.s.tree/mmol H2O)
    !---------------------------------------------------------------------

    r_whole = r_soil + r_root + r_sap;

    !---------------------------------------------------------------------
    ! Tree height potential (MPa)
    !---------------------------------------------------------------------

    height_psi = denh2o * grav * tree_h * 1.e-06

    !---------------------------------------------------------------------
    ! Maximum water uptake (mol H2O/tree/s)
    ! Multiplying 1.e-03 for (mmol H2O/tree/s) -> (mol H2O/tree/s)
    !---------------------------------------------------------------------

    sapmax = ((tot_psi - height_psi - minlp) / r_whole)  * 1.e-03

    end associate
  end subroutine plant_hydraulics

end module mod_plant_hydraulics
