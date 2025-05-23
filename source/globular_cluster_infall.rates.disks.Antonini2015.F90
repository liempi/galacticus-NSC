!! Copyright 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018,
!!           2019, 2020, 2021, 2022, 2023, 2024, 2025
!!    Andrew Benson <abenson@carnegiescience.edu>
!!
!! This file is part of Galacticus.
!!
!!    Galacticus is free software: you can redistribute it and/or modify
!!    it under the terms of the GNU General Public License as published by
!!    the Free Software Foundation, either version 3 of the License, or
!!    (at your option) any later version.
!!
!!    Galacticus is distributed in the hope that it will be useful,
!!    but WITHOUT ANY WARRANTY; without even the implied warranty of
!!    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!!    GNU General Public License for more details.
!!
!!    You should have received a copy of the GNU General Public License
!!    along with Galacticus.  If not, see <http://www.gnu.org/licenses/>.

  !!{
  Implementation of a globular cluster infall rate in galactic disks which computes the  rate over the disk.
  !!}

  use :: Dark_Matter_Halo_Scales, only : darkMatterHaloScale, darkMatterHaloScaleClass

  !![
  <globularClusterInfallRateDisks name="globularClusterInfallRateDisksAntonini2015">
   <description>
    A globular cluster infall rate in galactic disks which computes the rate by multiply the star formation rate of the disk by a factor. Specifically, the globular cluster formation rate is given by
   </description>
  </globularClusterInfallRateDisks>
  !!]
  type, extends(globularClusterInfallRateDisksClass) :: globularClusterInfallRateDisksAntonini2015
     !!{
     Implementation of a rate for globular cluster infall in galactic disks.
     !!}
     private
      class          (darkMatterHaloScaleClass), pointer :: darkMatterHaloScale_   => null()
     double precision                                    :: massMinimum                     , massMaximum
     integer                                             :: globularClusterStellarMassDiskID
   contains
     final     ::         globularClusterInfallDisksAntonini2015Destructor
     procedure :: rate => globularClusterInfallDisksAntonini2015Rate
  end type globularClusterInfallRateDisksAntonini2015

  interface globularClusterInfallRateDisksAntonini2015
     !!{
     Constructors for the {\normalfont \ttfamily globularClusterInfallDisksAntonini2015} infall rate in disks class.
     !!}
     module procedure globularClusterInfallDisksAntonini2015ConstructorParameters
     module procedure globularClusterInfallDisksAntonini2015ConstructorInternal
  end interface globularClusterInfallRateDisksAntonini2015
  
contains

  function globularClusterInfallDisksAntonini2015ConstructorParameters(parameters) result(self)
    !!{
    Constructor for the {\normalfont \ttfamily globularClusterInfallDisks} formation rate in disks class which takes a
    parameter set as input.
    !!}
    use :: Input_Parameters, only : inputParameter, inputParameters
    implicit none
    type            (globularClusterInfallRateDisksAntonini2015)                :: self
    type            (inputParameters                           ), intent(inout) :: parameters
    double precision                                                            :: massMinimum, massMaximum
    class           (darkMatterHaloScaleClass                  ), pointer       :: darkMatterHaloScale_

    !![
    <inputParameter>
      <name>massMinimum</name>
      <defaultValue>1.0d2</defaultValue>
      <description>Minimum mass of the globular clusters in the disk component.</description>
      <source>parameters</source>
    </inputParameter>
    <inputParameter>
      <name>massMaximum</name>
      <defaultValue>1.0d7</defaultValue>
      <description>Maximum mass of the globular clusters in the disk component.</description>
      <source>parameters</source>
    </inputParameter>
    <objectBuilder class="darkMatterHaloScale"  name="darkMatterHaloScale_"  source="parameters"/>
    !!]
    self=globularClusterInfallRateDisksAntonini2015(massMinimum,massMaximum, darkMatterHaloScale_)
    !![
    <inputParametersValidate source="parameters"/>
    <objectDestructor name="darkMatterHaloScale_" />
    !!]
    return
  end function globularClusterInfallDisksAntonini2015ConstructorParameters

  function globularClusterInfallDisksAntonini2015ConstructorInternal(massMinimum, massMaximum,darkMatterHaloScale_) result(self)
    !!{
    Internal constructor for the {\normalfont \ttfamily globularClusterInfallDisks} globular cluster infall rate in disks class.
    !!}
    implicit none
    type            (globularClusterInfallRateDisksAntonini2015)                        :: self
    double precision                                            , intent(in   )         :: massMinimum
    double precision                                            , intent(in   )         :: massMaximum
    class           (darkMatterHaloScaleClass                  ), intent(in   ), target :: darkMatterHaloScale_

    !![
    <constructorAssign variables="massMinimum, massMaximum,*darkMatterHaloScale_"/>
    !!]
    !![
    <addMetaProperty component="disk" name="globularClusterStellarMassDisk" id="self%globularClusterStellarMassDiskID" isEvolvable="yes" isCreator="no" />
    !!]
    return
  end function globularClusterInfallDisksAntonini2015ConstructorInternal

  subroutine globularClusterInfallDisksAntonini2015Destructor(self)
    !!{
    Internal constructor for the {\normalfont \ttfamily globularClusterInfallDisks} globular cluster infall rate in disks class.
    !!}
    implicit none
    type(globularClusterInfallRateDisksAntonini2015), intent(inout) :: self

    !![
    <objectDestructor name="self%darkMatterHaloScale_" />
    !!]
    return
  end subroutine globularClusterInfallDisksAntonini2015Destructor
  
  double precision function globularClusterInfallDisksAntonini2015Rate(self,node) result(rate)
    !!{
    Returns the globular infall rate (in $\mathrm{M}_\odot$ Gyr$^{-1}$) in the galactic spheroid of {\normalfont \ttfamily node}
    !!}
    use :: Galactic_Structure_Options, only : componentTypeDisk    , massTypeStellar          , massTypeGaseous
    use :: Galacticus_Nodes          , only : nodeComponentDisk    , nodeComponentDiskStandard, treeNode
    use :: Mass_Distributions        , only : massDistributionClass
    use :: Numerical_Constants_Math  , only : Pi
    use :: Numerical_Integration     , only : integrator
    implicit none
    class           (globularClusterInfallRateDisksAntonini2015), intent(inout), target  :: self
    type            (treeNode                                  ), intent(inout), target  :: node
    class           (massDistributionClass                     ), pointer                :: massDistributionStellar_        , massDistributionGaseous_
    class           (nodeComponentDisk                         ), pointer                :: disk
    double precision                                            , parameter              :: radiusInnerDimensionless=1.0d-13, radiusOuterDimensionless=10.0d0
    double precision                                                                     :: radiusDisk                      , massStellar                    , &
         &                                                                                  normalizationMassConstant       , massGlobularClusterDisk        , &
         &                                                                                  radiusInner                     , radiusOuter
    type            (integrator                                )                         :: integrator_

    ! Get the disk properties.
    disk        => node%disk       ()
    massStellar =  disk%massStellar()
    radiusDisk  =  disk%radius     ()

    select type(disk)
      class default 
        !Generic type, do nothing.
        rate= 0.0d0
      class is (nodeComponentDiskStandard)
        ! Check if it is a physical disk.
        if (massStellar <= 0.0d0 .or. radiusDisk <= 0.0d0) then
          ! It is not, so return zero rate.
          rate =+0.0d0
        else
          ! Here we use equation 10 from from F. Antonini, E. Barausse & J. Silk (2015; https://ui.adsabs.harvard.edu/abs/2015ApJ...812...72A).
          ! Specifically, we split the equations in two parts, the disk and spheroidal components. This allow us to track the mass of the globular
            ! clusters in each component. The integral over the globular cluster mass can be evaluated analyticaly. 

          normalizationMassConstant      = (self%massMaximum * self%massMinimum)         &
            &                             /(self%massMaximum - self%massMinimum)         ! M☉

          ! Compute suitable limits for the integration.
          radiusInner=radiusDisk*radiusInnerDimensionless
          radiusOuter=radiusDisk*radiusOuterDimensionless

          massGlobularClusterDisk = disk%floatRank0MetaPropertyGet(self%globularClusterStellarMassDiskID) ! M☉

          ! Find the rate of globSular cluster formation in the disk component.
          if (massGlobularClusterDisk <= 0.0d0) then
            rate   = 0.0d0
            return
          else
            massDistributionStellar_ => node%massDistribution(componentType=componentTypeDisk,massType=massTypeStellar)
            massDistributionGaseous_ => node%massDistribution(componentType=componentTypeDisk,massType=massTypeGaseous)

            integrator_ =  integrator(radialIntegrand,toleranceRelative=1.0d-3, hasSingularities=.true.)
            rate   = normalizationMassConstant                       &
             &      *massGlobularClusterDisk                         &
             &      /massStellar                                     &
             &      *integrator_%integrate(radiusInner,radiusOuter)  &
             &      *log(self%massMaximum/self%massMinimum)
             !![
              <objectDestructor name="massDistributionStellar_"/>
              <objectDestructor name="massDistributionGaseous_"/>
            !!]                                                
          end if
        end if
      end select 
    return

    contains
      double precision function radialIntegrand(radius)
        use :: Coordinates                     , only : coordinateCylindrical         , assignment(=)
        use :: Numerical_Constants_Astronomical, only : gravitationalConstant_internal, MpcPerKmPerSToGyr

        implicit none
        double precision                      , intent(in  ) :: radius
        type            (coordinateCylindrical)              :: coordinates
        double precision                                     :: surfaceDensityStellar            , surfaceDensityGas            , &
            &                                                   velocityDispersionStellar        , velocityDispersionGas        , &
            &                                                   rotationCurveMaximum             , virialVelocity               , &
            &                                                   rotationCurveDisk

        coordinates             = [radius,0.0d0,0.0d0]
        surfaceDensityStellar   = massDistributionStellar_ %surfaceDensity      (coordinates) ! M☉ Mpc⁻²
        surfaceDensityGas       = massDistributionGaseous_ %surfaceDensity      (coordinates) ! M☉ Mpc⁻²
        rotationCurveDisk       = massDistributionStellar_ %RotationCurve       (radius     ) ! km s⁻¹
        virialVelocity          = self%darkMatterHaloScale_%velocityVirial     (node       ) ! km s⁻¹
        rotationCurveMaximum    = massDistributionStellar_ %velocityRotationCurveMaximum(   ) ! km s⁻¹

        ! Use Eq. to approximate the velocity dispersion of stars of the gas () M. Kregel, P. C. van der Kruit, and K. C. Freeman. Structure and kinematics of edge- on galaxy discs - V. The dynamics of stellar discs. MNRAS, 358(2):503–520, Apr. 2005. doi: 10.1111/j.1365-2966.2005.08855.x.
        velocityDispersionStellar = 0.29d0*rotationCurveMaximum      ! km s⁻¹
        ! Use Eq. in  A. Dutton and F. C. van den Bosch. The impact of feedback on disc galaxy scaling relations. MNRAS, 396(1):141–164, June 2009. doi: 10.1111/j.1365-2966.2009. 14742.x.
        velocityDispersionGas     = 0.10d0*velocityDispersionStellar ! km s⁻¹

        ! Get stellar surface density.
        radialIntegrand= 2*Pi*radius*surfaceDensityStellar*(gravitationalConstant_internal**2.0d0/virialVelocity)*(surfaceDensityStellar/(velocityDispersionStellar**2.0d0*radius)+surfaceDensityGas/(velocityDispersionGas**2.0d0*radius))/MpcPerKmPerSToGyr
        return 
      end function radialIntegrand 

  end function globularClusterInfallDisksAntonini2015Rate
  