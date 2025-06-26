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
  Implementation of a globular cluster infall rate in galactic disks which computes the rate over the disk.
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
     double precision                                    :: massMinimumGlobularClusters     , massMaximumGlobularClusters
     integer                                             :: globularClusterStellarMassDiskID, globularClusterInfallTimescaleDiskID
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
    double precision                                                            :: massMinimumGlobularClusters, massMaximumGlobularClusters
    class           (darkMatterHaloScaleClass                  ), pointer       :: darkMatterHaloScale_

    !![
    <inputParameter>
      <name>massMinimumGlobularClusters</name>
      <defaultValue>1.0d2</defaultValue>
      <description>Minimum mass of the globular clusters in the disk component.</description>
      <source>parameters</source>
    </inputParameter>
    <inputParameter>
      <name>massMaximumGlobularClusters</name>
      <defaultValue>1.0d7</defaultValue>
      <description>Maximum mass of the globular clusters in the disk component.</description>
      <source>parameters</source>
    </inputParameter>
    <objectBuilder class="darkMatterHaloScale"  name="darkMatterHaloScale_"  source="parameters"/>
    !!]
    self=globularClusterInfallRateDisksAntonini2015(massMinimumGlobularClusters,massMaximumGlobularClusters, darkMatterHaloScale_)
    !![
    <inputParametersValidate source="parameters"/>
    <objectDestructor name="darkMatterHaloScale_" />
    !!]
    return
  end function globularClusterInfallDisksAntonini2015ConstructorParameters

  function globularClusterInfallDisksAntonini2015ConstructorInternal(massMinimumGlobularClusters,massMaximumGlobularClusters,darkMatterHaloScale_) result(self)
    !!{
    Internal constructor for the {\normalfont \ttfamily globularClusterInfallDisks} globular cluster infall rate in disks class.
    !!}
    implicit none
    type            (globularClusterInfallRateDisksAntonini2015)                        :: self
    double precision                                            , intent(in   )         :: massMinimumGlobularClusters
    double precision                                            , intent(in   )         :: massMaximumGlobularClusters
    class           (darkMatterHaloScaleClass                  ), intent(in   ), target :: darkMatterHaloScale_

    !![
    <constructorAssign variables="massMinimumGlobularClusters, massMaximumGlobularClusters,*darkMatterHaloScale_"/>
    <addMetaProperty component="disk" name="globularClusterStellarMassDisk"     id="self%globularClusterStellarMassDiskID"     isEvolvable="yes" isCreator="no" />
    <addMetaProperty component="disk" name="globularClusterInfallTimescaleDisk" id="self%globularClusterInfallTimescaleDiskID" isEvolvable="yes"  isCreator="yes"/>
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
    use :: Galactic_Structure_Options      , only : componentTypeDisk             , massTypeStellar          , massTypeGaseous
    use :: Galacticus_Nodes                , only : nodeComponentDisk             , nodeComponentDiskStandard, treeNode
    use :: Coordinates                     , only : coordinateCylindrical         , assignment(=)
    use :: Mass_Distributions              , only : massDistributionClass
    use :: Numerical_Integration           , only : integrator
    use :: Numerical_Constants_Astronomical, only : gravitationalConstant_internal, megaParsec               , gigayear       , MpcPerKmPerSToGyr
    implicit none
    class           (globularClusterInfallRateDisksAntonini2015), intent(inout), target  :: self
    type            (treeNode                                  ), intent(inout), target  :: node
    class           (massDistributionClass                     ), pointer                :: massDistributionStellar_        , massDistributionGaseous_
    class           (nodeComponentDisk                         ), pointer                :: disk
    type            (coordinateCylindrical                     )                         :: coordinates
    double precision                                            , parameter              :: radiusInnerDimensionless=1.0d-13, radiusOuterDimensionless=10.0d0
    double precision                                                                     :: radiusDisk                      , massStellar                      , &
         &                                                                                  normalizationIntegral           , massGlobularClusterDisk          , &
         &                                                                                  radiusInner                     , radiusOuter                      , &
         &                                                                                  radialIntegralResult            , globularClusterMassIntegralResult, &
         &                                                                                  virialVelocity                  , dynamicalFrictionTimescaleDisk
    type            (integrator                                )                         :: integrator_
    ! Get the disk properties.
    disk        => node%disk       ()
    massStellar =  disk%massStellar()
    radiusDisk  =  disk%radius     ()
    select type(disk)
      class default 
        !Generic type, do nothing.
        rate=+0.0d0
      class is (nodeComponentDiskStandard)
        ! Check if it is a physical disk.
        if (massStellar<=0.0d0.or.radiusDisk<=0.0d0) then
          ! It is not, so return zero rate.
          rate=+0.0d0
        else
          ! Here we use equation 10 from from F. Antonini, E. Barausse & J. Silk (2015; https://ui.adsabs.harvard.edu/abs/2015ApJ...812...72A).
          ! Specifically, we split the equations in two parts, the disk and spheroidal components. This allow us to track the mass of the globular
            ! clusters in each component. The integral over the globular cluster mass can be evaluated analyticaly. 
          massGlobularClusterDisk=disk%floatRank0MetaPropertyGet(self%globularClusterStellarMassDiskID) ! M☉
          ! Find the rate of globSular cluster formation in the disk component.
          if (massGlobularClusterDisk<= 0.0d0) then
            rate=+0.0d0
          else
            virialVelocity       = self%darkMatterHaloScale_%velocityVirial      (node       ) ! km s⁻¹
            normalizationIntegral=+(self%massMaximumGlobularClusters*self%massMinimumGlobularClusters) &
              &                   /(self%massMaximumGlobularClusters-self%massMinimumGlobularClusters)  ! M☉
            ! The integral over the mass of globular clusters can be evaluated analiticaly, let's say
            ! the mass of a globular cluster is "m", then the integral is
            ! ∫ m⁻¹dm = ln(mₘₐₓ/mₘᵢₙ)
            globularClusterMassIntegralResult=log(self%massMaximumGlobularClusters/self%massMinimumGlobularClusters)
            ! Compute suitable limits for the integration.
            radiusInner=radiusDisk*radiusInnerDimensionless
            radiusOuter=radiusDisk*radiusOuterDimensionless
            ! Get the mass distribution used by the radial integrand function.
            massDistributionStellar_ => node%massDistribution(componentType=componentTypeDisk,massType=massTypeStellar)
            massDistributionGaseous_ => node%massDistribution(componentType=componentTypeDisk,massType=massTypeGaseous)
            
            coordinates                      =[radiusDisk,0.0d0,0.0d0]

            dynamicalFrictionTimescaleDisk=+(gravitationalConstant_internal**2.0d0                                         &
              &                             /virialVelocity                                                                &
              &                             /radiusDisk                                                                    &
              &                             *(+         massDistributionStellar_%surfaceDensity              (coordinates)  &
              &                               /(0.290d0*massDistributionStellar_%velocityRotationCurveMaximum(           )) &
              &                               +         massDistributionGaseous_%surfaceDensity              (coordinates)  &
              &                               /(0.028d0*massDistributionStellar_%velocityRotationCurveMaximum(           )) &
              &                              )                                                                              &                               
              &                             /MpcPerKmPerSToGyr                                                             & ! Gyr M☉
              &                            )**(-1.0d0)                                                                     &
              &                             *(+0.5d0                                                                       &
              &                                     *(+self%massMaximumGlobularClusters**(-2.0d0)                          & ! Here we use the analytic solution
              &                                       -self%massMinimumGlobularClusters**(-2.0d0)                          & ! Units are M☉⁻¹
              &                                      )                                                                     &
              &                                     /(+self%massMaximumGlobularClusters**(-1.0d0)                          &
              &                                       -self%massMinimumGlobularClusters**(-1.0d0)                          & 
              &                                      )                                                                     &
              &                              )
            call disk%floatRank0MetaPropertySet(self%globularClusterInfallTimescaleDiskID, dynamicalFrictionTimescaleDisk)


            ! Integrate over the radius.
            integrator_         = integrator(radialIntegrand,toleranceRelative=1.0d-3, hasSingularities=.true.)
            radialIntegralResult= integrator_%integrate(radiusInner,radiusOuter)
            ! The rate is returned in units of M☉ Gyr⁻¹.
            rate=+normalizationIntegral             & ! M☉
             &   *massGlobularClusterDisk           & ! M☉
             &   /massStellar                       & ! M☉⁻¹
             &   *radialIntegralResult              & ! s⁻¹
             &   *globularClusterMassIntegralResult & ! Adimensional 
             &   *gigayear                            ! Convert from s⁻¹ to Gyr.
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
        use :: Numerical_Constants_Math        , only : Pi
        use :: Numerical_Constants_Prefixes    , only : kilo
        implicit none
        double precision                       , intent(in  ) :: radius
        double precision                                      :: surfaceDensityStellar    , surfaceDensityGas    , &
            &                                                    velocityDispersionStellar, velocityDispersionGas, &
            &                                                    rotationCurveMaximum

        coordinates             = [radius,0.0d0,0.0d0]

        ! Get stellar and gas surface density.
        surfaceDensityStellar   = massDistributionStellar_ %surfaceDensity      (coordinates) ! M☉ Mpc⁻²
        surfaceDensityGas       = massDistributionGaseous_ %surfaceDensity      (coordinates) ! M☉ Mpc⁻²
        ! Get rotation curves and virial velocity.
        rotationCurveMaximum    = massDistributionStellar_ %velocityRotationCurveMaximum(   ) ! km s⁻¹

        ! Use Eq. to approximate the velocity dispersion of stars of the gas ()
        ! M. Kregel, P. C. van der Kruit, and K. C. Freeman. Structure and kinematics of edge- on galaxy discs - V. The dynamics of stellar discs. MNRAS, 358(2):503–520, Apr. 2005. doi: 10.1111/j.1365-2966.2005.08855.x.
        velocityDispersionStellar = 0.29d0*rotationCurveMaximum      ! km s⁻¹
        ! Use Eq. in  A. Dutton and F. C. van den Bosch. The impact of feedback on disc galaxy scaling relations. MNRAS, 396(1):141–164, June 2009. doi: 10.1111/j.1365-2966.2009. 14742.x.
        velocityDispersionGas     = 0.10d0*velocityDispersionStellar ! km s⁻¹
        ! Define the integrand. The result it is originaly in km s⁻¹, but we transform from km to Mpc to be consistent with the units
        ! of dr. The final unit is s⁻¹.
        radialIntegrand=+2.0d0*Pi*radius                                                   & ! Mpc
           &            *surfaceDensityStellar                                             & ! M☉ Mpc⁻²
           &            *(gravitationalConstant_internal**2.0d0/virialVelocity)            & ! Mpc² M☉⁻²(km s⁻¹)³
           &            *(                                                                 &
           &              +surfaceDensityStellar/(velocityDispersionStellar**2.0d0*radius) & ! M☉ Mpc⁻¹ (km s⁻¹)⁻²
                          +surfaceDensityGas    /(velocityDispersionGas**2.0d0*radius    ) & ! M☉ Mpc⁻¹ (km s⁻¹)⁻²
           &             )                                                                 & 
           &            /megaParsec                                                        &
           &            *kilo
        return 
      end function radialIntegrand 

  end function globularClusterInfallDisksAntonini2015Rate
  