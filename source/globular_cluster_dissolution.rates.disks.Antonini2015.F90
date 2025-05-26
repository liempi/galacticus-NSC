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
  Implementation of a globular cluster Dissolution rate in galactic disks which computes the  rate over
  the disk.
  !!}
  !![
  <globularClusterDissolutionRateDisks name="globularClusterDissolutionRateDisksAntonini2015">
   <description>
    A globular cluster dissolution rate in galactic disks which computes the rate by multiply the star formation rate of the disk by a factor. Specifically, the globular cluster formation rate is given by
   </description>
  </globularClusterDissolutionRateDisks>
  !!]
  type, extends(globularClusterDissolutionRateDisksClass) :: globularClusterDissolutionRateDisksAntonini2015
     !!{
     Implementation of a rate for globular cluster dissolution in galactic disks.
     !!}
     private
     double precision :: massMinimumGlobularClusters     , massMaximumGlobularClusters
     integer          :: globularClusterStellarMassDiskID
   contains
     procedure :: rate => globularClusterDissolutionDisksAntonini2015Rate
  end type globularClusterDissolutionRateDisksAntonini2015

  interface globularClusterDissolutionRateDisksAntonini2015
     !!{
     Constructors for the {\normalfont \ttfamily globularClusterDissolutionDisksAntonini2015} dissolution rate in disks class.
     !!}
     module procedure globularClusterDssltnDisksAntonini2015ConstructorParameters
     module procedure globularClusterDssltnDisksAntonini2015ConstructorInternal
  end interface globularClusterDissolutionRateDisksAntonini2015
  
contains

  function globularClusterDssltnDisksAntonini2015ConstructorParameters(parameters) result(self)
    !!{
    Constructor for the {\normalfont \ttfamily globularClusterDissolutionDisks} dissolution rate in disks class which takes a
    parameter set as input.
    !!}
    use :: Input_Parameters, only : inputParameter, inputParameters
    implicit none
    type            (globularClusterDissolutionRateDisksAntonini2015)                :: self
    type            (inputParameters                                ), intent(inout) :: parameters
    double precision                                                                 :: massMinimumGlobularClusters, massMaximumGlobularClusters
  
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
    !!]
    self=globularClusterDissolutionRateDisksAntonini2015(massMinimumGlobularClusters,massMaximumGlobularClusters)
    !![
    <inputParametersValidate source="parameters"/>
    !!]
    return
  end function globularClusterDssltnDisksAntonini2015ConstructorParameters

  function globularClusterDssltnDisksAntonini2015ConstructorInternal(massMinimumGlobularClusters, massMaximumGlobularClusters) result(self)
    !!{
    Internal constructor for the {\normalfont \ttfamily globularClusterDissolutionDisks} dissolution rate in disks class.
    !!}
    implicit none
    type            (globularClusterDissolutionRateDisksAntonini2015)                :: self
    double precision                                                 , intent(in   ) :: massMinimumGlobularClusters
    double precision                                                 , intent(in   ) :: massMaximumGlobularClusters
    !![
    <constructorAssign variables="massMinimumGlobularClusters, massMaximumGlobularClusters"/>
    <addMetaProperty component="disk" name="globularClusterStellarMassDisk" id="self%globularClusterStellarMassDiskID" isEvolvable="yes" isCreator="no" />
    !!]
    return
  end function globularClusterDssltnDisksAntonini2015ConstructorInternal
   
  double precision function globularClusterDissolutionDisksAntonini2015Rate(self,node) result(rate)
    !!{
    Returns the globular dissolution rate (in $\mathrm{M}_\odot$ Gyr$^{-1}$) in the galactic disk of {\normalfont \ttfamily node}
    !!}
    use :: Galactic_Structure_Options, only : componentTypeDisk    , massTypeStellar          , massTypeGalactic
    use :: Galacticus_Nodes          , only : nodeComponentDisk    , nodeComponentDiskStandard, treeNode
    use :: Mass_Distributions        , only : massDistributionClass
    use :: Numerical_Integration     , only : integrator
    implicit none
    class           (globularClusterDissolutionRateDisksAntonini2015), intent(inout), target  :: self
    type            (treeNode                                       ), intent(inout), target  :: node
    class           (massDistributionClass                          ), pointer                :: massDistributionStellar_                  , massDistributionGalactic_
    class           (nodeComponentDisk                              ), pointer                :: disk
    double precision                                                 , parameter              :: radiusInnerDimensionless         =1.0d-13 , radiusOuterDimensionless=10.0d0
    double precision                                                 , parameter              :: rotationPeriodNormalization      =41.4d-3  ! Mpc⁻¹
    double precision                                                 , parameter              :: dissolutionTimescaleNormalization=10.0d0   ! Gyr
    double precision                                                 , parameter              :: globularClusterMassNormalization =2.0d5    ! M☉ 
    double precision                                                                          :: radiusDisk                                , massStellar                    , &
         &                                                                                       normalizationIntegral                     , radialIntegralResult           , &
         &                                                                                       massGlobularClusterDisk                   , globularClusterIntegralResult  , &
         &                                                                                       radiusInner                               , radiusOuter
    type            (integrator                                     )                         :: integrator_

    ! Get the disk properties.
    disk        => node%disk       ()
    massStellar =  disk%massStellar()
    radiusDisk  =  disk%radius     ()

    select type(disk)
      class default 
        !Generic type, do nothing.
        rate=+0.0d0
      class is (nodeComponentDiskStandard)
        if (massStellar<=0.0d0.or.radiusDisk<=0.0d0) then
          ! It is not, so return zero rate.
          rate=+0.0d0
        else
          ! Here we use equation 10 from from F. Antonini, E. Barausse & J. Silk (2015; https://ui.adsabs.harvard.edu/abs/2015ApJ...812...72A).
          ! Specifically, we split the equations in two parts, the disk and spheroidal components. This allow us to track the mass of the globular
          ! clusters in each component. The integral over the globular cluster mass can be evaluated analyticaly. 
          !  Ṁ_dᵢₛₛᵍᶜ = M_dᵢₛₖᵍᶜ ∫ [π_dᵢₛₖᵍᶜ(r, m_cl) / t_dᵢₛₛ(m_cl)]  2πr dr dm_cl
          massGlobularClusterDisk = disk%floatRank0MetaPropertyGet(self%globularClusterStellarMassDiskID) ! M☉

          if (massGlobularClusterDisk<=0.0d0) then
            ! If there are no globular clusters, it does not make sense to compute their dissolution.
            rate=+0.0d0
          else
            ! First, compute the normalization constant in units of M☉.
            normalizationIntegral= (self%massMaximumGlobularClusters*self%massMinimumGlobularClusters) &
              &                   /(self%massMaximumGlobularClusters-self%massMinimumGlobularClusters)   ! M☉
            ! Evaluate the analytic part of the integral for the mass of the globular clusters.
            ! ∫  m_cl⁻⁸/³ d m_cl = -3/5 m_cl⁻⁵/³+C
            globularClusterIntegralResult=-3.0d0/5.0d0                                        &
             &                            *(                                                  &
             &                              +self%massMaximumGlobularClusters**(-5.0d0/3.0d0) &
             &                              -self%massMinimumGlobularClusters**(-5.0d0/3.0d0) &
             &                             ) ! M☉⁻⁵/³
            ! Compute suitable limits for the integration.
            radiusInner=radiusDisk*radiusInnerDimensionless
            radiusOuter=radiusDisk*radiusOuterDimensionless
            ! Get the mass distributions to use in the radial integrand function.
            massDistributionStellar_ => node%massDistribution(componentType=componentTypeDisk,massType=massTypeStellar )
            massDistributionGalactic_=> node%massDistribution(                                massType=massTypeGalactic)
            ! Integrate over the radius.
            integrator_         =integrator(radialIntegrand,toleranceRelative=1.0d-3, hasSingularities=.true.)
            radialIntegralResult=integrator_%integrate(radiusInner,radiusOuter) ! M☉ Mpc⁻¹
            ! The rate is given in units of M☉ Gyr⁻¹.
            rate=+normalizationIntegral                           & ! M☉
             &   *massGlobularClusterDisk                         & ! M☉
             &   /massStellar                                     & ! M☉⁻¹
             &   *radialIntegralResult                            & ! M☉ Mpc⁻¹
             &   *globularClusterMassNormalization**(2.0d0/3.0d0) & ! M☉²/³                         
             &   /dissolutionTimescaleNormalization               & ! Gyr⁻¹
             &   /rotationPeriodNormalization                     & ! Mpc
             &   *globularClusterIntegralResult                     ! M☉⁻⁵/³  
            !![
              <objectDestructor name="massDistributionStellar_"/>
              <objectDestructor name="massDistributionGalactic_"/>
            !!]                                                
          end if
        end if
      end select 
    return

    contains
      double precision function radialIntegrand(radius)
        use :: Coordinates             , only : coordinateCylindrical, assignment(=)
        use :: Numerical_Constants_Math, only : Pi
        implicit none
        double precision                       , intent(in  ) :: radius
        type            (coordinateCylindrical)               :: coordinates
        double precision                                      :: surfaceDensity             , velocityRotation
        double precision                       , parameter    :: velocityNormalization=1.0d0 !km s⁻¹

        ! Define coordinates.
        coordinates    = [radius,0.0d0,0.0d0]
        ! Get the galactic rotation curve at the radius.
        velocityRotation=massDistributionGalactic_%rotationCurve (     radius) ! km s⁻¹
        surfaceDensity  =massDistributionStellar_ %surfaceDensity(coordinates) ! M☉ Mpc⁻²
        ! The result of the integral is M☉ Mpc⁻¹.
        radialIntegrand =+2.0d0*Pi                                 & ! Adimensional
          &              *surfaceDensity                           & ! M☉ Mpc⁻²
          &              *(velocityRotation/velocityNormalization) & ! Adimensional
        return 
      end function radialIntegrand 

  end function globularClusterDissolutionDisksAntonini2015Rate
  

  
