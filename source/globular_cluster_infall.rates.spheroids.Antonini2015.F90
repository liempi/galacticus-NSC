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
  Implementation of a globular cluster infall rate in galactic spheroids which computes the  rate over
  the spheroid.
  !!}
  use :: Dark_Matter_Halo_Scales, only : darkMatterHaloScale, darkMatterHaloScaleClass
  !![
  <globularClusterInfallRateSpheroids name="globularClusterInfallRateSpheroidsAntonini2015">
   <description>
    A globular cluster dissolution rate in galactic spheroids which computes the rate by multiply the star formation rate of the spheroid by a factor. Specifically, the globular cluster formation rate is given by
   </description>
  </globularClusterInfallRateSpheroids>
  !!]
  type, extends(globularClusterInfallRateSpheroidsClass) :: globularClusterInfallRateSpheroidsAntonini2015
     !!{
     Implementation of a rate for globular cluster Dissolution in galactic spheroids.
     !!}
     private
     class           (darkMatterHaloScaleClass), pointer :: darkMatterHaloScale_ => null()
     double precision                                    :: massMinimumGlobularClusters         , massMaximumGlobularClusters
     integer                                             :: globularClusterStellarMassSpheroidID
   contains
     final     ::         globularClusterInfallSpheroidsAntonini2015Destructor
     procedure :: rate => globularClusterInfallSpheroidsAntonini2015Rate
  end type globularClusterInfallRateSpheroidsAntonini2015

  interface globularClusterInfallRateSpheroidsAntonini2015
     !!{
     Constructors for the {\normalfont \ttfamily globularClusterInfallSpheroidsAntonini2015} Dissolution rate in spheroids class.
     !!}
     module procedure globularClusterInfallSpheAntonini2015ConstructorParameters
     module procedure globularClusterInfallSpheAntonini2015ConstructorInternal
  end interface globularClusterInfallRateSpheroidsAntonini2015
  
contains

  function globularClusterInfallSpheAntonini2015ConstructorParameters(parameters) result(self)
    !!{
    Constructor for the {\normalfont \ttfamily globularClusterDissolutionSpheroids} formation rate in spheroids class which takes a
    parameter set as input.
    !!}
    use :: Input_Parameters, only : inputParameter, inputParameters
    implicit none
    type            (globularClusterInfallRateSpheroidsAntonini2015)                :: self
    type            (inputParameters                               ), intent(inout) :: parameters
    class           (darkMatterHaloScaleClass                      ), pointer       :: darkMatterHaloScale_
    double precision                                                                :: massMinimumGlobularClusters, massMaximumGlobularClusters
  
    !![
    <inputParameter>
      <name>massMinimumGlobularClusters</name>
      <defaultValue>1.0d2</defaultValue>
      <description>Minimum mass of the globular clusters in the spheroid component.</description>
      <source>parameters</source>
    </inputParameter>
    <inputParameter>
      <name>massMaximumGlobularClusters</name>
      <defaultValue>1.0d7</defaultValue>
      <description>Maximum mass of the globular clusters in the spheroid component.</description>
      <source>parameters</source>
    </inputParameter>
    <objectBuilder class="darkMatterHaloScale" name="darkMatterHaloScale_" source="parameters"/>
    !!]
    self=globularClusterInfallRateSpheroidsAntonini2015(massMinimumGlobularClusters,massMaximumGlobularClusters,darkMatterHaloScale_)
    !![
    <inputParametersValidate source="parameters"/>
    <objectDestructor name="darkMatterHaloScale_"/>
    !!]
    return
  end function globularClusterInfallSpheAntonini2015ConstructorParameters

  function globularClusterInfallSpheAntonini2015ConstructorInternal(massMinimumGlobularClusters, massMaximumGlobularClusters,darkMatterHaloScale_) result(self)
    !!{
    Internal constructor for the {\normalfont \ttfamily globularClusterDissolutionspheroids} globular cluster Dissolution rate in spheroids class.
    !!}
    implicit none
    type            (globularClusterInfallRateSpheroidsAntonini2015)                        :: self
    double precision                                                , intent(in   )         :: massMinimumGlobularClusters
    double precision                                                , intent(in   )         :: massMaximumGlobularClusters
    class           (darkMatterHaloScaleClass                      ), intent(in   ), target :: darkMatterHaloScale_

    !![
    <constructorAssign variables="massMinimumGlobularClusters, massMaximumGlobularClusters, *darkMatterHaloScale_"/>

    <addMetaProperty component="spheroid" name="globularClusterStellarMassSpheroid" id="self%globularClusterStellarMassSpheroidID" isEvolvable="yes" isCreator="no"/>
    !!]
    return
  end function globularClusterInfallSpheAntonini2015ConstructorInternal
  
  subroutine globularClusterInfallSpheroidsAntonini2015Destructor(self)
    !!{
    Destructor for the cut off cooling rate class.
    !!}
    implicit none
    type(globularClusterInfallRateSpheroidsAntonini2015), intent(inout) :: self

    !![
    <objectDestructor name="self%darkMatterHaloScale_"/>
    !!]
    return
  end subroutine globularClusterInfallSpheroidsAntonini2015Destructor
 
  double precision function globularClusterInfallSpheroidsAntonini2015Rate(self,node) result(rate)
    !!{
    Returns the globular formation rate (in $\mathrm{M}_\odot$ Gyr$^{-1}$) in the galactic spheroid of {\normalfont \ttfamily node}
    !!}
    use :: Galactic_Structure_Options  , only : componentTypeSpheroid, massTypeStellar              , massTypeGaseous
    use :: Galacticus_Nodes            , only : nodeComponentSpheroid, nodeComponentSpheroidStandard, treeNode
    use :: Mass_Distributions          , only : massDistributionClass
    use :: Numerical_Integration       , only : integrator
    implicit none
    class           (globularClusterInfallRateSpheroidsAntonini2015), intent(inout), target  :: self
    type            (treeNode                                      ), intent(inout), target  :: node
    class           (massDistributionClass                         ), pointer                :: massDistributionStellar_                  , massDistributionGaseous_
    class           (nodeComponentSpheroid                         ), pointer                :: spheroid
    double precision                                                , parameter              :: radiusInnerDimensionless          =1.0d-13, radiusOuterDimensionless=10.0d0
    double precision                                                , parameter              :: dynamicalFrictionTimescale        =15.0d0  ! Gyr
    double precision                                                , parameter              :: dynamicalFrictionMassNormalization=1.0d7   ! M
    double precision                                                                         :: radiusSpheroid                            , massStellarSpheroid              , &
         &                                                                                      normalizationIntegral                     , velocityVirial                   , &
         &                                                                                      massGlobularClusterSpheroid               , radiusInner                      , &
         &                                                                                      radiusOuter                               , globularClusterMassIntegralResult, &
         &                                                                                      radialIntegralResult
    type            (integrator                                    )                         :: integrator_

    ! Get the spheroid properties.
    spheroid           => node    %spheroid   ()
    massStellarSpheroid=  spheroid%massStellar()
    radiusSpheroid     =  spheroid%radius     ()

    select type(spheroid)
      class default 
        !Generic type, do nothing.
        rate=+0.0d0
      class is (nodeComponentSpheroidStandard)
        if (massStellar<=0.0d0.or.radiusSpheroid<=0.0d0) then
          ! It is not, so return zero rate.
          rate=+0.0d0
        else
          ! Here we use equation 10 from from F. Antonini, E. Barausse & J. Silk (2015; https://ui.adsabs.harvard.edu/abs/2015ApJ...812...72A).
          ! Specifically, we split the equations in two parts, the spheroid and spheroid components. This allow us to track the mass of the globular
          ! clusters in each component. 
          massGlobularClusterSpheroid=spheroid%floatRank0MetaPropertyGet(self%globularClusterStellarMassSpheroidID) ! M☉
          
          if (massGlobularClusterSpheroid<=0.0d0) then
            rate=+0.0d0
          else
            velocityVirial       =+self%darkMatterHaloScale_%velocityVirial(node)
            normalizedVelocity   = 100.0d0/(0.65d0*velocityVirial)
            normalizationIntegral=+(self%massMaximumGlobularClusters*self%massMinimumGlobularClusters) &
              &                   /(self%massMaximumGlobularClusters-self%massMinimumGlobularClusters)  ! M☉
            ! The integral over the mass of globular clusters can be evaluated analiticaly, let's say
            ! the mass of a globular cluster is "m", then the integral is
            ! ∫ m⁻¹dm = ln(mₘₐₓ/mₘᵢₙ)
            globularClusterMassIntegralResult=log(self%massMaximumGlobularClusters/self%massMinimumGlobularClusters)
            ! Compute suitable limits for the integration.
            radiusInner=radiusSpheroid*radiusInnerDimensionless
            radiusOuter=radiusSpheroid*radiusOuterDimensionless
            ! Get the mass distributions to parse to the integrand function.
            massDistributionStellar_ => node%massDistribution(componentType=componentTypeSpheroid,massType=massTypeStellar)
            massDistributionGaseous_ => node%massDistribution(componentType=componentTypeSpheroid,massType=massTypeGaseous)
            ! Integrate over the radius.
            integrator_              =  integrator(radialIntegrand,toleranceRelative=1.0d-3, hasSingularities=.true.)
            radialIntegralResult     =  integrator_%integrate(radiusInner,radiusOuter)
            
            ! Compute the rate.
            rate   =+normalizationIntegral              & ! M☉
             &      *massGlobularClusterSpheroid        & ! M☉
             &      *globularClusterMassIntegralResult  & ! Adimensional
             &      *radialIntegralResult               & ! M☉
             &      *normalizedVelocity                 & ! Adimensional
             &      /massStellarSpheroid                & ! M☉⁻¹
             &      /dynamicalFrictionMassNormalization & ! M☉⁻¹
             &      /dynamicalFrictionTimescale           ! Gyr⁻¹

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
        use :: Coordinates                 , only : coordinateSpherical, assignment(=)
        use :: Numerical_Constants_Math    , only : Pi
        use :: Numerical_Constants_Prefixes, only : kilo               , mega
        implicit none
        double precision                     , intent(in  ) :: radius
        type            (coordinateSpherical)               :: coordinates
        double precision                                    :: density    , normalizationRadius

        coordinates        = [radius,0.0d0,0.0d0]
        normalizationRadius= (5*kilo/mega)
        ! Get stellar density.
        density     = massDistributionStellar_%density(coordinates)
        ! Return the integral in units of M☉.
        radialIntegrand=+4.0d0*Pi                            & ! Adimensional
          &             *density                             & ! M☉ Mpc⁻³
          &             *radius**2.0d0                       & ! Mpc²
          &             /(radius/normalizationRadius)**2.0d0   ! Adimensional
        return 
      end function radialIntegrand 

  end function globularClusterInfallSpheroidsAntonini2015Rate
  