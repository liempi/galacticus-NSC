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
  Implementation of a globular cluster dissolution rate in galactic spheroids.
  !!}

  !![
  <globularClusterDissolutionRateSpheroids name="globularClusterDissolutionRateSpheroidsAntonini2015">
   <description>
    A globular cluster dissolution rate in galactic spheroids which computes the rate by multiply the star formation rate of the spheroid by a factor. Specifically, the globular cluster formation rate is given by
   </description>
  </globularClusterDissolutionRateSpheroids>
  !!]
  type, extends(globularClusterDissolutionRateSpheroidsClass) :: globularClusterDissolutionRateSpheroidsAntonini2015
     !!{
     Implementation of a rate for globular cluster Dissolution in galactic spheroids.
     !!}
     private
     double precision :: massMinimum                         , massMaximum
     integer          :: globularClusterStellarMassSpheroidID
   contains
     procedure :: rate => globularClusterDissolutionSpheroidsAntonini2015Rate
  end type globularClusterDissolutionRateSpheroidsAntonini2015

  interface globularClusterDissolutionRateSpheroidsAntonini2015
     !!{
     Constructors for the {\normalfont \ttfamily globularClusterDissolutionSpheroidsAntonini2015} Dissolution rate in spheroids class.
     !!}
     module procedure globularClusterDssltnSpheAntonini2015ConstructorParameters
     module procedure globularClusterDssltnSpheAntonini2015ConstructorInternal
  end interface globularClusterDissolutionRateSpheroidsAntonini2015
  
contains

  function globularClusterDssltnSpheAntonini2015ConstructorParameters(parameters) result(self)
    !!{
    Constructor for the {\normalfont \ttfamily globularClusterDissolutionSpheroids} rate in spheroids class which takes a
    parameter set as input.
    !!}
    use :: Input_Parameters, only : inputParameter, inputParameters
    implicit none
    type            (globularClusterDissolutionRateSpheroidsAntonini2015)                :: self
    type            (inputParameters                                    ), intent(inout) :: parameters
    double precision                                                                     :: massMinimum, massMaximum
  
    !![
    <inputParameter>
      <name>massMinimum</name>
      <defaultValue>1.0d2</defaultValue>
      <description>Minimum mass of the globular clusters in the spheroid component.</description>
      <source>parameters</source>
    </inputParameter>
    <inputParameter>
      <name>massMaximum</name>
      <defaultValue>1.0d7</defaultValue>
      <description>Maximum mass of the globular clusters in the spheroid component.</description>
      <source>parameters</source>
    </inputParameter>
    !!]
    self=globularClusterDissolutionRateSpheroidsAntonini2015(massMinimum,massMaximum)
    !![
    <inputParametersValidate source="parameters"/>
    !!]
    return
  end function globularClusterDssltnSpheAntonini2015ConstructorParameters

  function globularClusterDssltnSpheAntonini2015ConstructorInternal(massMinimum, massMaximum) result(self)
    !!{
    Internal constructor for the {\normalfont \ttfamily globularClusterDissolutionSpheroids} globular cluster Dissolution rate in spheroids class.
    !!}
    implicit none
    type            (globularClusterDissolutionRateSpheroidsAntonini2015)                :: self
    double precision                                                     , intent(in   ) :: massMinimum
    double precision                                                     , intent(in   ) :: massMaximum

    !![
    <constructorAssign variables="massMinimum, massMaximum"/>
    <addMetaProperty component="spheroid" name="globularClusterStellarMassSpheroid" id="self%globularClusterStellarMassSpheroidID" isEvolvable="yes" isCreator="no" />
    !!]
    return
  end function globularClusterDssltnSpheAntonini2015ConstructorInternal
  
  double precision function globularClusterDissolutionSpheroidsAntonini2015Rate(self,node) result(rate)
    !!{
    Returns the globular dissolution rate (in $\mathrm{M}_\odot$ Gyr$^{-1}$) in the galactic spheroid of {\normalfont \ttfamily node}
    !!}
    use :: Galactic_Structure_Options, only : componentTypeSpheroid, massTypeStellar              , massTypeGalactic
    use :: Galacticus_Nodes          , only : nodeComponentSpheroid, nodeComponentSpheroidStandard, treeNode
    use :: Mass_Distributions        , only : massDistributionClass
    use :: Numerical_Integration     , only : integrator
    implicit none
    class           (globularClusterDissolutionRateSpheroidsAntonini2015), intent(inout), target  :: self
    type            (treeNode                                           ), intent(inout), target  :: node
    class           (massDistributionClass                              ), pointer                :: massDistributionStellar_        , massDistributionGalactic_
    class           (nodeComponentSpheroid                              ), pointer                :: spheroid 
    double precision                                                     , parameter              :: radiusInnerDimensionless=1.0d-13, radiusOuterDimensionless=10.0d0
    double precision                                                                              :: radiusSpheroid                  , massStellar                    , &
         &                                                                                           normalizationIntegral           , integralEvaluated              , &
         &                                                                                           massGlobularClusterSpheroid     , radiusInner                    , &
         &                                                                                           radiusOuter
    double precision                                                     , parameter              :: rotationPeriodNormalization      =41.4d-3  ! Mpc⁻¹
    double precision                                                     , parameter              :: dissolutionTimescaleNormalization=10.0d0   ! Gyr
    double precision                                                     , parameter              :: globularClusterMassNormalization =2.0d5    ! M☉ 
    type            (integrator                                         )                         :: integrator_

    ! Get the spheroid properties.
    spheroid       => node    %spheroid   ()
    massStellar    =  spheroid%massStellar()
    radiusSpheroid =  spheroid%radius     ()

    select type(spheroid)
      class default 
        !Generic type, do nothing.
        rate= 0.0d0
      class is (nodeComponentSpheroidStandard)
        if (massStellar <= 0.0d0 .or. radiusSpheroid <= 0.0d0) then
          ! It is not, so return zero rate.
          rate =+0.0d0
          return
        else
          ! Here we use equation 10 from from F. Antonini, E. Barausse & J. Silk (2015; https://ui.adsabs.harvard.edu/abs/2015ApJ...812...72A).
          ! Specifically, we split the equations in two parts, the disk and spheroidal components. This allow us to track the mass of the globular
            ! clusters in each component. The integral over the globular cluster mass can be evaluated analyticaly. 
          massGlobularClusterSpheroid = spheroid%floatRank0MetaPropertyGet(self%globularClusterStellarMassSpheroidID) ! M☉

          ! Find the rate of globSular cluster formation in the spheroid component.
          if (massGlobularClusterSpheroid <= 0.0d0) then
            rate   = 0.0d0
            return
          else

            normalizationIntegral = (self%massMaximum*self%massMinimum) &
               &                   /(self%massMaximum-self%massMinimum) ! M☉
            
            ! Compute suitable limits for the integration.
            radiusInner=radiusSpheroid*radiusInnerDimensionless
            radiusOuter=radiusSpheroid*radiusOuterDimensionless

            massDistributionStellar_ => node%massDistribution(componentType=componentTypeSpheroid,massType=massTypeStellar )
            massDistributionGalactic_=> node%massDistribution(                                    massType=massTypeGalactic)

            integrator_      = integrator(radialIntegrand,toleranceRelative=1.0d-3, hasSingularities=.true.)
            integralEvaluated= integrator_%integrate(radiusInner,radiusOuter)
            rate   =+normalizationIntegral                           &
             &      *massGlobularClusterSpheroid                     &
             &      /massStellar                                     &
             &      *integralEvaluated                               &
             &      *globularClusterMassNormalization**(2.0d0/3.0d0) &                         
             &      *dissolutionTimescaleNormalization**(-1.0d0)     &
             &      *rotationPeriodNormalization**(-1.0d0)           &
             &      *(-3.0d0/5.0d0                                   &
             &      *(+1.0d0/self%massMaximum**(5.0d0/3.0d0)         &
             &        -1.0d0/self%massMinimum**(5.0d0/3.0d0)         &
             &         )                                             &
             &        )
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
        use :: Coordinates             , only : coordinateSpherical, assignment(=)
        use :: Numerical_Constants_Math, only : Pi
        implicit none
        double precision                     , intent(in  ) :: radius
        type            (coordinateSpherical)               :: coordinates
        double precision                                    :: density                    , velocityRotation
        double precision                     , parameter    :: velocityNormalization=1.0d0 !km s⁻¹
        ! Define coordinates.
        coordinates    = [radius,0.0d0,0.0d0]
        ! Get the galactic rotation curve at the radius.
        velocityRotation=massDistributionGalactic_%rotationCurve (     radius)
        density         =massDistributionStellar_ %density       (coordinates)
        radialIntegrand =4.0d0*Pi*density*(velocityRotation/velocityNormalization)*radius**3.0d0
        return 
      end function radialIntegrand 

  end function globularClusterDissolutionSpheroidsAntonini2015Rate
  

  
