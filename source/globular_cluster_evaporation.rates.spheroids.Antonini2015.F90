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
  Implementation of a globular cluster evaporation rate in galactic spheroids which computes the  rate over
  the spheroid.
  !!}

  !![
  <globularClusterEvaporationRateSpheroids name="globularClusterEvaporationRateSpheroidsAntonini2015">
   <description>
    A globular cluster evaporation rate in galactic spheroids which computes the rate by multiply the star formation rate of the spheroid by a factor. Specifically, the globular cluster formation rate is given by
    \begin{equation}
     \dot{M}_\mathrm{disk}^\mathrm{gc} = f_\mathrm{gc} \dot{\Sigma}_\star(r) \mathrm{d}r,
    \end{equation}
    where $\dot{\Sigma}_\star(r)$ is the surface density of star formation rate.
   </description>
  </globularClusterEvaporationRateSpheroids>
  !!]
  type, extends(globularClusterEvaporationRateSpheroidsClass) :: globularClusterEvaporationRateSpheroidsAntonini2015
     !!{
     Implementation of a rate for globular cluster evaporation in galactic spheroids.
     !!}
     private
     double precision :: massMinimum                         , massMaximum                               
     integer          :: globularClusterStellarMassSpheroidID
   contains
     procedure        :: rate => globularClusterEvaporationSpheroidsAntonini2015Rate
  end type globularClusterEvaporationRateSpheroidsAntonini2015

  interface globularClusterEvaporationRateSpheroidsAntonini2015
     !!{
     Constructors for the {\normalfont \ttfamily globularClusterEvaporationSpheroidsAntonini2015} evaporation rate in spheroids class.
     !!}
     module procedure globularClusterEvaprtnSpheAntonini2015ConstructorParameters
     module procedure globularClusterEvaprtnSpheAntonini2015ConstructorInternal
  end interface globularClusterEvaporationRateSpheroidsAntonini2015
    
contains

  function globularClusterEvaprtnSpheAntonini2015ConstructorParameters(parameters) result(self)
    !!{
    Constructor for the {\normalfont \ttfamily globularClusterEvaporationSpheroids} formation rate in spheroids class which takes a
    parameter set as input.
    !!}
    use :: Input_Parameters, only : inputParameter, inputParameters
    implicit none
    type            (globularClusterEvaporationRateSpheroidsAntonini2015)                :: self
    type            (inputParameters                                    ), intent(inout) :: parameters
    double precision                                                                     :: massMinimum, massMaximum
  
    !![
    <inputParameter>
      <name>massMinimum</name>
      <defaultValue>1.0d2</defaultValue>
      <description>Minimum mass of the globular clusters in the spheroids component.</description>
      <source>parameters</source>
    </inputParameter>
    <inputParameter>
      <name>massMaximum</name>
      <defaultValue>1.0d7</defaultValue>
      <description>Maximum mass of the globular clusters in the spheroids component.</description>
      <source>parameters</source>
    </inputParameter>
    !!]
    self=globularClusterEvaporationRateSpheroidsAntonini2015(massMinimum,massMaximum)
    !![
    <inputParametersValidate source="parameters"/>
    !!]
    return
  end function globularClusterEvaprtnSpheAntonini2015ConstructorParameters

  function globularClusterEvaprtnSpheAntonini2015ConstructorInternal(massMinimum, massMaximum) result(self)
    !!{
    Internal constructor for the {\normalfont \ttfamily globularClusterEvaporationSpheroids} globular cluster evaporation rate in spheroids class.
    !!}
    implicit none
    type            (globularClusterEvaporationRateSpheroidsAntonini2015)                :: self
    double precision                                                     , intent(in   ) :: massMinimum
    double precision                                                     , intent(in   ) :: massMaximum

    !![
    <constructorAssign variables="massMinimum, massMaximum"/>
    !!]
    !![
    <addMetaProperty component="spheroid" name="globularClusterStellarMassSpheroid" id="self%globularClusterStellarMassSpheroidID" isEvolvable="yes" isCreator="no" />
    !!]
    return
  end function globularClusterEvaprtnSpheAntonini2015ConstructorInternal

  double precision function globularClusterEvaporationSpheroidsAntonini2015Rate(self,node) result(rate)
    !!{
    Returns the globular formation rate (in $\mathrm{M}_\odot$ Gyr$^{-1}$) in the galactic spheroid of {\normalfont \ttfamily node}
    !!}
    use :: Galactic_Structure_Options, only : componentTypeSpheroid, massTypeStellar
    use :: Galacticus_Nodes          , only : nodeComponentSpheroid, treeNode
    use :: Mass_Distributions        , only : massDistributionClass
    use :: Numerical_Constants_Math  , only : Pi
    use :: Numerical_Integration     , only : integrator
    
    implicit none
    class           (globularClusterEvaporationRateSpheroidsAntonini2015), intent(inout), target  :: self
    type            (treeNode                                           ), intent(inout), target  :: node
    class           (massDistributionClass                              ), pointer                :: massDistributionStellar_
    double precision                                                     , parameter              :: radiusInnerDimensionless=1.0d-13, radiusOuterDimensionless=1.0d0
    class           (nodeComponentSpheroid                              ), pointer                :: spheroid 
    double precision                                                                              :: radiusSpheroid             , massStellar         , &
         &                                                                                           normalizationMass          , evaporationTimescale, &
         &                                                                                           massGlobularClusterSpheroid, radiusInner         , &
         &                                                                                           radiusOuter
    type            (integrator                                         )                         :: integrator_

    ! Get the spheroid properties.
    spheroid       => node    %spheroid   ()
    massStellar    =  spheroid%massStellar()
    radiusSpheroid =  spheroid%radius     ()
    ! Check if the spheroid is physical.
    if (massStellar <= 0.0d0 .or. radiusSpheroid <= 0.0d0) then
        ! It is not, so return zero rate.
        rate=+0.0d0
    else
        ! Here we use equation 15 from from F. Antonini, E. Barausse & J. Silk (2015; https://ui.adsabs.harvard.edu/abs/2015ApJ...812...72A).
        ! Specifically, we split the equations in two parts, the disk and spheroidal components. This allow us to track the mass of the globular
        ! clusters in each component. The integral, can be evaluated analyticaly. 

        normalizationMass           = (self%massMaximum * self%massMinimum)  &
          &                          /(self%massMaximum - self%massMinimum)  ! M☉
        evaporationTimescale        = 17.0d0 / 2.0d5 ! Gyr M☉⁻¹
        massGlobularClusterSpheroid = spheroid%floatRank0MetaPropertyGet(self%globularClusterStellarMassSpheroidID    ) ! M☉

        ! Find the rate of globular cluster formation in the spheroid component.
        if (massGlobularClusterSpheroid <= 0.0d0) then
          rate   = 0.0d0
          return
        else
          radiusInner = radiusSpheroid*radiusInnerDimensionless
          radiusOuter = radiusSpheroid*radiusOuterDimensionless

            massDistributionStellar_ => node%massDistribution(componentType=componentTypeSpheroid,massType=massTypeStellar)
            integrator_              =  integrator(radialIntegrand,toleranceRelative=1.0d-3, hasSingularities=.true.)

            rate   =+4.0d0*Pi                    &                    
             &      *normalizationMass           &
             &      *massGlobularClusterSpheroid &
             &      *(-1.0d0/2.0d0            &
             &       *(1.0d0/                 &
             &         self%massMaximum**2.0d0&
             &        -1.0d0/self%massMinimum**2.0d0&
             &          )                     &
             &         )                      &
             &       *integrator_%integrate(radiusInner,radiusOuter)  &
             &      /evaporationTimescale      ! M☉ Gyr⁻¹
            !![
              <objectDestructor name="massDistributionStellar_"/>
            !!]     
        end if 
      end if
    return

    contains
      double precision function radialIntegrand(radius)
        use :: Coordinates, only : coordinateSpherical, assignment(=)
        implicit none
        double precision                     , intent(in  ) :: radius
        type            (coordinateSpherical)               :: coordinates
        double precision                                    :: density

        coordinates = [radius,0.0d0,0.0d0]
        density     = massDistributionStellar_%density(coordinates)
        ! Get stellar surface density.
        radialIntegrand =density*radius**2.0d0
        return 
      end function radialIntegrand 
  end function globularClusterEvaporationSpheroidsAntonini2015Rate

  
