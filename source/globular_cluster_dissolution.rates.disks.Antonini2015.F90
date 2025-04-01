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
  use :: Dark_Matter_Halo_Scales, only : darkMatterHaloScale, darkMatterHaloScaleClass

  !![
  <globularClusterDissolutionRateDisks name="globularClusterDissolutionRateDisksAntonini2015">
   <description>
    A globular cluster dissolution rate in galactic disks which computes the rate by multiply the star formation rate of the disk by a factor. Specifically, the globular cluster formation rate is given by
   </description>
  </globularClusterDissolutionRateDisks>
  !!]
  type, extends(globularClusterDissolutionRateDisksClass) :: globularClusterDissolutionRateDisksAntonini2015
     !!{
     Implementation of a rate for globular cluster Dissolution in galactic disks.
     !!}
     private
     class           (darkMatterHaloScaleClass ), pointer :: darkMatterHaloScale_ => null()
     double precision                                     :: massMinimum                     , massMaximum
     integer                                              :: globularClusterStellarMassDiskID
   contains
     final     ::         globularClusterDissolutionDisksAntonini2015Destructor
     procedure :: rate => globularClusterDissolutionDisksAntonini2015Rate
  end type globularClusterDissolutionRateDisksAntonini2015

  interface globularClusterDissolutionRateDisksAntonini2015
     !!{
     Constructors for the {\normalfont \ttfamily globularClusterDissolutionDisksAntonini2015} Dissolution rate in disks class.
     !!}
     module procedure globularClusterDssltnSpheAntonini2015ConstructorParameters
     module procedure globularClusterDssltnSpheAntonini2015ConstructorInternal
  end interface globularClusterDissolutionRateDisksAntonini2015
  
contains

  function globularClusterDssltnSpheAntonini2015ConstructorParameters(parameters) result(self)
    !!{
    Constructor for the {\normalfont \ttfamily globularClusterDissolutionDisks} formation rate in disks class which takes a
    parameter set as input.
    !!}
    use :: Input_Parameters, only : inputParameter, inputParameters
    implicit none
    type            (globularClusterDissolutionRateDisksAntonini2015)                :: self
    type            (inputParameters                                    ), intent(inout) :: parameters
    class           (darkMatterHaloScaleClass                           ), pointer       :: darkMatterHaloScale_
    double precision                                                                     :: massMinimum         , massMaximum
  
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
    <objectBuilder class="darkMatterHaloScale" name="darkMatterHaloScale_" source="parameters"/>
    !!]
    self=globularClusterDissolutionRateDisksAntonini2015(massMinimum,massMaximum, darkMatterHaloScale_)
    !![
    <inputParametersValidate source="parameters"/>
    <objectDestructor name="darkMatterHaloScale_"/>
    !!]
    return
  end function globularClusterDssltnSpheAntonini2015ConstructorParameters

  function globularClusterDssltnSpheAntonini2015ConstructorInternal(massMinimum, massMaximum,darkMatterHaloScale_) result(self)
    !!{
    Internal constructor for the {\normalfont \ttfamily globularClusterDissolutionDisks} globular cluster Dissolution rate in disks class.
    !!}
    implicit none
    type            (globularClusterDissolutionRateDisksAntonini2015)                        :: self
    double precision                                                     , intent(in   )         :: massMinimum
    double precision                                                     , intent(in   )         :: massMaximum
    class           (darkMatterHaloScaleClass                           ), intent(in   ), target :: darkMatterHaloScale_

    !![
    <constructorAssign variables="massMinimum, massMaximum, *darkMatterHaloScale_"/>
    !!]
    !![
    <addMetaProperty component="disk" name="globularClusterStellarMassDisk" id="self%globularClusterStellarMassDiskID" isEvolvable="yes" isCreator="no" />
    !!]
    return
  end function globularClusterDssltnSpheAntonini2015ConstructorInternal
  
  subroutine globularClusterDissolutionDisksAntonini2015Destructor(self)
    !!{
    Destructor for the cut off cooling rate class.
    !!}
    implicit none
    type(globularClusterDissolutionRateDisksAntonini2015), intent(inout) :: self

    !![
    <objectDestructor name="self%darkMatterHaloScale_"/>
    !!]
    return
  end subroutine globularClusterDissolutionDisksAntonini2015Destructor
 
  double precision function globularClusterDissolutionDisksAntonini2015Rate(self,node) result(rate)
    !!{
    Returns the globular formation rate (in $\mathrm{M}_\odot$ Gyr$^{-1}$) in the galactic disk of {\normalfont \ttfamily node}
    !!}
    use :: Galactic_Structure_Options, only : componentTypeDisk    , massTypeStellar
    use :: Galacticus_Nodes          , only : nodeComponentDisk    , nodeComponentDiskStandard, treeNode
    use :: Mass_Distributions        , only : massDistributionClass
    use :: Numerical_Constants_Math  , only : Pi
    use :: Numerical_Integration     , only : integrator
    implicit none
    class           (globularClusterDissolutionRateDisksAntonini2015), intent(inout), target  :: self
    type            (treeNode                                       ), intent(inout), target  :: node
    class           (massDistributionClass                          ), pointer                :: massDistributionStellar_
    class           (nodeComponentDisk                              ), pointer                :: disk
    double precision                                                 , parameter              :: radiusInnerDimensionless=1.0d-13, radiusOuterDimensionless=10.0d0
    double precision                                                                          :: radiusDisk                      , massStellar                    , &
         &                                                                                       normalizationMass               , velocityVirial                 , &
         &                                                                                       massGlobularClusterDisk         , radiusInner                    , &
         &                                                                                       radiusOuter
    type            (integrator                                     )                         :: integrator_

    ! Get the disk properties.
    disk        => node%disk       ()
    massStellar =  disk%massStellar()
    radiusDisk  =  disk%radius     ()

    select type(disk)
      class default 
        !Generic type, do nothing.
        rate= 0.0d0
      class is (nodeComponentDiskStandard)
        if (massStellar <= 0.0d0 .or. radiusDisk <= 0.0d0) then
          ! It is not, so return zero rate.
          rate =+0.0d0
        else
          ! Here we use equation 10 from from F. Antonini, E. Barausse & J. Silk (2015; https://ui.adsabs.harvard.edu/abs/2015ApJ...812...72A).
          ! Specifically, we split the equations in two parts, the disk and spheroidal components. This allow us to track the mass of the globular
            ! clusters in each component. The integral over the globular cluster mass can be evaluated analyticaly. 

          velocityVirial                 =self%darkMatterHaloScale_%velocityVirial(node)
          normalizationMass              = (self%massMaximum * self%massMinimum)         &
            &                             /(self%massMaximum - self%massMinimum)         ! M☉
              ! Compute suitable limits for the integration.
          radiusInner=radiusDisk*radiusInnerDimensionless
          radiusOuter=radiusDisk*radiusOuterDimensionless
          massGlobularClusterDisk = disk%floatRank0MetaPropertyGet(self%globularClusterStellarMassDiskID) ! M☉

          ! Find the rate of globSular cluster formation in the disk component.
          if (massGlobularClusterDisk <= 0.0d0) then
            rate   = 0.0d0
          else
            massDistributionStellar_ => node%massDistribution(componentType=componentTypeDisk,massType=massTypeStellar)
            integrator_              =  integrator(radialIntegrand,toleranceRelative=1.0d-3, hasSingularities=.true.)
            rate   = 4.0d0*Pi                                        &      
             &      *normalizationMass                               &
             &      *massGlobularClusterDisk                         &
             &      /massStellar                                     &
             &      *(2.0d+5)**(-2.0d0/3.0d0)                        &
             &      *(10.0d0)**-1.0d0                                &
             &      *(41.4d3)**-1.0d0                                &
             &      *velocityVirial                                  &
             &      *integrator_%integrate(radiusInner,radiusOuter)  &
             &      *(-3.0d0/5.0d0                                   &
             &       *(+1.0d0/self%massMaximum**(5.0d0/3.0d0)        &
             &         -1.0d0/self%massMinimum**(5.0d0/3.0d0)        &
             &          )                                            &
             &         )
             !![
              <objectDestructor name="massDistributionStellar_"/>
            !!]                                                
          end if
        end if
      end select 
    return

    contains
      double precision function radialIntegrand(radius)
        use :: Coordinates, only : coordinateCylindrical, assignment(=)
        implicit none
        double precision                      , intent(in  ) :: radius
        type            (coordinateCylindrical)              :: coordinates
        double precision                                     :: surfaceDensity

        coordinates    = [radius,0.0d0,0.0d0]
        surfaceDensity = massDistributionStellar_%surfaceDensity(coordinates)
        ! Get stellar surface density.
        radialIntegrand  = surfaceDensity*radius
        return 
      end function radialIntegrand 

  end function globularClusterDissolutionDisksAntonini2015Rate
  

  
