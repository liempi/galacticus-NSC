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
  Implements a property extractor class for the density at a set of radii.
  !!}

  !![
  <nodePropertyExtractor name="nodePropertyExtractorSurfaceDensityDisk">
   <description>A property extractor class for the surface density at a set of radii.</description>
  </nodePropertyExtractor>
  !!]
  type, extends(nodePropertyExtractorArray) :: nodePropertyExtractorSurfaceDensityDisk
     !!{
     A property extractor class for the density at a set of radii.
     !!}
     private
   contains
     procedure :: elementCount       => surfaceDensityProfileElementCount
     procedure :: extract            => surfaceDensityProfileExtract
     procedure :: names              => surfaceDensityProfileNames
     procedure :: descriptions       => surfaceDensityProfileDescriptions
     procedure :: columnDescriptions => surfaceDensityProfileColumnDescriptions
     procedure :: size               => surfaceDensityProfileSize
     procedure :: unitsInSI          => surfaceDensityProfileUnitsInSI
  end type nodePropertyExtractorSurfaceDensityDisk

  interface nodePropertyExtractorSurfaceDensityDisk
     !!{
     Constructors for the {\normalfont \ttfamily surfaceDensityProfile} output analysis class.
     !!}
     module procedure surfaceDensityProfileConstructorParameters
  end interface nodePropertyExtractorSurfaceDensityDisk

contains

  function surfaceDensityProfileConstructorParameters(parameters) result(self)
    !!{
    Constructor for the {\normalfont \ttfamily surfaceDensityProfile} property extractor class which takes a parameter set as input.
    !!}
    use :: Input_Parameters, only : inputParameter, inputParameters
    implicit none
    type(nodePropertyExtractorSurfaceDensityDisk)                :: self
    type(inputParameters                        ), intent(inout) :: parameters

    self=nodePropertyExtractorSurfaceDensityDisk()
    !![
    <inputParametersValidate source="parameters"/>
    !!]
    return
  end function surfaceDensityProfileConstructorParameters

  integer function surfaceDensityProfileElementCount(self,time)
    !!{
    Return the number of elements in the {\normalfont \ttfamily surfaceDensityProfile} property extractors.
    !!}
    implicit none
    class           (nodePropertyExtractorSurfaceDensityDisk), intent(inout) :: self
    double precision                                         , intent(in   ) :: time
    !$GLC attributes unused :: self, time

    surfaceDensityProfileElementCount=2
    return
  end function surfaceDensityProfileElementCount

   function surfaceDensityProfileSize(self,time)
    !!{
    Return the number of array elements in the {\normalfont \ttfamily massProfile} property extractors.
    !!}
    implicit none
    integer         (c_size_t                               )                :: surfaceDensityProfileSize
    class           (nodePropertyExtractorSurfaceDensityDisk), intent(inout) :: self
    double precision                                         , intent(in   ) :: time
    !$GLC attributes unused :: time

    surfaceDensityProfileSize=20
    return
  end function surfaceDensityProfileSize

  function surfaceDensityProfileExtract(self,node,time,instance)
    !!{
    Implement a {\normalfont \ttfamilysurfaceDensityProfile} property extractor.
    !!}
    use :: Galactic_Structure_Options, only : componentTypeDisk    , massTypeStellar
    use :: Galacticus_Nodes          , only : nodeComponentDisk    , treeNode
    use :: Mass_Distributions        , only : massDistributionClass
    use :: Numerical_Ranges          , only : Make_Range
    use :: Coordinates               , only : coordinateCylindrical, assignment(=)
    use :: Error                     , only : Error_Report
    use :: Numerical_Constants_Math  , only : Pi
    implicit none
    double precision                                         , dimension(:,:), allocatable :: surfaceDensityProfileExtract
    class           (nodePropertyExtractorSurfaceDensityDisk), intent(inout) , target      :: self
    type            (treeNode                               ), intent(inout) , target      :: node
    double precision                                         , intent(in   )               :: time
    type            (multiCounter                           ), intent(inout) , optional    :: instance
    double precision                                         , allocatable   , dimension(:):: radii
    class           (nodeComponentDisk                      ), pointer                     :: disk
    class           (massDistributionClass                  ), pointer                     :: massDistribution_
    type            (coordinateCylindrical                  )                              :: coordinates
    double precision                                                                       :: radius
    integer                                                                                :: i
    !$GLC attributes unused :: instance

    disk   => node%disk  ()
    radius =  disk%radius()
    radii  =  Make_Range(radius*1.0d-13,radius,20)
    
    allocate(surfaceDensityProfileExtract(20,2))

    select type(disk)
    type is (nodeComponentDisk)
      ! Disk does not yet exist - nothin to do here.
      do i=1, 20

        surfaceDensityProfileExtract(i,1) =  0.0d0
        surfaceDensityProfileExtract(i,2) =  0.0d0
      end do
    class default
      do i=1, 20
        coordinates      =[radii(i),Pi/2.0d0,0.0d0]
        massDistribution_=> node%massDistribution(                                 &
               &                                 componentType=componentTypeDisk, &
               &                                 massType     =massTypeStellar    &
                &                                 )
        surfaceDensityProfileExtract(i,1) =  massDistribution_%surfaceDensity(coordinates)
        surfaceDensityProfileExtract(i,2) =  radii(i)
      end do
    end select
    !![
        <objectDestructor name="massDistribution_"/>
    !!]
    return
  end function surfaceDensityProfileExtract

  subroutine surfaceDensityProfileNames(self,names,time)
    !!{
    Return the name of the surface density and radii at which it is extracted for disk component.
    !!}
    implicit none
    class           (nodePropertyExtractorSurfaceDensityDisk), intent(inout)                             :: self
    double precision                                         , intent(in   ), optional                   :: time
    type            (varying_string                         ), intent(inout), dimension(:) , allocatable :: names
    !$GLC attributes unused :: self, time
    
    allocate(names(2))
    names(1)=var_str(      'surfaceDensityProfileDisk')
    names(2)=var_str('surfaceDensityProfileRadiusDisk')
    return
  end subroutine surfaceDensityProfileNames

  subroutine surfaceDensityProfileDescriptions(self,descriptions,time)
    !!{
    Return descriptions of the {\normalfont \ttfamily surfaceDensityProfile} property.
    !!}
    implicit none
    class           (nodePropertyExtractorSurfaceDensityDisk), intent(inout)                             :: self
    double precision                                         , intent(in   ), optional                   :: time
    type            (varying_string                         ), intent(inout), dimension(:) , allocatable :: descriptions
    !$GLC attributes unused :: self, time

    allocate(descriptions(2))
    descriptions(1)="Surface density at a given radius [M☉/Mpc⁻²]."
    descriptions(2)="Radius at which surface density is output [Mpc]."
    return
  end subroutine surfaceDensityProfileDescriptions

  subroutine surfaceDensityProfileColumnDescriptions(self,descriptions,values,valuesDescription,valuesUnitsInSI,time)
    !!{
    Return column descriptions of the {\normalfont \ttfamily massProfile} property.
    !!}
    implicit none
    class           (nodePropertyExtractorSurfaceDensityDisk), intent(inout)                            :: self
    double precision                                         , intent(in   ), optional                  :: time
    type            (varying_string                         ), intent(inout), dimension(:), allocatable :: descriptions
    double precision                                         , intent(inout), dimension(:), allocatable :: values 
    type            (varying_string                         ), intent(  out)                            :: valuesDescription
    double precision                                         , intent(  out)                            :: valuesUnitsInSI
    !$GLC attributes unused :: time

    allocate(descriptions(20))
    allocate(values      (0))
    valuesDescription=var_str('')
    valuesUnitsInSI  =0.0d0
    descriptions     ='diskRadii'
    return
  end subroutine surfaceDensityProfileColumnDescriptions

  function surfaceDensityProfileUnitsInSI(self,time)
    !!{
    Return the units of the {\normalfont \ttfamily surfaceDensityProfile} properties in the SI system.
    !!}
    use :: Numerical_Constants_Astronomical, only : massSolar, megaParsec
    implicit none
    double precision                                         , allocatable  , dimension(:) :: surfaceDensityProfileUnitsInSI
    class           (nodePropertyExtractorSurfaceDensityDisk), intent(inout)               :: self
    double precision                                         , intent(in   ), optional               :: time
    !$GLC attributes unused :: self, time

    allocate(surfaceDensityProfileUnitsInSI(2))
    surfaceDensityProfileUnitsInSI(1)=massSolar/megaParsec**2
    surfaceDensityProfileUnitsInSI(2)=          megaParsec
    return
  end function surfaceDensityProfileUnitsInSI

