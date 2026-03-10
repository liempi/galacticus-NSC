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
  Implements a node property extractor which reports the formation redshift for the central black hole of a given node.
  !!}

  !![
  <nodePropertyExtractor name="nodePropertyExtractorBlackHoleSeedMass">
   <description>
    A node property extractor class which extracts mass of the black hole seed.
   </description>
  </nodePropertyExtractor>
  !!]
  type, extends(nodePropertyExtractorList) :: nodePropertyExtractorBlackHoleSeedMass
     !!{
     A node property extractor class which extracts the mass of the black hole seeds.
     !!}
     private
     integer :: blackHoleSeedMassID
   contains
     procedure :: elementCount => blackHoleSeedElementCount
     procedure :: extract      => blackHoleSeedMassExtract
     procedure :: names        => blackHoleSeedMassNames
     procedure :: descriptions => blackHoleSeedMassDescriptions
     procedure :: unitsInSI    => blackHoleSeedMassUnitsInSI
  end type nodePropertyExtractorBlackHoleSeedMass

  interface nodePropertyExtractorBlackHoleSeedMass
     !!{
     Constructors for the {\normalfont \ttfamily blackHoleSeedMass} node property extractor class.
     !!}
     module procedure blackHoleSeedMassConstructorParameters
     module procedure blackHoleSeedMassConstructorInternal
  end interface nodePropertyExtractorBlackHoleSeedMass

contains

  function blackHoleSeedMassConstructorParameters(parameters) result(self)
    !!{
    Constructor for the {\normalfont \ttfamily blackHoleSeedMass} node property extractor class which takes a parameter set as input.
    !!}
    use :: Input_Parameters, only : inputParameters
    implicit none
    type(nodePropertyExtractorBlackHoleSeedMass)                :: self
    type(inputParameters                       ), intent(inout) :: parameters

    self=nodePropertyExtractorBlackHoleSeedMass()
    !![
    <inputParametersValidate source="parameters"/>
    !!]
    return
  end function blackHoleSeedMassConstructorParameters

  function blackHoleSeedMassConstructorInternal() result(self)
    !!{
    Internal constructor for the {\normalfont \ttfamily blackHoleSeedMass} node property extractor class.
    !!}
    implicit none
    type(nodePropertyExtractorBlackHoleSeedMass) :: self
    
    !![
    <addMetaProperty component="blackHole" name="blackHoleSeedMass" type="float" id="self%blackHoleSeedMassID" isCreator="no"/>
    !!]
    return
  end function blackHoleSeedMassConstructorInternal

  integer function blackHoleSeedElementCount(self)
    !!{
    Return a count of the number of properties extracted.
    !!}
    implicit none
    class(nodePropertyExtractorBlackHoleSeedMass), intent(inout) :: self

    blackHoleSeedElementCount=1
    return
  end function blackHoleSeedElementCount

  function blackHoleSeedMassExtract(self,node,instance) result(seedMass)
    !!{
    Implement a {\normalfont \ttfamily blackHoleSeedMass} node property extractor.
    !!}
    use :: Galacticus_Nodes, only : nodeComponentBlackHole, nodeComponentBlackHoleStandard
    implicit none
    double precision                                        , dimension(:,:), allocatable :: seedMass
    class           (nodePropertyExtractorBlackHoleSeedMass), intent(inout)               :: self
    type            (treeNode                              ), intent(inout)               :: node
    type            (multiCounter                          ), intent(inout) , optional    :: instance
    class           (nodeComponentBlackHole                )                , pointer     :: blackHole
    integer                                                                               :: i        , countBlackHoles
    !$GLC attributes unused :: instance
    
    countBlackHoles=node%blackHoleCount()
    allocate(seedMass(countBlackHoles,1))
    do i=1,countBlackHoles
       blackHole      => node     %                blackHole(instance=i              )
       seedMass (i,1) =  blackHole%floatRank0MetaPropertyGet(self%blackHoleSeedMassID)
    end do
    return
  end function blackHoleSeedMassExtract

  subroutine blackHoleSeedMassNames(self,names)
    !!{
    Return the name of the blackHoleSeedMass property.
    !!}
    implicit none
    class(nodePropertyExtractorBlackHoleSeedMass), intent(inout)                             :: self
    type (varying_string                        ), intent(inout), dimension(:) , allocatable :: names
    !$GLC attributes unused :: self
    
    allocate(names(1))
    names(1)=var_str('blackHoleSeedMass')
    return
  end subroutine blackHoleSeedMassNames
  
  subroutine blackHoleSeedMassDescriptions(self,descriptions) 
    !!{
    Return a description of the blackHoleSeedMass property.
    !!}
    implicit none
    class(nodePropertyExtractorBlackHoleSeedMass), intent(inout) :: self
    type (varying_string                        ), intent(inout), dimension(:) , allocatable :: descriptions
    !$GLC attributes unused :: self
    
    allocate(descriptions(1))
    descriptions(1)=var_str('Indicates the mass of the black hole seed (M☉).')
    return
  end subroutine blackHoleSeedMassDescriptions

  function blackHoleSeedMassUnitsInSI(self) result(unitsInSI)
    !!{
    Return the units of the bound mass radius property in the SI system.
    !!}
    use :: Numerical_Constants_Astronomical, only : massSolar
    implicit none
    class           (nodePropertyExtractorBlackHoleSeedMass), intent(inout)              :: self
    double precision                                        , dimension(:) , allocatable :: unitsInSI

    !$GLC attributes unused :: self

    allocate(unitsInSI(1))
    unitsInSI(1)=massSolar
    return
  end function blackHoleSeedMassUnitsInSI

