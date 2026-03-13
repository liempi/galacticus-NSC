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
  <nodePropertyExtractor name="nodePropertyExtractorBlackHoleSeedRedshift">
   <description>
    A node property extractor class which extracts the formation redshift for black hole seeds.
   </description>
  </nodePropertyExtractor>
  !!]
  type, extends(nodePropertyExtractorList) :: nodePropertyExtractorBlackHoleSeedRedshift
     !!{
     A node property extractor class which extracts the formation redshift for black hole seeds.
     !!}
     private
     integer :: blackHoleSeedRedshiftID
   contains
     procedure :: elementCount => blackHoleSeedRedshiftElementCount 
     procedure :: extract      => blackHoleSeedRedshiftExtract
     procedure :: names        => blackHoleSeedRedshiftNames
     procedure :: descriptions => blackHoleSeedRedshiftDescriptions
     procedure :: unitsInSI    => blackHoleSeedRedshiftUnitsInSI
  end type nodePropertyExtractorBlackHoleSeedRedshift

  interface nodePropertyExtractorBlackHoleSeedRedshift
     !!{
     Constructors for the {\normalfont \ttfamily BlackHoleSeedRedshift} node property extractor class.
     !!}
     module procedure blackHoleSeedRedshiftConstructorParameters
     module procedure blackHoleSeedRedshiftConstructorInternal
  end interface nodePropertyExtractorBlackHoleSeedRedshift

contains

  function blackHoleSeedRedshiftConstructorParameters(parameters) result(self)
    !!{
    Constructor for the {\normalfont \ttfamily BlackHoleSeedRedshift} node property extractor class which takes a parameter set as input.
    !!}
    use :: Input_Parameters, only : inputParameters
    implicit none
    type(nodePropertyExtractorBlackHoleSeedRedshift)                :: self
    type(inputParameters                           ), intent(inout) :: parameters

    self=nodePropertyExtractorBlackHoleSeedRedshift()
    !![
    <inputParametersValidate source="parameters"/>
    !!]
    return
  end function blackHoleSeedRedshiftConstructorParameters

  function blackHoleSeedRedshiftConstructorInternal() result(self)
    !!{
    Internal constructor for the {\normalfont \ttfamily BlackHoleSeedRedshift} node property extractor class.
    !!}
    implicit none
    type(nodePropertyExtractorBlackHoleSeedRedshift) :: self
    
    !![
    <addMetaProperty component="blackHole" name="blackHoleSeedRedshift" type="float" id="self%blackHoleSeedRedshiftID" isCreator="no"/>
    !!]
    return
  end function blackHoleSeedRedshiftConstructorInternal

  integer function blackHoleSeedRedshiftElementCount(self)
    !!{
    Return a count of the number of properties extracted.
    !!}
    implicit none
    class(nodePropertyExtractorBlackHoleSeedRedshift), intent(inout) :: self

    blackHoleSeedRedshiftElementCount=1
    return
  end function blackHoleSeedRedshiftElementCount

  function blackHoleSeedRedshiftExtract(self,node,instance) result(seedRedshift)
    !!{
    Implement a {\normalfont \ttfamily BlackHoleSeedRedshift} node property extractor.
    !!}
    use :: Galacticus_Nodes, only : nodeComponentBlackHole
    implicit none
    double precision                                            , dimension(:,:), allocatable  :: seedRedshift
    class           (nodePropertyExtractorBlackHoleSeedRedshift), intent(inout)                :: self
    type            (treeNode                                  ), intent(inout)                :: node
    type            (multiCounter                              ), intent(inout) , optional     :: instance
    class           (nodeComponentBlackHole                    )                , pointer      :: blackHole
    integer                                                                                    :: i           , countBlackHoles
    !$GLC attributes unused :: instance

    countBlackHoles=node%blackHoleCount()
    allocate(seedRedshift(countBlackHoles,1))
    do i=1,countBlackHoles
       blackHole         => node     %                blackHole(instance=i                  )
       seedRedshift(i,1) =  blackHole%floatRank0MetaPropertyGet(self%blackHoleSeedRedshiftID)
    end do
    return
  end function blackHoleSeedRedshiftExtract

  subroutine blackHoleSeedRedshiftNames(self,names)
    !!{
    Return the name of the BlackHoleSeedRedshift property.
    !!}
    implicit none
    class(nodePropertyExtractorBlackHoleSeedRedshift), intent(inout) :: self
    type (varying_string                            ), intent(inout), dimension(:) , allocatable :: names
    !$GLC attributes unused :: self
    
    allocate(names(1))
    names(1)=var_str('blackHoleSeedRedshift')
    return
  end subroutine blackHoleSeedRedshiftNames
  
  subroutine blackHoleSeedRedshiftDescriptions(self,descriptions)
    !!{
    Return a description of the BlackHoleSeedRedshift property.
    !!}
    implicit none
    class(nodePropertyExtractorBlackHoleSeedRedshift), intent(inout)                             :: self
    type (varying_string                            ), intent(inout), dimension(:) , allocatable :: descriptions

    !$GLC attributes unused :: self
    allocate(descriptions(1))
    descriptions(1)=var_str('Indicates the redshift of the black hole seed formation.')
    return
  end subroutine blackHoleSeedRedshiftDescriptions

  function blackHoleSeedRedshiftUnitsInSI(self) result(unitsInSI)
    !!{
    Return the units of the bound mass radius property in the SI system.
    !!}
    implicit none
    class           (nodePropertyExtractorBlackHoleSeedRedshift), intent(inout)              :: self
    double precision                                            , dimension(:) , allocatable :: unitsInSI
    !$GLC attributes unused :: self

    allocate(unitsInSI(1))
    unitsInSI(1)=1.0d0
    return
  end function blackHoleSeedRedshiftUnitsInSI

