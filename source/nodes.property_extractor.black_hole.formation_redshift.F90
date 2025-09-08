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
  <nodePropertyExtractor name="nodePropertyExtractorBlackHoleFormationRedshift">
   <description>
    A node property extractor class which extracts the formation redshift for black hole seeds.
   </description>
  </nodePropertyExtractor>
  !!]
  type, extends(nodePropertyExtractorScalar) :: nodePropertyExtractorBlackHoleFormationRedshift
     !!{
     A node property extractor class which extracts the formation redshift for black hole seeds.
     !!}
     private
     integer :: blackHoleSeedsFormationRedshiftID
   contains
     procedure :: extract     => blackHoleFormationRedshiftExtract
     procedure :: name        => blackHoleFormationRedshiftName
     procedure :: description => blackHoleFormationRedshiftDescription
     procedure :: unitsInSI   => blackHoleFormationRedshiftUnitsInSI
  end type nodePropertyExtractorBlackHoleFormationRedshift

  interface nodePropertyExtractorBlackHoleFormationRedshift
     !!{
     Constructors for the {\normalfont \ttfamily blackHoleFormationRedshift} node property extractor class.
     !!}
     module procedure blackHoleFormationRedshiftConstructorParameters
     module procedure blackHoleFormationRedshiftConstructorInternal
  end interface nodePropertyExtractorBlackHoleFormationRedshift

contains

  function blackHoleFormationRedshiftConstructorParameters(parameters) result(self)
    !!{
    Constructor for the {\normalfont \ttfamily blackHoleFormationRedshift} node property extractor class which takes a parameter set as input.
    !!}
    use :: Input_Parameters, only : inputParameters
    implicit none
    type(nodePropertyExtractorBlackHoleFormationRedshift)                :: self
    type(inputParameters                                ), intent(inout) :: parameters

    self=nodePropertyExtractorBlackHoleFormationRedshift()
    !![
    <inputParametersValidate source="parameters"/>
    !!]
    return
  end function blackHoleFormationRedshiftConstructorParameters

  function blackHoleFormationRedshiftConstructorInternal() result(self)
    !!{
    Internal constructor for the {\normalfont \ttfamily blackHoleFormationRedshift} node property extractor class.
    !!}
    implicit none
    type(nodePropertyExtractorBlackHoleFormationRedshift) :: self
    
    !![
    <addMetaProperty component="blackHole" name="blackHoleSeedsFormationRedshift" type="float" id="self%blackHoleSeedsFormationRedshiftID" isCreator="no"/>
    !!]
    return
  end function blackHoleFormationRedshiftConstructorInternal

  function blackHoleFormationRedshiftExtract(self,node,instance)
    !!{
    Implement a {\normalfont \ttfamily blackHoleFormationRedshift} node property extractor.
    !!}
    use :: Galacticus_Nodes, only : nodeComponentBlackHole, nodeComponentBlackHoleStandard
    implicit none
    double precision                                                                           :: blackHoleFormationRedshiftExtract
    class           (nodePropertyExtractorBlackHoleFormationRedshift), intent(inout), target   :: self
    type            (treeNode                                       ), intent(inout), target   :: node
    type            (multiCounter                                   ), intent(inout), optional :: instance
    class           (nodeComponentBlackHole                         )               , pointer  :: blackHole
    !$GLC attributes unused :: instance

    blackHole => node%blackHole()
    select type (blackHole)
    class is (nodeComponentBlackHole)
       ! No black hole exists - the formation redshift is undetermined.
       blackHoleFormationRedshiftExtract=-1.0d0
    class default
       ! Extract the formation redshift.
       blackHoleFormationRedshiftExtract=blackHole%floatRank0MetaPropertyGet(self%blackHoleSeedsFormationRedshiftID)
    end select
    return
  end function blackHoleFormationRedshiftExtract

  function blackHoleFormationRedshiftName(self)
    !!{
    Return the name of the blackHoleFormationRedshift property.
    !!}
    implicit none
    type (varying_string                                 )                :: blackHoleFormationRedshiftName
    class(nodePropertyExtractorBlackHoleFormationRedshift), intent(inout) :: self
    !$GLC attributes unused :: self
  
    blackHoleFormationRedshiftName=var_str('blackHoleSeedsFormationRedshift')
    return
  end function blackHoleFormationRedshiftName
  
  function blackHoleFormationRedshiftDescription(self) result(description)
    !!{
    Return a description of the blackHoleFormationRedshift property.
    !!}
    implicit none
    type     (varying_string                                 )                :: description
    class    (nodePropertyExtractorBlackHoleFormationRedshift), intent(inout) :: self
    !$GLC attributes unused :: self

    description='Indicates the redshift of the black hole seed formation.'
    return
  end function blackHoleFormationRedshiftDescription

  double precision function blackHoleFormationRedshiftUnitsInSI(self)
    !!{
    Return the units of the bound mass radius property in the SI system.
    !!}
    implicit none
    class(nodePropertyExtractorBlackHoleFormationRedshift), intent(inout) :: self
    !$GLC attributes unused :: self

    blackHoleFormationRedshiftUnitsInSI=1.0d0
    return
  end function blackHoleFormationRedshiftUnitsInSI

