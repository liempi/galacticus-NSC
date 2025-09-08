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
  <nodePropertyExtractor name="nodePropertyExtractorBlackHoleFormationMass">
   <description>
    A node property extractor class which extracts the formation redshift for black hole seeds.
   </description>
  </nodePropertyExtractor>
  !!]
  type, extends(nodePropertyExtractorScalar) :: nodePropertyExtractorBlackHoleFormationMass
     !!{
     A node property extractor class which extracts the mass of the black hole seeds.
     !!}
     private
     integer :: blackHoleSeedsFormationMassID
   contains
     procedure :: extract     => blackHoleFormationMassExtract
     procedure :: name        => blackHoleFormationMassName
     procedure :: description => blackHoleFormationMassDescription
     procedure :: unitsInSI   => blackHoleFormationMassUnitsInSI
  end type nodePropertyExtractorBlackHoleFormationMass

  interface nodePropertyExtractorBlackHoleFormationMass
     !!{
     Constructors for the {\normalfont \ttfamily blackHoleFormationMass} node property extractor class.
     !!}
     module procedure blackHoleFormationMassConstructorParameters
     module procedure blackHoleFormationMassConstructorInternal
  end interface nodePropertyExtractorBlackHoleFormationMass

contains

  function blackHoleFormationMassConstructorParameters(parameters) result(self)
    !!{
    Constructor for the {\normalfont \ttfamily blackHoleFormationMass} node property extractor class which takes a parameter set as input.
    !!}
    use :: Input_Parameters, only : inputParameters
    implicit none
    type(nodePropertyExtractorBlackHoleFormationMass)                :: self
    type(inputParameters                                ), intent(inout) :: parameters

    self=nodePropertyExtractorBlackHoleFormationMass()
    !![
    <inputParametersValidate source="parameters"/>
    !!]
    return
  end function blackHoleFormationMassConstructorParameters

  function blackHoleFormationMassConstructorInternal() result(self)
    !!{
    Internal constructor for the {\normalfont \ttfamily blackHoleFormationMass} node property extractor class.
    !!}
    implicit none
    type(nodePropertyExtractorBlackHoleFormationMass) :: self
    
    !![
    <addMetaProperty component="blackHole" name="blackHoleSeedsFormationMass" type="float" id="self%blackHoleSeedsFormationMassID" isCreator="no"/>
    !!]
    return
  end function blackHoleFormationMassConstructorInternal

  function blackHoleFormationMassExtract(self,node,instance)
    !!{
    Implement a {\normalfont \ttfamily blackHoleFormationMass} node property extractor.
    !!}
    use :: Galacticus_Nodes, only : nodeComponentBlackHole, nodeComponentBlackHoleStandard
    implicit none
    double precision                                                                       :: blackHoleFormationMassExtract
    class           (nodePropertyExtractorBlackHoleFormationMass), intent(inout), target   :: self
    type            (treeNode                                   ), intent(inout), target   :: node
    type            (multiCounter                               ), intent(inout), optional :: instance
    class           (nodeComponentBlackHole                     )               , pointer  :: blackHole
    !$GLC attributes unused :: instance

    blackHole => node%blackHole()
    select type (blackHole)
    class is (nodeComponentBlackHole)
       ! No black hole exists - the formation redshift is undetermined.
       blackHoleFormationMassExtract=-1.0d0
    class default
       ! Extract the formation redshift.
       blackHoleFormationMassExtract=blackHole%floatRank0MetaPropertyGet(self%blackHoleSeedsFormationMassID)
    end select
    return
  end function blackHoleFormationMassExtract

  function blackHoleFormationMassName(self)
    !!{
    Return the name of the blackHoleFormationMass property.
    !!}
    implicit none
    type (varying_string                             )                :: blackHoleFormationMassName
    class(nodePropertyExtractorBlackHoleFormationMass), intent(inout) :: self
    !$GLC attributes unused :: self
  
    blackHoleFormationMassName=var_str('blackHoleSeedsMass')
    return
  end function blackHoleFormationMassName
  
  function blackHoleFormationMassDescription(self) result(description)
    !!{
    Return a description of the blackHoleFormationMass property.
    !!}
    implicit none
    type     (varying_string                                 )                :: description
    class    (nodePropertyExtractorBlackHoleFormationMass), intent(inout) :: self
    !$GLC attributes unused :: self

    description='Indicates the mass of the black hole seed (M☉).'
    return
  end function blackHoleFormationMassDescription

  double precision function blackHoleFormationMassUnitsInSI(self)
    !!{
    Return the units of the bound mass radius property in the SI system.
    !!}
    implicit none
    class(nodePropertyExtractorBlackHoleFormationMass), intent(inout) :: self
    !$GLC attributes unused :: self

    blackHoleFormationMassUnitsInSI=1.0d0
    return
  end function blackHoleFormationMassUnitsInSI

