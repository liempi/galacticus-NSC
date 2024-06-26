!! Copyright 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018,
!!           2019, 2020, 2021, 2022, 2023, 2024
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
Contains a module which implements a stellar mass property extractor class.
!!}

  use :: Galactic_Structure, only : galacticStructureClass
  
  !![
  <nodePropertyExtractor name="nodePropertyExtractorMassStellar">
   <description>A stellar mass output analysis property extractor class.</description>
  </nodePropertyExtractor>
  !!]
  type, extends(nodePropertyExtractorScalar) :: nodePropertyExtractorMassStellar
     !!{
     A stellar mass property extractor class.
     !!}
     private
     class(galacticStructureClass), pointer :: galacticStructure_ => null()
   contains
     final     ::                massStellarDestructor
     procedure :: extract     => massStellarExtract
     procedure :: name        => massStellarName
     procedure :: description => massStellarDescription
     procedure :: unitsInSI   => massStellarUnitsInSI
     procedure :: quantity    => massStellarQuantity
  end type nodePropertyExtractorMassStellar

  interface nodePropertyExtractorMassStellar
     !!{
     Constructors for the ``massStellar'' property extractor class.
     !!}
     module procedure massStellarConstructorParameters
     module procedure massStellarConstructorInternal
  end interface nodePropertyExtractorMassStellar

contains

  function massStellarConstructorParameters(parameters) result(self)
    !!{
    Constructor for the ``massStellar'' output analysis property extractor class which takes a parameter set as input.
    !!}
    use :: Input_Parameters, only : inputParameters
    implicit none
    type (nodePropertyExtractorMassStellar)                :: self
    type (inputParameters                 ), intent(inout) :: parameters
    class(galacticStructureClass          ), pointer       :: galacticStructure_

    !![
    <objectBuilder class="galacticStructure" name="galacticStructure_" source="parameters"/>
    !!]
    self=nodePropertyExtractorMassStellar(galacticStructure_)
    !![
    <inputParametersValidate source="parameters"/>
    <objectDestructor name="galacticStructure_"/>
    !!]
    return
  end function massStellarConstructorParameters

  function massStellarConstructorInternal(galacticStructure_) result(self)
    !!{
    Internal constructor for the ``massStellar'' output analysis property extractor class.
    !!}
    implicit none
    type (nodePropertyExtractorMassStellar)                        :: self
    class(galacticStructureClass          ), intent(in   ), target :: galacticStructure_
    !![
    <constructorAssign variables="*galacticStructure_"/>
    !!]

    return
  end function massStellarConstructorInternal
  
  subroutine massStellarDestructor(self)
    !!{
    Destructor for the ``massStellar'' output analysis property extractor class.
    !!}
    implicit none
    type(nodePropertyExtractorMassStellar), intent(inout) :: self
    
    !![
    <objectDestructor name="self%galacticStructure_"/>
    !!]
    return
  end subroutine massStellarDestructor

  double precision function massStellarExtract(self,node,instance)
    !!{
    Implement a massStellar output analysis.
    !!}
    use :: Galactic_Structure_Options, only : massTypeStellar, radiusLarge
    implicit none
    class(nodePropertyExtractorMassStellar), intent(inout), target   :: self
    type (treeNode                        ), intent(inout), target   :: node
    type (multiCounter                    ), intent(inout), optional :: instance
    !$GLC attributes unused :: self, instance

    massStellarExtract=self%galacticStructure_%massEnclosed(node,radiusLarge,massType=massTypeStellar)
    return
  end function massStellarExtract

  function massStellarName(self)
    !!{
    Return the name of the massStellar property.
    !!}
    implicit none
    type (varying_string                  )                :: massStellarName
    class(nodePropertyExtractorMassStellar), intent(inout) :: self
    !$GLC attributes unused :: self

    massStellarName=var_str('massStellarTotal')
    return
  end function massStellarName

  function massStellarDescription(self)
    !!{
    Return a description of the massStellar property.
    !!}
    implicit none
    type (varying_string                  )                :: massStellarDescription
    class(nodePropertyExtractorMassStellar), intent(inout) :: self
    !$GLC attributes unused :: self

    massStellarDescription=var_str('The total mass of stars in this node [M☉].')
    return
  end function massStellarDescription

  double precision function massStellarUnitsInSI(self)
    !!{
    Return the units of the massStellar property in the SI system.
    !!}
    use :: Numerical_Constants_Astronomical, only : massSolar
    implicit none
    class(nodePropertyExtractorMassStellar), intent(inout) :: self
    !$GLC attributes unused :: self

    massStellarUnitsInSI=massSolar
    return
  end function massStellarUnitsInSI


  function massStellarQuantity(self)
    !!{
    Return the class of the stellar luminosity property.
    !!}
    use :: Output_Analyses_Options, only : outputAnalysisPropertyQuantityMass
    implicit none
    type (enumerationOutputAnalysisPropertyQuantityType)                :: massStellarQuantity
    class(nodePropertyExtractorMassStellar             ), intent(inout) :: self
    !$GLC attributes unused :: self

    massStellarQuantity=outputAnalysisPropertyQuantityMass
    return
  end function massStellarQuantity
