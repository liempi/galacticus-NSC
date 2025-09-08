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
  <nodePropertyExtractor name="nodePropertyExtractorBlackHoleFormationTimescale">
   <description>
    A node property extractor class which extracts the formation redshift for black hole seeds.
   </description>
  </nodePropertyExtractor>
  !!]
  type, extends(nodePropertyExtractorScalar) :: nodePropertyExtractorBlackHoleFormationTimescale
     !!{
     A node property extractor class which extracts the formation timescale for black hole seeds.
     !!}
     private
     integer :: blackHoleSeedsFormationTimescaleID
   contains
     procedure :: extract     => blackHoleFormationTimescaleExtract
     procedure :: name        => blackHoleFormationTimescaleName
     procedure :: description => blackHoleFormationTimescaleDescription
     procedure :: unitsInSI   => blackHoleFormationTimescaleUnitsInSI
  end type nodePropertyExtractorBlackHoleFormationTimescale

  interface nodePropertyExtractorBlackHoleFormationTimescale
     !!{
     Constructors for the {\normalfont \ttfamily blackHoleFormationTimescale} node property extractor class.
     !!}
     module procedure blackHoleFormationTimescaleConstructorParameters
     module procedure blackHoleFormationTimescaleConstructorInternal
  end interface nodePropertyExtractorBlackHoleFormationTimescale

contains

  function blackHoleFormationTimescaleConstructorParameters(parameters) result(self)
    !!{
    Constructor for the {\normalfont \ttfamily blackHoleFormationTimescale} node property extractor class which takes a parameter set as input.
    !!}
    use :: Input_Parameters, only : inputParameters
    implicit none
    type(nodePropertyExtractorBlackHoleFormationTimescale)                :: self
    type(inputParameters                                ), intent(inout) :: parameters

    self=nodePropertyExtractorBlackHoleFormationTimescale()
    !![
    <inputParametersValidate source="parameters"/>
    !!]
    return
  end function blackHoleFormationTimescaleConstructorParameters

  function blackHoleFormationTimescaleConstructorInternal() result(self)
    !!{
    Internal constructor for the {\normalfont \ttfamily blackHoleFormationTimescale} node property extractor class.
    !!}
    implicit none
    type(nodePropertyExtractorBlackHoleFormationTimescale) :: self
    
    !![
    <addMetaProperty component="blackHole" name="blackHoleSeedsFormationTimescale" type="float" id="self%blackHoleSeedsFormationTimescaleID" isCreator="no"/>
    !!]
    return
  end function blackHoleFormationTimescaleConstructorInternal

  function blackHoleFormationTimescaleExtract(self,node,instance)
    !!{
    Implement a {\normalfont \ttfamily blackHoleFormationTimescale} node property extractor.
    !!}
    use :: Galacticus_Nodes, only : nodeComponentBlackHole, nodeComponentBlackHoleStandard
    implicit none
    double precision                                                                           :: blackHoleFormationTimescaleExtract
    class           (nodePropertyExtractorBlackHoleFormationTimescale), intent(inout), target   :: self
    type            (treeNode                                        ), intent(inout), target   :: node
    type            (multiCounter                                    ), intent(inout), optional :: instance
    class           (nodeComponentBlackHole                          )               , pointer  :: blackHole
    !$GLC attributes unused :: instance

    blackHole => node%blackHole()
    select type (blackHole)
    class is (nodeComponentBlackHole)
       ! No black hole exists - the formation redshift is undetermined.
       blackHoleFormationTimescaleExtract=-1.0d0
    class default
       ! Extract the formation redshift.
       blackHoleFormationTimescaleExtract=blackHole%floatRank0MetaPropertyGet(self%blackHoleSeedsFormationTimescaleID)
    end select
    return
  end function blackHoleFormationTimescaleExtract

  function blackHoleFormationTimescaleName(self)
    !!{
    Return the name of the blackHoleFormationTimescale property.
    !!}
    implicit none
    type (varying_string                                  )                :: blackHoleFormationTimescaleName
    class(nodePropertyExtractorBlackHoleFormationTimescale), intent(inout) :: self
    !$GLC attributes unused :: self
  
    blackHoleFormationTimescaleName=var_str('blackHoleSeedsFormationTimescale')
    return
  end function blackHoleFormationTimescaleName
  
  function blackHoleFormationTimescaleDescription(self) result(description)
    !!{
    Return a description of the blackHoleFormationTimescale property.
    !!}
    implicit none
    type     (varying_string                                  )                :: description
    class    (nodePropertyExtractorBlackHoleFormationTimescale), intent(inout) :: self
    !$GLC attributes unused :: self

    description='Indicates the timescale (in Gyr) of the black hole seed formation.'
    return
  end function blackHoleFormationTimescaleDescription

  double precision function blackHoleFormationTimescaleUnitsInSI(self)
    !!{
    Return the units of the bound mass radius property in the SI system.
    !!}
    use :: Numerical_Constants_Astronomical, only : gigaYear
    implicit none
    class(nodePropertyExtractorBlackHoleFormationTimescale), intent(inout) :: self
    !$GLC attributes unused :: self

    blackHoleFormationTimescaleUnitsInSI=gigaYear
    return
  end function blackHoleFormationTimescaleUnitsInSI

