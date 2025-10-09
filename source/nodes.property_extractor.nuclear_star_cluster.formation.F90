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
  Implements a property extractor class the properties of nuclear star cluster when it reaches a minimum mass.
  !!}
  
  !![
  <nodePropertyExtractor name="nodePropertyExtractorNuclearStarClusterFormation">
   <description>
    A property extractor class for the properties of the nuclear star cluster when it reaches a minimum mass.
   </description>
  </nodePropertyExtractor>
  !!]
  type, extends(nodePropertyExtractorTuple) :: nodePropertyExtractorNuclearStarClusterFormation
     !!{
     A property extractor class for the velocity dispersion at a set of radii.
     !!}
     private
    integer :: nuclearStarClusterFormationTimeID       , nuclearStarClusterFormationStarFormationRateID, &
         &     nuclearStarClusterFormationStellarMassID
   contains
     procedure :: elementCount       => nuclearStarClusterFormationElementCount
     procedure :: extract            => nuclearStarClusterFormationExtract
     procedure :: names              => nuclearStarClusterFormationNames
     procedure :: descriptions       => nuclearStarClusterFormationDescriptions
     procedure :: unitsInSI          => nuclearStarClusterFormationUnitsInSI
  end type nodePropertyExtractorNuclearStarClusterFormation

  interface nodePropertyExtractorNuclearStarClusterFormation
     !!{
     Constructors for the ``nuclearStarClusterFormation'' output analysis class.
     !!}
     module procedure nuclearStarClusterFormationConstructorParameters
     module procedure nuclearStarClusterFormationConstructorInternal
  end interface nodePropertyExtractorNuclearStarClusterFormation

contains

  function nuclearStarClusterFormationConstructorParameters(parameters) result(self)
    !!{
    Constructor for the {\normalfont \ttfamily nuclearStarClusterFormation} property extractor class.
    !!}
    use :: Input_Parameters, only : inputParameters
    implicit none
    type (nodePropertyExtractorNuclearStarClusterFormation)                :: self
    type (inputParameters                                  ), intent(inout) :: parameters

    self=nodePropertyExtractorNuclearStarClusterFormation()
    !![
    <inputParametersValidate source="parameters"/>
    !!]
    return
  end function nuclearStarClusterFormationConstructorParameters

  function nuclearStarClusterFormationConstructorInternal() result(self)
    !!{
    Internal constructor for the {\normalfont \ttfamily nuclearStarClusterFormation} property extractor class.
    !!}
    implicit none
    type (nodePropertyExtractorNuclearStarClusterFormation) :: self
    !![
    <addMetaProperty component="NSC" name="nuclearStarClusterFormationTime"              id="self%nuclearStarClusterFormationTimeID"              isEvolvable="no" isCreator="no"/>
    <addMetaProperty component="NSC" name="nuclearStarClusterFormationStarFormationRate" id="self%nuclearStarClusterFormationStarFormationRateID" isEvolvable="no" isCreator="no"/>
    <addMetaProperty component="NSC" name="nuclearStarClusterFormationStellarMass"       id="self%nuclearStarClusterFormationStellarMassID"       isEvolvable="no" isCreator="no"/>
    !!]
    return
  end function nuclearStarClusterFormationConstructorInternal

  integer function nuclearStarClusterFormationElementCount(self,time)
    !!{
    Return the number of elements in the {\normalfont \ttfamily nuclearStarClusterFormation} property extractors.
    !!}
    implicit none
    class           (nodePropertyExtractorNuclearStarClusterFormation), intent(inout) :: self
    double precision                                                  , intent(in   ) :: time
    !$GLC attributes unused :: time

    nuclearStarClusterFormationElementCount=3
    return
  end function nuclearStarClusterFormationElementCount

  function nuclearStarClusterFormationExtract(self,node,time,instance)
    !!{
    Implement a {\normalfont \ttfamily nuclearStarClusterFormation} property extractor.
    !!}
    use :: Galacticus_Nodes, only : nodeComponentNSC
    implicit none
    double precision                                                  , dimension(:) , allocatable :: nuclearStarClusterFormationExtract
    class           (nodePropertyExtractorNuclearStarClusterFormation), intent(inout), target      :: self
    type            (treeNode                                        ), intent(inout), target      :: node
    double precision                                                  , intent(in   )              :: time
    type            (multiCounter                                    ), intent(inout), optional    :: instance
    class           (nodeComponentNSC                                )               , pointer     :: nuclearStarCluster

    !$GLC attributes unused :: time, instance

    allocate(nuclearStarClusterFormationExtract(3))
    nuclearStarCluster => node%NSC()
    select type (nuclearStarCluster)
    type is (nodeComponentNSC)
      ! Nuclear star cluster does not yet exist.
      nuclearStarClusterFormationExtract=[       & 
        &                                 0.0d0, &
        &                                 0.0d0, &
        &                                 0.0d0  &
        &                                ]
    class default
      nuclearStarClusterFormationExtract=[                                                                                                   &
       &                                  nuclearStarCluster%floatRank0MetaPropertyGet(self%nuclearStarClusterFormationTimeID             ), &
       &                                  nuclearStarCluster%floatRank0MetaPropertyGet(self%nuclearStarClusterFormationStarFormationRateID), &
       &                                  nuclearStarCluster%floatRank0MetaPropertyGet(self%nuclearStarClusterFormationStellarMassID      )  &
       &                                 ]
    end select
    return
  end function nuclearStarClusterFormationExtract

  subroutine nuclearStarClusterFormationNames(self,time,names)
    !!{
    Return the names of the {\normalfont \ttfamily nuclearStarClusterFormation} properties.
    !!}
    implicit none
    class           (nodePropertyExtractorNuclearStarClusterFormation), intent(inout)                             :: self
    double precision                                                  , intent(in   )                             :: time
    type            (varying_string                                  ), intent(inout), dimension(:) , allocatable :: names
    !$GLC attributes unused :: self, time
    
    allocate(names(3))
    names(1)=var_str('nuclearStarClusterFormationTime'             )
    names(2)=var_str('nuclearStarClusterFormationStarFormationRate')
    names(3)=var_str('nuclearStarClusterFormationStellarMass'      )
    return
  end subroutine nuclearStarClusterFormationNames

  subroutine nuclearStarClusterFormationDescriptions(self,time,descriptions)
    !!{
    Return descriptions of the {\normalfont \ttfamily nuclearStarClusterFormation} property.
    !!}
    implicit none
    class           (nodePropertyExtractorNuclearStarClusterFormation), intent(inout)                             :: self
    double precision                                                  , intent(in   )                             :: time
    type            (varying_string                                  ), intent(inout), dimension(:) , allocatable :: descriptions
    !$GLC attributes unused :: time

    allocate(descriptions(3))
    descriptions(1)=var_str('Hubble time at the moment that the nuclear star cluster reach a minium mass [Gyr].')
    descriptions(2)=var_str('Star formation rate of the nuclear star cluster [M⊙/Gyr].'                         )
    descriptions(3)=var_str('Stellar mass of the nuclear star cluster [M⊙].'                                    )
    return
  end subroutine nuclearStarClusterFormationDescriptions

  function nuclearStarClusterFormationUnitsInSI(self,time)
    !!{
    Return the units of the {\normalfont \ttfamily nuclearStarClusterFormation} properties in the SI system.
    !!}
    use :: Numerical_Constants_Astronomical, only : massSolar, megaParsec, gigayear
    use :: Numerical_Constants_Prefixes    , only : kilo
    implicit none
    double precision                                                  , allocatable  , dimension(:) :: nuclearStarClusterFormationUnitsInSI
    class           (nodePropertyExtractorNuclearStarClusterFormation), intent(inout)               :: self
    double precision                                                  , intent(in   )               :: time
    !$GLC attributes unused :: time

    allocate(nuclearStarClusterFormationUnitsInSI(3))
    nuclearStarClusterFormationUnitsInSI=[                    &
     &                                    gigayear          , &
     &                                    massSolar/gigayear, &
     &                                    massSolar           &
     &                                   ]
    return
  end function nuclearStarClusterFormationUnitsInSI
