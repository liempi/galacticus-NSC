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

  !+    Contributions to this file made by: Matías Liempi

  !!{
  Implements a property extractor class the properties of nuclear star cluster when a black hole seed is formed using the model of \cite{vergara_global_2023}.
  !!}
  
  !![
  <nodePropertyExtractor name="nodePropertyExtractorBlackHoleSeedingDarkCores">
   <description>
    A property extractor class for the properties of the nuclear star cluster at the moment of the black hole formation.
   </description>
  </nodePropertyExtractor>
  !!]
  type, extends(nodePropertyExtractorTuple) :: nodePropertyExtractorBlackHoleSeedingDarkCores
     !!{
     A property extractor class for the velocity dispersion at a set of radii.
     !!}
     private
     integer  :: radiusNuclearStarClustersID          , blackHoleSeedMassID              , &
         &       velocityNuclearStarClustersID        , ageNuclearStarClustersID         , &
         &       gasMassNuclearStarClustersID         , criticalMassNuclearStarClustersID, &
         &       redshiftBlackHoleSeedFormationID     , stellarMassNuclearStarClustersID , &                                                                    
         &       mergerTreeWeightNuclearStarClustersID   
   contains
     procedure :: elementCount       => blackHoleSeedingDarkCoresElementCount
     procedure :: extract            => blackHoleSeedingDarkCoresExtract
     procedure :: names              => blackHoleSeedingDarkCoresNames
     procedure :: descriptions       => blackHoleSeedingDarkCoresDescriptions
     procedure :: unitsInSI          => blackHoleSeedingDarkCoresUnitsInSI
  end type nodePropertyExtractorBlackHoleSeedingDarkCores

  interface nodePropertyExtractorBlackHoleSeedingDarkCores
     !!{
     Constructors for the ``BlackHoleSeedingVergara2023'' output analysis class.
     !!}
     module procedure blackHoleSeedingDarkCoresConstructorParameters
     module procedure blackHoleSeedingDarkCoresConstructorInternal
  end interface nodePropertyExtractorBlackHoleSeedingDarkCores

contains

  function blackHoleSeedingDarkCoresConstructorParameters(parameters) result(self)
    !!{
    Constructor for the {\normalfont \ttfamily BlackHoleSeedingVergara2023} property extractor class.
    !!}
    use :: Input_Parameters, only : inputParameters
    implicit none
    type (nodePropertyExtractorBlackHoleSeedingDarkCores)                :: self
    type (inputParameters                                 ), intent(inout) :: parameters

    self=nodePropertyExtractorBlackHoleSeedingDarkCores()
    !![
    <inputParametersValidate source="parameters"/>
    !!]
    return
  end function blackHoleSeedingDarkCoresConstructorParameters

  function blackHoleSeedingDarkCoresConstructorInternal() result(self)
    !!{
    Internal constructor for the {\normalfont \ttfamily BlackHoleSeedingVergara2023} property extractor class.
    !!}
    implicit none
    type          (nodePropertyExtractorBlackHoleSeedingDarkCores) :: self
    !![
    <addMetaProperty   component="NSC"  name="blackHoleSeedMassFormed"             id="self%blackHoleSeedMassID"                   isEvolvable="no"  isCreator="no"/>
    <addMetaProperty   component="NSC"  name="redshiftBlackHoleSeedFormation"      id="self%redshiftBlackHoleSeedFormationID"      isEvolvable="no"  isCreator="no"/>
    !!]
    return
  end function blackHoleSeedingDarkCoresConstructorInternal

  integer function blackHoleSeedingDarkCoresElementCount(self,time)
    !!{
    Return the number of elements in the {\normalfont \ttfamily BlackHoleSeedingVergara2023} property extractors.
    !!}
    implicit none
    class           (nodePropertyExtractorBlackHoleSeedingDarkCores), intent(inout) :: self
    double precision                                                , intent(in   ) :: time
    !$GLC attributes unused :: time

    blackHoleSeedingDarkCoresElementCount=2
    return
  end function blackHoleSeedingDarkCoresElementCount

  function blackHoleSeedingDarkCoresExtract(self,node,time,instance)
    !!{
    Implement a {\normalfont \ttfamily BlackHoleSeedingVergara2023} property extractor.
    !!}
    use :: Galacticus_Nodes, only : nodeComponentNSC
    implicit none
    double precision                                                , dimension(:) , allocatable :: blackHoleSeedingDarkCoresExtract
    class           (nodePropertyExtractorBlackHoleSeedingDarkCores), intent(inout), target      :: self
    type            (treeNode                                      ), intent(inout), target      :: node
    double precision                                                , intent(in   )              :: time
    type            (multiCounter                                  ), intent(inout), optional    :: instance
    class           (nodeComponentNSC                              )               , pointer     :: nuclearStarCluster

    !$GLC attributes unused :: time, instance

    allocate(blackHoleSeedingDarkCoresExtract(2))
    nuclearStarCluster => node%NSC()
    select type (nuclearStarCluster)
    type is (nodeComponentNSC)
      ! Nuclear star cluster does not yet exist.
      blackHoleSeedingDarkCoresExtract=[0.0d0, 0.0d0]
    class default
      blackHoleSeedingDarkCoresExtract=[                                                                                        &
       &                                  nuclearStarCluster%floatRank0MetaPropertyGet(self%redshiftBlackHoleSeedFormationID ), &
       &                                  nuclearStarCluster%floatRank0MetaPropertyGet(self%blackHoleSeedMassID              )  &
       &                                ]
    end select
    return
  end function blackHoleSeedingDarkCoresExtract

  subroutine blackHoleSeedingDarkCoresNames(self,time,names)
    !!{
    Return the names of the {\normalfont \ttfamily BlackHoleSeedingVergara2023} properties.
    !!}
    implicit none
    class           (nodePropertyExtractorBlackHoleSeedingDarkCores), intent(inout)                             :: self
    double precision                                                , intent(in   )                             :: time
    type            (varying_string                                ), intent(inout), dimension(:) , allocatable :: names
    !$GLC attributes unused :: self, time
    
    allocate(names(2))
    names(1)=var_str('blackHoleFormationRedshift')
    names(2)=var_str('blackHoleSeedMass'         )
    return
  end subroutine blackHoleSeedingDarkCoresNames

  subroutine blackHoleSeedingDarkCoresDescriptions(self,time,descriptions)
    !!{
    Return descriptions of the {\normalfont \ttfamily blackHoleSeedingVergara2023} property.
    !!}
    implicit none
    class           (nodePropertyExtractorBlackHoleSeedingDarkCores), intent(inout)                             :: self
    double precision                                                , intent(in   )                             :: time
    type            (varying_string                                ), intent(inout), dimension(:) , allocatable :: descriptions
    !$GLC attributes unused :: time

    allocate(descriptions(2))
    descriptions(1)=var_str('Redshift at the formation of the black hole seed in dark cores.')
    descriptions(2)=var_str('Black hole seed mass in dark cores[M⊙].'                       )
    return
  end subroutine blackHoleSeedingDarkCoresDescriptions

  function blackHoleSeedingDarkCoresUnitsInSI(self,time)
    !!{
    Return the units of the {\normalfont \ttfamily BlackHoleSeedingVergara2023} properties in the SI system.
    !!}
    use :: Numerical_Constants_Astronomical, only : massSolar, megaParsec, gigayear
    use :: Numerical_Constants_Prefixes    , only : kilo
    implicit none
    double precision                                                , allocatable  , dimension(:) :: blackHoleSeedingDarkCoresUnitsInSI
    class           (nodePropertyExtractorBlackHoleSeedingDarkCores), intent(inout)               :: self
    double precision                                                , intent(in   )               :: time
    !$GLC attributes unused :: time

    allocate(blackHoleSeedingDarkCoresUnitsInSI(2))
    blackHoleSeedingDarkCoresUnitsInSI = [1.0d0, massSolar]
    return
  end function blackHoleSeedingDarkCoresUnitsInSI
