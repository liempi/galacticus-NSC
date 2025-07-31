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
  type, extends(nodePropertyExtractorScalar) :: nodePropertyExtractorBlackHoleSeedingDarkCores
     !!{
     A property extractor class for the velocity dispersion at a set of radii.
     !!}
     private
     integer  :: radiusNuclearStarClustersID     , blackHoleSeedMassID              , &
         &       velocityNuclearStarClustersID   , ageNuclearStarClustersID         , &
         &       gasMassNuclearStarClustersID    , criticalMassNuclearStarClustersID, &
         &       stellarMassNuclearStarClustersID                                                      
   contains
     procedure :: extract      => blackHoleSeedingDarkCoresExtract
     procedure :: name         => blackHoleSeedingDarkCoresName
     procedure :: description  => blackHoleSeedingDarkCoresDescription
     procedure :: unitsInSI    => blackHoleSeedingDarkCoresUnitsInSI
  end type nodePropertyExtractorBlackHoleSeedingDarkCores

  interface nodePropertyExtractorBlackHoleSeedingDarkCores
     !!{
     Constructors for the ``BlackHoleSeedingDarkCores'' output analysis class.
     !!}
     module procedure blackHoleSeedingDarkCoresConstructorParameters
     module procedure blackHoleSeedingDarkCoresConstructorInternal
  end interface nodePropertyExtractorBlackHoleSeedingDarkCores

contains

  function blackHoleSeedingDarkCoresConstructorParameters(parameters) result(self)
    !!{
    Constructor for the {\normalfont \ttfamily BlackHoleSeedingDarkCores} property extractor class.
    !!}
    use :: Input_Parameters, only : inputParameters
    implicit none
    type (nodePropertyExtractorBlackHoleSeedingDarkCores)                :: self
    type (inputParameters                               ), intent(inout) :: parameters

    self=nodePropertyExtractorBlackHoleSeedingDarkCores()
    !![
    <inputParametersValidate source="parameters"/>
    !!]
    return
  end function blackHoleSeedingDarkCoresConstructorParameters

  function blackHoleSeedingDarkCoresConstructorInternal() result(self)
    !!{
    Internal constructor for the {\normalfont \ttfamily BlackHoleSeedingDarkCores} property extractor class.
    !!}
    implicit none
    type          (nodePropertyExtractorBlackHoleSeedingDarkCores) :: self
    !![
    <addMetaProperty component="NSC" name="blackHoleSeedMassFormed" id="self%blackHoleSeedMassID" isEvolvable="no" isCreator="no"/>
    !!]
    return
  end function blackHoleSeedingDarkCoresConstructorInternal

  function blackHoleSeedingDarkCoresExtract(self,node,instance)
    !!{
    Implement a {\normalfont \ttfamily BlackHoleSeedingDarkCores} property extractor.
    !!}
    use :: Galacticus_Nodes, only : nodeComponentNSC
    implicit none
    double precision                                                                             :: blackHoleSeedingDarkCoresExtract
    class           (nodePropertyExtractorBlackHoleSeedingDarkCores), intent(inout), target      :: self
    type            (treeNode                                      ), intent(inout), target      :: node
    type            (multiCounter                                  ), intent(inout), optional    :: instance
    class           (nodeComponentNSC                              )               , pointer     :: nuclearStarCluster
    !$GLC attributes unused :: instance

    nuclearStarCluster => node%NSC()
    select type (nuclearStarCluster)
    type is (nodeComponentNSC)
      ! Nuclear star cluster does not yet exist.
      blackHoleSeedingDarkCoresExtract=0.0d0
    class default
      blackHoleSeedingDarkCoresExtract=nuclearStarCluster%floatRank0MetaPropertyGet(self%blackHoleSeedMassID)
    end select
    return
  end function blackHoleSeedingDarkCoresExtract

  function blackHoleSeedingDarkCoresName(self)
    !!{
    Return the names of the {\normalfont \ttfamily BlackHoleSeedingDarkCores} properties.
    !!}
    implicit none
    class (nodePropertyExtractorBlackHoleSeedingDarkCores), intent(inout) :: self
    type  (varying_string                                )                :: blackHoleSeedingDarkCoresName
    !$GLC attributes unused :: self
    
    blackHoleSeedingDarkCoresName=var_str('blackHoleSeedMass')
    return
  end function blackHoleSeedingDarkCoresName

  function blackHoleSeedingDarkCoresDescription(self)
    !!{
    Return descriptions of the {\normalfont \ttfamily BlackHoleSeedingDarkCores} property.
    !!}
    implicit none
    class (nodePropertyExtractorBlackHoleSeedingDarkCores), intent(inout) :: self
    type  (varying_string                                )                :: blackHoleSeedingDarkCoresDescription

    blackHoleSeedingDarkCoresDescription=var_str('Black hole seed mass [M⊙].')
    return
  end function blackHoleSeedingDarkCoresDescription

  double precision function blackHoleSeedingDarkCoresUnitsInSI(self)
    !!{
    Return the units of the {\normalfont \ttfamily BlackHoleSeedingDarkCores} properties in the SI system.
    !!}
    use :: Numerical_Constants_Astronomical, only : massSolar
    implicit none
    class (nodePropertyExtractorBlackHoleSeedingDarkCores), intent(inout):: self
    !$GLC attributes unused :: self
    blackHoleSeedingDarkCoresUnitsInSI = massSolar
    return
  end function blackHoleSeedingDarkCoresUnitsInSI
