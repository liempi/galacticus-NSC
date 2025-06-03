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
  Implements a property extractor class the properties of nuclear star cluster when a black hole seed is formed using the model of \cite{vergara_global_2023}.
  !!}
  
  !![
  <nodePropertyExtractor name="nodePropertyExtractorDarkCores">
   <description>
    A property extractor class for the properties of the nuclear star cluster at the moment of the black hole formation.
   </description>
  </nodePropertyExtractor>
  !!]
  type, extends(nodePropertyExtractorTuple) :: nodePropertyExtractorDarkCores
     !!{
     A property extractor class for the velocity dispersion at a set of radii.
     !!}
     private
    integer   :: darkCoreRadiusID                 , darkCoreGasMassID          , &
       &         darkCoreVelocityDispersionID     , darkCoreTimescaleID        , &
       &         nuclearStarClusterNumberOfStarsID, nuclearStarClusterDensityID
   
   contains
     procedure :: elementCount       => darkCoresElementCount
     procedure :: extract            => darkCoresExtract
     procedure :: names              => darkCoresNames
     procedure :: descriptions       => darkCoresDescriptions
     procedure :: unitsInSI          => darkCoresUnitsInSI
  end type nodePropertyExtractorDarkCores

  interface nodePropertyExtractorDarkCores
     !!{
     Constructors for the ``DarkCores'' output analysis class.
     !!}
     module procedure darkCoresConstructorParameters
     module procedure darkCoresConstructorInternal
  end interface nodePropertyExtractorDarkCores

contains

  function darkCoresConstructorParameters(parameters) result(self)
    !!{
    Constructor for the {\normalfont \ttfamily DarkCores} property extractor class.
    !!}
    use :: Input_Parameters, only : inputParameters
    implicit none
    type (nodePropertyExtractorDarkCores)                :: self
    type (inputParameters               ), intent(inout) :: parameters

    self=nodePropertyExtractorDarkCores()
    !![
    <inputParametersValidate source="parameters"/>
    !!]
    return
  end function darkCoresConstructorParameters

  function darkCoresConstructorInternal() result(self)
    !!{
    Internal constructor for the {\normalfont \ttfamily DarkCores} property extractor class.
    !!}
    implicit none
    type          (nodePropertyExtractorDarkCores) :: self
    !![
    <addMetaProperty component="NSC" name="darkCoreRadius"                  id="self%darkCoreRadiusID"                  isEvolvable="no" isCreator="no"/>
    <addMetaProperty component="NSC" name="darkCoreGasMass"                 id="self%darkCoreGasMassID"                 isEvolvable="no" isCreator="no"/>
    <addMetaProperty component="NSC" name="darkCoreVelocityDispersion"      id="self%darkCoreVelocityDispersionID"      isEvolvable="no" isCreator="no"/>
    <addMetaProperty component="NSC" name="darkCoreTimescale"               id="self%darkCoreTimeScaleID"               isEvolvable="no" isCreator="no"/>
    <addMetaProperty component="NSC" name="nuclearStarClusterNumberOfStars" id="self%nuclearStarClusterNumberOfStarsID" isEvolvable="no" isCreator="no"/>
    <addMetaProperty component="NSC" name="nuclearStarClusterDensity"       id="self%nuclearStarClusterDensityID"       isEvolvable="no" isCreator="no"/> 
    !!]
    return
  end function darkCoresConstructorInternal

  integer function darkCoresElementCount(self,time)
    !!{
    Return the number of elements in the {\normalfont \ttfamily darkCores} property extractors.
    !!}
    implicit none
    class           (nodePropertyExtractorDarkCores), intent(inout) :: self
    double precision                                , intent(in   ) :: time
    !$GLC attributes unused :: time

    darkCoresElementCount=6
    return
  end function darkCoresElementCount

  function darkCoresExtract(self,node,time,instance)
    !!{
    Implement a {\normalfont \ttfamily darkCores} property extractor.
    !!}
    use :: Galacticus_Nodes, only : nodeComponentNSC
    implicit none
    double precision                                , dimension(:) , allocatable :: darkCoresExtract
    class           (nodePropertyExtractorDarkCores), intent(inout), target      :: self
    type            (treeNode                      ), intent(inout), target      :: node
    double precision                                , intent(in   )              :: time
    type            (multiCounter                  ), intent(inout), optional    :: instance
    class           (nodeComponentNSC              )               , pointer     :: nuclearStarCluster

    !$GLC attributes unused :: time, instance

    allocate(darkCoresExtract(6))
    nuclearStarCluster => node%NSC()
    select type (nuclearStarCluster)
    type is (nodeComponentNSC)
      ! Nuclear star cluster does not yet exist.
      darkCoresExtract=[ & 
        &                  0.0d0, &
        &                  0.0d0, &
        &                  0.0d0, &
        &                  0.0d0, &
        &                  0.0d0, &
        &                  0.0d0  &
        &              ]
    class default
      darkCoresExtract=[                                                                                       &
       &                 nuclearStarCluster%floatRank0MetaPropertyGet(self%darkCoreRadiusID                 ), &
       &                 nuclearStarCluster%floatRank0MetaPropertyGet(self%darkCoreGasMassID                ), &
       &                 nuclearStarCluster%floatRank0MetaPropertyGet(self%darkCoreVelocityDispersionID     ), &
       &                 nuclearStarCluster%floatRank0MetaPropertyGet(self%darkCoreTimescaleID              ), &
       &                 nuclearStarCluster%floatRank0MetaPropertyGet(self%nuclearStarClusterNumberOfStarsID), &
       &                 nuclearStarCluster%floatRank0MetaPropertyGet(self%nuclearStarClusterDensityID)        &
       &               ]
    end select
    return
  end function darkCoresExtract

  subroutine darkCoresNames(self,time,names)
    !!{
    Return the names of the {\normalfont \ttfamily DarkCores} properties.
    !!}
    implicit none
    class           (nodePropertyExtractorDarkCores), intent(inout)                             :: self
    double precision                                , intent(in   )                             :: time
    type            (varying_string                ), intent(inout), dimension(:) , allocatable :: names
    !$GLC attributes unused :: self, time
    
    allocate(names(6))
    names(1)=var_str('darkCoreRadius'                 )
    names(2)=var_str('darkCoreGasMass'                )
    names(3)=var_str('darkCoreVelocityDispersion'     )
    names(4)=var_str('darkCoreTimescale'              )
    names(5)=var_str('nuclearStarClusterNumberOfStars')
    names(6)=var_str('nuclearStarClusterDensityGas'   )

    return
  end subroutine darkCoresNames

  subroutine darkCoresDescriptions(self,time,descriptions)
    !!{
    Return descriptions of the {\normalfont \ttfamily DarkCores} property.
    !!}
    implicit none
    class           (nodePropertyExtractorDarkCores), intent(inout)                             :: self
    double precision                                                  , intent(in   )                             :: time
    type            (varying_string                                  ), intent(inout), dimension(:) , allocatable :: descriptions
    !$GLC attributes unused :: time

    allocate(descriptions(6))
    descriptions(1)=var_str('Radius of the dark core [Mpc].'                         )
    descriptions(2)=var_str('Gas mass of the nuclear star cluster enclose in the dark core radius [M⊙].'                                                )
    descriptions(3)=var_str('Velocity dispersion of the dark core[km/s].'                      )
    descriptions(4)=var_str('Dynamical friction timescale of the nuclear star cluster [Gyr].')
    descriptions(5)=var_str('Number of stars of the nuclear star cluster.')
    descriptions(6)=var_str('Gas density of the nuclear star cluster at the radius [M⊙ Mpc-3].')
    return
  end subroutine darkCoresDescriptions

  function darkCoresUnitsInSI(self,time)
    !!{
    Return the units of the {\normalfont \ttfamily DarkCores} properties in the SI system.
    !!}
    use :: Numerical_Constants_Astronomical, only : massSolar, megaParsec, gigayear
    use :: Numerical_Constants_Prefixes    , only : kilo
    implicit none
    double precision                                , allocatable  , dimension(:) :: darkCoresUnitsInSI
    class           (nodePropertyExtractorDarkCores), intent(inout)               :: self
    double precision                                , intent(in   )               :: time
    !$GLC attributes unused :: time

    allocate(darkCoresUnitsInSI(6))
    DarkCoresUnitsInSI=[                       &
     &                  megaParsec           , &
     &                  massSolar            , &
     &                  kilo                 , &
     &                  gigayear             , &
     &                  1.0d0                , &
     &                  1.0d0/megaParsec**3d0  &
     &                 ]
    return
  end function darkCoresUnitsInSI
