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
    integer   :: darkCoreRadiusID            , darkCoreGasMassID  , &
       &         darkCoreVelocityDispersionID, darkCoreTimeScaleID
   
   contains
     procedure :: elementCount       => DarkCoresElementCount
     procedure :: extract            => DarkCoresExtract
     procedure :: names              => DarkCoresNames
     procedure :: descriptions       => DarkCoresDescriptions
     procedure :: unitsInSI          => DarkCoresUnitsInSI
  end type nodePropertyExtractorDarkCores

  interface nodePropertyExtractorDarkCores
     !!{
     Constructors for the ``DarkCores'' output analysis class.
     !!}
     module procedure DarkCoresConstructorParameters
     module procedure DarkCoresConstructorInternal
  end interface nodePropertyExtractorDarkCores

contains

  function DarkCoresConstructorParameters(parameters) result(self)
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
  end function DarkCoresConstructorParameters

  function DarkCoresConstructorInternal() result(self)
    !!{
    Internal constructor for the {\normalfont \ttfamily DarkCores} property extractor class.
    !!}
    implicit none
    type          (nodePropertyExtractorDarkCores) :: self
    !![
    <addMetaProperty component="NSC" name="darkCoreRadius"             id="self%darkCoreRadiusID"             isEvolvable="no" isCreator="no"/>
    <addMetaProperty component="NSC" name="darkCoreGasMass"            id="self%darkCoreGasMassID"            isEvolvable="no" isCreator="no"/>
    <addMetaProperty component="NSC" name="darkCoreVelocityDispersion" id="self%darkCoreVelocityDispersionID" isEvolvable="no" isCreator="no"/>
    <addMetaProperty component="NSC" name="darkCoreTimeScale"          id="self%darkCoreTimeScaleID"          isEvolvable="no" isCreator="no"/>
    !!]
    return
  end function DarkCoresConstructorInternal

  integer function DarkCoresElementCount(self,time)
    !!{
    Return the number of elements in the {\normalfont \ttfamily DarkCores} property extractors.
    !!}
    implicit none
    class           (nodePropertyExtractorDarkCores), intent(inout) :: self
    double precision                                , intent(in   ) :: time
    !$GLC attributes unused :: time

    DarkCoresElementCount=4
    return
  end function DarkCoresElementCount

  function DarkCoresExtract(self,node,time,instance)
    !!{
    Implement a {\normalfont \ttfamily DarkCores} property extractor.
    !!}
    use :: Galacticus_Nodes, only : nodeComponentNSC
    implicit none
    double precision                                , dimension(:) , allocatable :: DarkCoresExtract
    class           (nodePropertyExtractorDarkCores), intent(inout), target      :: self
    type            (treeNode                      ), intent(inout), target      :: node
    double precision                                , intent(in   )              :: time
    type            (multiCounter                  ), intent(inout), optional    :: instance
    class           (nodeComponentNSC              )               , pointer     :: nuclearStarCluster

    !$GLC attributes unused :: time, instance

    allocate(DarkCoresExtract(4))
    nuclearStarCluster => node%NSC()
    select type (nuclearStarCluster)
    type is (nodeComponentNSC)
      ! Nuclear star cluster does not yet exist.
      DarkCoresExtract=[ & 
        &                  0.0d0, &
        &                  0.0d0, &
        &                  0.0d0, &
        &                  0.0d0  &
        &              ]
    class default
      DarkCoresExtract=[                                                                                  &
       &                 nuclearStarCluster%floatRank0MetaPropertyGet(self%darkCoreRadiusID            ), &
       &                 nuclearStarCluster%floatRank0MetaPropertyGet(self%darkCoreGasMassID           ), &
       &                 nuclearStarCluster%floatRank0MetaPropertyGet(self%darkCoreVelocityDispersionID), &
       &                 nuclearStarCluster%floatRank0MetaPropertyGet(self%darkCoreTimeScaleID         )  &
       &               ]
    end select
    return
  end function DarkCoresExtract

  subroutine DarkCoresNames(self,time,names)
    !!{
    Return the names of the {\normalfont \ttfamily DarkCores} properties.
    !!}
    implicit none
    class           (nodePropertyExtractorDarkCores), intent(inout)                             :: self
    double precision                                                  , intent(in   )                             :: time
    type            (varying_string                                  ), intent(inout), dimension(:) , allocatable :: names
    !$GLC attributes unused :: self, time
    
    allocate(names(4))
    names(1)=var_str('darkCoreRadius'            )
    names(2)=var_str('darkCoreGasMass'           )
    names(3)=var_str('darkCoreVelocityDispersion')
    names(4)=var_str('darkCoreTimescale'         )
    return
  end subroutine DarkCoresNames

  subroutine DarkCoresDescriptions(self,time,descriptions)
    !!{
    Return descriptions of the {\normalfont \ttfamily DarkCores} property.
    !!}
    implicit none
    class           (nodePropertyExtractorDarkCores), intent(inout)                             :: self
    double precision                                                  , intent(in   )                             :: time
    type            (varying_string                                  ), intent(inout), dimension(:) , allocatable :: descriptions
    !$GLC attributes unused :: time

    allocate(descriptions(4))
    descriptions(1)=var_str('Radius of the dark core [Mpc].'                         )
    descriptions(2)=var_str('Gas mass of the nuclear star cluster enclose in the dark core radius [M⊙].'                                                )
    descriptions(3)=var_str('Velocity dispersion of the dark core[km/s].'                      )
    descriptions(4)=var_str('Dynamical friction timescale of the nuclear star cluster [Gyr].')
    return
  end subroutine DarkCoresDescriptions

  function DarkCoresUnitsInSI(self,time)
    !!{
    Return the units of the {\normalfont \ttfamily DarkCores} properties in the SI system.
    !!}
    use :: Numerical_Constants_Astronomical, only : massSolar, megaParsec, gigayear
    use :: Numerical_Constants_Prefixes    , only : kilo
    implicit none
    double precision                                , allocatable  , dimension(:) :: DarkCoresUnitsInSI
    class           (nodePropertyExtractorDarkCores), intent(inout)               :: self
    double precision                                , intent(in   )               :: time
    !$GLC attributes unused :: time

    allocate(DarkCoresUnitsInSI(4))
    DarkCoresUnitsInSI = [  &
     &                        megaParsec, &
     &                        massSolar , &
     &                        kilo      , &
     &                        gigayear    &
     &                    ]
    return
  end function DarkCoresUnitsInSI
