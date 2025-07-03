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

  !![
  <nodePropertyExtractor name="nodePropertyExtractorGlobularClusterMass">
   <description>
     A node property extractor which extracts stellar mass for disk and spheroid components.
   </description>
  </nodePropertyExtractor>
  !!]
  type, extends(nodePropertyExtractorTuple) :: nodePropertyExtractorGlobularClusterMass
     !!{
     A property extractor which extracts globular cluster mass of disk and spheroid.
     !!}
     private
     integer   :: globularClusterStellarMassDiskID, globularClusterStellarMassSpheroidID
   contains
     procedure :: elementCount => globularClusterMassElementCount
     procedure :: extract      => globularClusterMassExtract
     procedure :: names        => globularClusterMassNames
     procedure :: descriptions => globularClusterMassDescriptions
     procedure :: unitsInSI    => globularClusterMassUnitsInSI
  end type nodePropertyExtractorGlobularClusterMass

  interface nodePropertyExtractorGlobularClusterMass
     !!{
     Constructors for the ``GlobularClusterMass'' output extractor class.
     !!}
     module procedure globularClusterMassConstructorParameters
     module procedure globularClusterMassConstructorInternal
  end interface nodePropertyExtractorGlobularClusterMass

contains

  function globularClusterMassConstructorParameters(parameters) result(self)
    !!{
    Constructor for the ``GlobularClusterMass'' property extractor class which takes a parameter set as input.
    !!}
    use :: Input_Parameters, only : inputParameter, inputParameters
    implicit none
    type (nodePropertyExtractorGlobularClusterMass)                :: self
    type (inputParameters                         ), intent(inout) :: parameters

    self=nodePropertyExtractorGlobularClusterMass()
    !![
    <inputParametersValidate source="parameters"/>
    !!]
    return
  end function globularClusterMassConstructorParameters

  function globularClusterMassConstructorInternal() result(self)
    !!{
    Internal constructor for the ``GlobularClusterMass'' output extractor property extractor class.
    !!}
    implicit none
    type(nodePropertyExtractorGlobularClusterMass) :: self
    
    !![
    <addMetaProperty component="disk"     name="globularClusterStellarMassDisk"     id="self%globularClusterStellarMassDiskID"     isEvolvable="yes" isCreator="no" />
    <addMetaProperty component="spheroid" name="globularClusterStellarMassSpheroid" id="self%globularClusterStellarMassSpheroidID" isEvolvable="yes" isCreator="no" />
    !!]
    return
  end function globularClusterMassConstructorInternal

  integer function globularClusterMassElementCount(self,time)
    !!{
    Return the number of elements in the {\normalfont \ttfamily GlobularClusterMass} property extractors.
    !!}
    implicit none
    class           (nodePropertyExtractorGlobularClusterMass), intent(inout) :: self
    double precision                                          , intent(in   ) :: time
    !$GLC attributes unused :: self, time

    globularClusterMassElementCount=2
    return
  end function globularClusterMassElementCount

  function globularClusterMassExtract(self,node,time,instance)
    !!{
    Implement a GlobularClusterMass output extractor.
    !!}
    use :: Galacticus_Nodes, only : nodeComponentDisk, nodeComponentSpheroid
    implicit none
    double precision                                          , dimension(:) , allocatable :: globularClusterMassExtract
    class           (nodePropertyExtractorGlobularClusterMass), intent(inout), target      :: self
    type            (treeNode                                ), intent(inout), target      :: node
    double precision                                          , intent(in   )              :: time
    type            (multiCounter                            ), intent(inout), optional    :: instance
    class           (nodeComponentDisk                       )               , pointer     :: disk
    class           (nodeComponentSpheroid                   )               , pointer     :: spheroid
    double precision                                                                       :: massStellarDisk           , massStellarSpheroid
            
    !$GLC attributes unused :: time, instance

    ! Extract required quantities.
    disk               => node%disk    ()
    spheroid           => node%spheroid()

    select type (disk)
    type is (nodeComponentDisk    )
       ! Disk does not yet exist.
       massStellarDisk                   =0.0d0
    class default
       massStellarDisk                   =disk    %floatRank0MetaPropertyGet(self%globularClusterStellarMassDiskID    )
    end select
    select type (spheroid)
    type is (nodeComponentSpheroid)
       ! Spheroid does not yet exist.
       massStellarSpheroid               =0.0d0
    class default
       massStellarSpheroid               =spheroid%floatRank0MetaPropertyGet(self%globularClusterStellarMassSpheroidID)
    end select   

    ! Set return results.
    allocate(globularClusterMassExtract(2))
    globularClusterMassExtract=[                      &
         &                       massStellarDisk    , &
         &                       massStellarSpheroid  &
         &                     ]
    return
  end function globularClusterMassExtract

  subroutine globularClusterMassNames(self,time,names)
    !!{
    Return the names of the {\normalfont \ttfamily GlobularClusterMass} properties.
    !!}
    implicit none
    class           (nodePropertyExtractorGlobularClusterMass), intent(inout)                             :: self
    double precision                                          , intent(in   )                             :: time
    type            (varying_string                          ), intent(inout), dimension(:) , allocatable :: names
    !$GLC attributes unused :: self, time

    allocate(names(2))
    names(1)=var_str('diskGlobularClusterMass'    )
    names(2)=var_str('spheroidGlobularClusterMass')
    return
  end subroutine globularClusterMassNames

  subroutine globularClusterMassDescriptions(self,time,descriptions)
    !!{
    Return the descriptions of the {\normalfont \ttfamily GlobularClusterMass} properties.
    !!}
    implicit none
    class           (nodePropertyExtractorGlobularClusterMass), intent(inout)                             :: self
    double precision                                          , intent(in   )                             :: time
    type            (varying_string                          ), intent(inout), dimension(:) , allocatable :: descriptions
    !$GLC attributes unused :: self, time

    allocate(descriptions(2))
    descriptions(1)=var_str('Stellar mass of globular clusters in the disk [M☉].'    )
    descriptions(2)=var_str('Stellar mass of globular clusters in the spheroid [M☉].')
    return
  end subroutine globularClusterMassDescriptions

  function globularClusterMassUnitsInSI(self,time)
    !!{
    Return the units of the {\normalfont \ttfamily GlobularClusterMass} properties in the SI system.
    !!}
    use :: Numerical_Constants_Astronomical, only : massSolar
    implicit none
    double precision                                          , dimension(:) , allocatable :: globularClusterMassUnitsInSI
    class           (nodePropertyExtractorGlobularClusterMass), intent(inout)              :: self
    double precision                                          , intent(in   )              :: time
   !$GLC attributes unused :: self, time

    allocate(globularClusterMassUnitsInSI(2))
    globularClusterMassUnitsInSI=[massSolar, massSolar]
    return
  end function globularClusterMassUnitsInSI

