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
Implements a star Evaporation rate property extractor class.
!!}

  use :: Globular_Cluster_Evaporation_Rates_Disks    , only : globularClusterEvaporationRateDisksClass
  use :: Globular_Cluster_Evaporation_Rates_Spheroids, only : globularClusterEvaporationRateSpheroidsClass

  !![
  <nodePropertyExtractor name="nodePropertyExtractorGlobularClusterEvaporationRate">
   <description>
    A node property extractor which extracts the globular cluster evaporation rate in a galaxy. The type globular evaporation rate is controlled by
    the {\normalfont \ttfamily [component]} parameter, which can be either ``{\normalfont \ttfamily disk}'', ``{\normalfont
    \ttfamily spheroid}'' or ``{\normalfont \ttfamily total}''. The corresponding globular cluster evaporation
    rate is extracted as {\normalfont \ttfamily \textless\ component\textgreater\ globularClusterEvaporationRate} in units of $M_\odot$/Gyr.
   </description>
  </nodePropertyExtractor>
  !!]
  type, extends(nodePropertyExtractorTuple) :: nodePropertyExtractorGlobularClusterEvaporationRate
     !!{
     A star Evaporation rate property extractor class.
     !!}
     private
     class(globularClusterEvaporationRateDisksClass    ), pointer :: globularClusterEvaporationRateDisks_    => null()
     class(globularClusterEvaporationRateSpheroidsClass), pointer :: globularClusterEvaporationRateSpheroids_=> null()
   contains
     final     ::                globularClusterEvaporationRateDestructor
     procedure :: elementCount=> globularClusterEvaporationRateElementCount
     procedure :: extract     => globularClusterEvaporationRateExtract
     procedure :: names       => globularClusterEvaporationRateNames
     procedure :: descriptions=> globularClusterEvaporationRateDescriptions
     procedure :: unitsInSI   => globularClusterEvaporationRateUnitsInSI
  end type nodePropertyExtractorGlobularClusterEvaporationRate

  interface nodePropertyExtractorGlobularClusterEvaporationRate
     !!{
     Constructors for the {\normalfont \ttfamily globularClusterEvaporationRate} output analysis class.
     !!}
     module procedure globularClusterEvaporationRateConstructorParameters
     module procedure globularClusterEvaporationRateConstructorInternal
  end interface nodePropertyExtractorGlobularClusterEvaporationRate

contains

  function globularClusterEvaporationRateConstructorParameters(parameters) result(self)
    !!{
    Constructor for the {\normalfont \ttfamily globularClusterEvaporationRate} property extractor class which takes a parameter set as input.
    !!}
    use :: Input_Parameters, only : inputParameter, inputParameters
    implicit none
    type (nodePropertyExtractorGlobularClusterEvaporationRate)                :: self
    type (inputParameters                                    ), intent(inout) :: parameters
    class(globularClusterEvaporationRateDisksClass           ), pointer       :: globularClusterEvaporationRateDisks_
    class(globularClusterEvaporationRateSpheroidsClass       ), pointer       :: globularClusterEvaporationRateSpheroids_

    !![
    <objectBuilder class="globularClusterEvaporationRateDisks"     name="globularClusterEvaporationRateDisks_"     source="parameters"/>
    <objectBuilder class="globularClusterEvaporationRateSpheroids" name="globularClusterEvaporationRateSpheroids_" source="parameters"/>
    !!]
    self=nodePropertyExtractorGlobularClusterEvaporationRate(globularClusterEvaporationRateDisks_,globularClusterEvaporationRateSpheroids_)
    !![
    <inputParametersValidate source="parameters"/>
    <objectDestructor name="globularClusterEvaporationRateDisks_"    />
    <objectDestructor name="globularClusterEvaporationRateSpheroids_"/>
    !!]
    return
  end function globularClusterEvaporationRateConstructorParameters

  integer function globularClusterEvaporationRateElementCount(self,time)
    !!{
    Return the number of elements in the {\normalfont \ttfamily radiiHalfLightProperties} property extractor class.
    !!}
    implicit none
    class           (nodePropertyExtractorGlobularClusterEvaporationRate), intent(inout) :: self
    double precision                                                     , intent(in   ) :: time
    !$GLC attributes unused :: self

    globularClusterEvaporationRateElementCount=2
    return
  end function globularClusterEvaporationRateElementCount

  function globularClusterEvaporationRateConstructorInternal(globularClusterEvaporationRateDisks_,globularClusterEvaporationRateSpheroids_) result(self)
    !!{
    Internal constructor for the {\normalfont \ttfamily globularClusterEvaporationRate} property extractor class.
    !!}
    implicit none
    type (nodePropertyExtractorGlobularClusterEvaporationRate)                        :: self
    class(globularClusterEvaporationRateDisksClass           ), intent(in   ), target :: globularClusterEvaporationRateDisks_
    class(globularClusterEvaporationRateSpheroidsClass       ), intent(in   ), target :: globularClusterEvaporationRateSpheroids_
    !![
    <constructorAssign variables="*globularClusterEvaporationRateDisks_, *globularClusterEvaporationRateSpheroids_"/>
    !!]
    return
  end function globularClusterEvaporationRateConstructorInternal

  subroutine globularClusterEvaporationRateDestructor(self)
    !!{
    Destructor for the {\normalfont \ttfamily globularClusterEvaporationRate} property extractor class.
    !!}
    implicit none
    type(nodePropertyExtractorGlobularClusterEvaporationRate), intent(inout) :: self
  
    !![
    <objectDestructor name="self%globularClusterEvaporationRateDisks_"    />
    <objectDestructor name="self%globularClusterEvaporationRateSpheroids_"/>
    !!]
    return
  end subroutine globularClusterEvaporationRateDestructor

  double precision function globularClusterEvaporationRateExtract(self,node,time,instance)
    !!{
    Implement a star Evaporation rate output analysis property extractor.
    !!}
    implicit none
    double precision                                                     , dimension(:) , allocatable :: globularClusterEvaporationRateExtract
    class           (nodePropertyExtractorGlobularClusterEvaporationRate), intent(inout), target      :: self
    type            (treeNode                                           ), intent(inout), target      :: node
    double precision                                                     , intent(in   )              :: time
    type            (multiCounter                                       ), intent(inout), optional    :: instance
    !$GLC attributes unused ::self, instance
    
    allocate(globularClusterEvaporationRateExtract(2))
    globularClusterEvaporationRateExtract(0)=+self%globularClusterEvaporationRateDisks_    %rate(node)
    globularClusterEvaporationRateExtract(1)=+self%globularClusterEvaporationRateSpheroids_%rate(node)
    return 
  end function globularClusterEvaporationRateExtract

  subroutine globularClusterEvaporationRateNames(self,time,names)
    !!{
    Return the name of the globularClusterEvaporationRate property.
    !!}
    implicit none
    class           (nodePropertyExtractorGlobularClusterEvaporationRate), intent(inout)                             :: self
    double precision                                                     , intent(in   )                             :: time
    type            (varying_string                                     ), intent(inout), dimension(:) , allocatable :: names
    !$GLC attributes unused :: self, time
    allocate(names(2))
    names(0)=var_str('diskGlobularClusterEvaporationRate'    )
    names(1)=var_str('spheroidGlobularClusterEvapotationRate')
    return
  end subroutine globularClusterEvaporationRateNames

  subroutine globularClusterEvaporationRateDescriptions(self,time,descriptions)
    !!{
    Return a description of the globularClusterEvaporationRate property.
    !!}
    implicit none
    class(nodePropertyExtractorGlobularClusterEvaporationRate), intent(inout)                             :: self
    double precision                                          , intent(in   )                             :: time
    type            (varying_string                          ), intent(inout), dimension(:) , allocatable :: descriptions
    !$GLC attributes unused :: self, time
    allocate(descriptions(2))
    descriptions(0)=var_str('Disk globular cluster evaporation rate [M☉ Gyr⁻¹].')
    descriptions(1)=var_str('Spheroidal globular cluster evaporation rate [M☉ Gyr⁻¹].')
    return
  end subroutine globularClusterEvaporationRateDescriptions

  double precision function globularClusterEvaporationRateUnitsInSI(self, time)
    !!{
    Return the units of the globularClusterEvaporationRate property in the SI system.
    !!}
    use :: Numerical_Constants_Astronomical, only : massSolar, gigaYear
    implicit none
    double precision                                                     , allocatable  , dimension(:) :: globularClusterEvaporationRateUnitsInSI
    class           (nodePropertyExtractorGlobularClusterEvaporationRate), intent(inout)               :: self
    double precision                                                     , intent(in   )               :: time
    !$GLC attributes unused :: self, time
    allocate(globularClusterEvaporationRateUnitsInSI(2))
    globularClusterEvaporationRateUnitsInSI=[massSolar/gigaYear, massSolar/gigaYear]
    return
  end function globularClusterEvaporationRateUnitsInSI
