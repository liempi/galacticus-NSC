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
  type, extends(nodePropertyExtractorScalar) :: nodePropertyExtractorGlobularClusterEvaporationRate
     !!{
     A star Evaporation rate property extractor class.
     !!}
     private
     class(globularClusterEvaporationRateDisksClass    ), pointer :: globularClusterEvaporationRateDisks_    => null()
     class(globularClusterEvaporationRateSpheroidsClass), pointer :: globularClusterEvaporationRateSpheroids_=> null()
     type (varying_string                              )          :: name_                                            , description_, &
          &                                                          component
   contains
     final     ::                globularClusterEvaporationRateDestructor
     procedure :: extract     => globularClusterEvaporationRateExtract
     procedure :: name        => globularClusterEvaporationRateName
     procedure :: description => globularClusterEvaporationRateDescription
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
    type (varying_string                                     )                :: component
    type (enumerationGalacticComponentType                   )                :: component_

    !![
    <inputParameter>
      <name>component</name>
      <source>parameters</source>
      <description>The component from which to extract globular cluster evaporation rate.</description>
    </inputParameter>
    !!]
    component_=enumerationGalacticComponentEncode(char(component),includesPrefix=.false.)
    select case (component_%ID)
    case (galacticComponentDisk    %ID)
       !![
       <objectBuilder class="globularClusterEvaporationRateDisks"     name="globularClusterEvaporationRateDisks_"               source="parameters"/>
       !!]
       globularClusterEvaporationRateSpheroids_ => null()
    case (galacticComponentSpheroid%ID)
       !![
       <objectBuilder class="globularClusterEvaporationRateSpheroids" name="globularClusterEvaporationRateSpheroids_"           source="parameters"/>
       !!]
       globularClusterEvaporationRateDisks_     => null()
    case (galacticComponentTotal   %ID)
       !![
       <objectBuilder class="globularClusterEvaporationRateDisks"     name="globularClusterEvaporationRateDisks_"               source="parameters"/>
       <objectBuilder class="globularClusterEvaporationRateSpheroids" name="globularClusterEvaporationRateSpheroids_"           source="parameters"/>
       !!]
    end select
    self=nodePropertyExtractorGlobularClusterEvaporationRate(globularClusterEvaporationRateDisks_,globularClusterEvaporationRateSpheroids_)
    !![
    <inputParametersValidate source="parameters"/>
    <objectDestructor name="globularClusterEvaporationRateDisks_"    />
    <objectDestructor name="globularClusterEvaporationRateSpheroids_"/>
    !!]
    return
  end function globularClusterEvaporationRateConstructorParameters

  function globularClusterEvaporationRateConstructorInternal(globularClusterEvaporationRateDisks_,globularClusterEvaporationRateSpheroids_) result(self)
    !!{
    Internal constructor for the {\normalfont \ttfamily globularClusterEvaporationRate} property extractor class.
    !!}
    use :: Error, only : Error_Report
    implicit none
    type (nodePropertyExtractorGlobularClusterEvaporationRate)                        :: self
    class(globularClusterEvaporationRateDisksClass           ), intent(in   ), target :: globularClusterEvaporationRateDisks_
    class(globularClusterEvaporationRateSpheroidsClass       ), intent(in   ), target :: globularClusterEvaporationRateSpheroids_
    !![
    <constructorAssign variables="*globularClusterEvaporationRateDisks_, *globularClusterEvaporationRateSpheroids_"/>
    !!]

    if      (associated(self%globularClusterEvaporationRateDisks_).and.associated(self%globularClusterEvaporationRateSpheroids_)) then
       self%name_       ="totalglobularClusterEvaporationRate"
       self%description_="Total (disk + spheroid) globular cluster evaporation rate [M☉ Gyr⁻¹]."
       self%component   ="total"
    else if (associated(self%globularClusterEvaporationRateDisks_)                                                                                                            ) then
       self%name_       ="diskglobularClusterEvaporationRate"
       self%description_="Disk globular cluster evaporation rate [M☉ Gyr⁻¹]."
       self%component   ="disk"
    else if (                                                          associated(self%globularClusterEvaporationRateSpheroids_)                                                           ) then
       self%name_       ="spheroidglobularClusterEvaporationRate"
       self%description_="Spheroid globular cluster evaporation rate [M☉ Gyr⁻¹]."
       self%component   ="spheroid"
    else
       call Error_Report('No globular clusuter evaporation rate specified.'//{introspection:location})
    end if
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

  double precision function globularClusterEvaporationRateExtract(self,node,instance)
    !!{
    Implement a star Evaporation rate output analysis property extractor.
    !!}
    implicit none
    class(nodePropertyExtractorGlobularClusterEvaporationRate), intent(inout), target   :: self
    type (treeNode                                           ), intent(inout), target   :: node
    type (multiCounter                                       ), intent(inout), optional :: instance
    !$GLC attributes unused :: instance

    globularClusterEvaporationRateExtract=0.0d0
    if (associated(self%globularClusterEvaporationRateDisks_    )) globularClusterEvaporationRateExtract=globularClusterEvaporationRateExtract+self%globularClusterEvaporationRateDisks_    %rate(node)
    if (associated(self%globularClusterEvaporationRateSpheroids_)) globularClusterEvaporationRateExtract=globularClusterEvaporationRateExtract+self%globularClusterEvaporationRateSpheroids_%rate(node)
    return
  end function globularClusterEvaporationRateExtract

  function globularClusterEvaporationRateName(self)
    !!{
    Return the name of the globularClusterEvaporationRate property.
    !!}
    implicit none
    type (varying_string                                     )                :: globularClusterEvaporationRateName
    class(nodePropertyExtractorGlobularClusterEvaporationRate), intent(inout) :: self

    globularClusterEvaporationRateName=self%name_
    return
  end function globularClusterEvaporationRateName

  function globularClusterEvaporationRateDescription(self)
    !!{
    Return a description of the globularClusterEvaporationRate property.
    !!}
    implicit none
    type (varying_string                                     )                :: globularClusterEvaporationRateDescription
    class(nodePropertyExtractorGlobularClusterEvaporationRate), intent(inout) :: self

    globularClusterEvaporationRateDescription=self%description_
    return
  end function globularClusterEvaporationRateDescription

  double precision function globularClusterEvaporationRateUnitsInSI(self)
    !!{
    Return the units of the globularClusterEvaporationRate property in the SI system.
    !!}
    use :: Numerical_Constants_Astronomical, only : massSolar, gigaYear
    implicit none
    class(nodePropertyExtractorGlobularClusterEvaporationRate), intent(inout) :: self
    !$GLC attributes unused :: self

    globularClusterEvaporationRateUnitsInSI=massSolar/gigaYear
    return
  end function globularClusterEvaporationRateUnitsInSI
