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
Implements a star formation rate property extractor class.
!!}

  use :: Globular_Cluster_Formation_Rates_Disks    , only : globularClusterFormationRateDisksClass
  use :: Globular_Cluster_Formation_Rates_Spheroids, only : globularClusterFormationRateSpheroidsClass

  !![
  <nodePropertyExtractor name="nodePropertyExtractorGlobularClusterFormationRate">
   <description>
    A node property extractor which extracts the globular cluster formation rate in a galaxy. The type of star formation rate is controlled by
    the {\normalfont \ttfamily [component]} parameter, which can be either ``{\normalfont \ttfamily disk}'', ``{\normalfont
    \ttfamily spheroid}'' or ``{\normalfont \ttfamily total}''. The corresponding globular cluster formation
    rate is extracted as {\normalfont \ttfamily \textless\ component\textgreater\ globularClusterFormationRate} in units of $M_\odot$/Gyr.
   </description>
  </nodePropertyExtractor>
  !!]
  type, extends(nodePropertyExtractorScalar) :: nodePropertyExtractorGlobularClusterFormationRate
     !!{
     A star formation rate property extractor class.
     !!}
     private
     class(globularClusterFormationRateDisksClass    ), pointer :: globularClusterFormationRateDisks_    => null()
     class(globularClusterFormationRateSpheroidsClass), pointer :: globularClusterFormationRateSpheroids_=> null()
     type (varying_string                            )          :: name_                                          , description_, &
          &                                                        component
   contains
     final     ::                globularClusterFormationRateDestructor
     procedure :: extract     => globularClusterFormationRateExtract
     procedure :: name        => globularClusterFormationRateName
     procedure :: description => globularClusterFormationRateDescription
     procedure :: unitsInSI   => globularClusterFormationRateUnitsInSI
  end type nodePropertyExtractorGlobularClusterFormationRate

  interface nodePropertyExtractorGlobularClusterFormationRate
     !!{
     Constructors for the {\normalfont \ttfamily globularClusterFormationRate} output analysis class.
     !!}
     module procedure globularClusterFormationRateConstructorParameters
     module procedure globularClusterFormationRateConstructorInternal
  end interface nodePropertyExtractorGlobularClusterFormationRate

contains

  function globularClusterFormationRateConstructorParameters(parameters) result(self)
    !!{
    Constructor for the {\normalfont \ttfamily globularClusterFormationRate} property extractor class which takes a parameter set as input.
    !!}
    use :: Input_Parameters, only : inputParameter, inputParameters
    implicit none
    type (nodePropertyExtractorGlobularClusterFormationRate)                :: self
    type (inputParameters                                  ), intent(inout) :: parameters
    class(globularClusterFormationRateDisksClass           ), pointer       :: globularClusterFormationRateDisks_
    class(globularClusterFormationRateSpheroidsClass       ), pointer       :: globularClusterFormationRateSpheroids_
    type (varying_string                                   )                :: component
    type (enumerationGalacticComponentType                 )                :: component_

    !![
    <inputParameter>
      <name>component</name>
      <source>parameters</source>
      <description>The component from which to extract star formation rate.</description>
    </inputParameter>
    !!]
    component_=enumerationGalacticComponentEncode(char(component),includesPrefix=.false.)
    select case (component_%ID)
    case (galacticComponentDisk              %ID)
       !![
       <objectBuilder class="globularClusterFormationRateDisks"               name="globularClusterFormationRateDisks_"               source="parameters"/>
       !!]
       globularClusterFormationRateSpheroids_           => null()
    case (galacticComponentSpheroid          %ID)
       !![
       <objectBuilder class="globularClusterFormationRateSpheroids"           name="globularClusterFormationRateSpheroids_"           source="parameters"/>
       !!]
       globularClusterFormationRateDisks_               => null()
    case (galacticComponentTotal            %ID)
       !![
       <objectBuilder class="globularClusterFormationRateDisks"               name="globularClusterFormationRateDisks_"               source="parameters"/>
       <objectBuilder class="globularClusterFormationRateSpheroids"           name="globularClusterFormationRateSpheroids_"           source="parameters"/>
       !!]
    end select
    self=nodePropertyExtractorGlobularClusterFormationRate(globularClusterFormationRateDisks_,globularClusterFormationRateSpheroids_)
    !![
    <inputParametersValidate source="parameters"/>
    <objectDestructor name="globularClusterFormationRateDisks_"              />
    <objectDestructor name="globularClusterFormationRateSpheroids_"          />
    !!]
    return
  end function globularClusterFormationRateConstructorParameters

  function globularClusterFormationRateConstructorInternal(globularClusterFormationRateDisks_,globularClusterFormationRateSpheroids_) result(self)
    !!{
    Internal constructor for the {\normalfont \ttfamily globularClusterFormationRate} property extractor class.
    !!}
    use :: Error, only : Error_Report
    implicit none
    type (nodePropertyExtractorGlobularClusterFormationRate   )                        :: self
    class(globularClusterFormationRateDisksClass              ), intent(in   ), target :: globularClusterFormationRateDisks_
    class(globularClusterFormationRateSpheroidsClass          ), intent(in   ), target :: globularClusterFormationRateSpheroids_
    !![
    <constructorAssign variables="*globularClusterFormationRateDisks_, *globularClusterFormationRateSpheroids_"/>
    !!]

    if      (associated(self%globularClusterFormationRateDisks_).and.associated(self%globularClusterFormationRateSpheroids_)) then
       self%name_       ="totalglobularClusterFormationRate"
       self%description_="Total (disk + spheroid) star formation rate [M☉ Gyr⁻¹]."
       self%component   ="total"
    else if (associated(self%globularClusterFormationRateDisks_)                                                                                                            ) then
       self%name_       ="diskglobularClusterFormationRate"
       self%description_="Disk globular cluster formation rate [M☉ Gyr⁻¹]."
       self%component   ="disk"
    else if (                                                        associated(self%globularClusterFormationRateSpheroids_)                                                           ) then
       self%name_       ="spheroidglobularClusterFormationRate"
       self%description_="Spheroid globular cluster formation rate [M☉ Gyr⁻¹]."
       self%component   ="spheroid"
    else
       call Error_Report('No star formation rate specified.'//{introspection:location})
    end if
    return
  end function globularClusterFormationRateConstructorInternal

  subroutine globularClusterFormationRateDestructor(self)
    !!{
    Destructor for the {\normalfont \ttfamily globularClusterFormationRate} property extractor class.
    !!}
    implicit none
    type(nodePropertyExtractorGlobularClusterFormationRate), intent(inout) :: self
  
    !![
    <objectDestructor name="self%globularClusterFormationRateDisks_"              />
    <objectDestructor name="self%globularClusterFormationRateSpheroids_"          />
    !!]
    return
  end subroutine globularClusterFormationRateDestructor

  double precision function globularClusterFormationRateExtract(self,node,instance)
    !!{
    Implement a star formation rate output analysis property extractor.
    !!}
    implicit none
    class(nodePropertyExtractorGlobularClusterFormationRate), intent(inout), target   :: self
    type (treeNode                                         ), intent(inout), target   :: node
    type (multiCounter                                     ), intent(inout), optional :: instance
    !$GLC attributes unused :: instance

    globularClusterFormationRateExtract=0.0d0
    if (associated(self%globularClusterFormationRateDisks_              )) globularClusterFormationRateExtract=globularClusterFormationRateExtract+self%globularClusterFormationRateDisks_              %rate(node)
    if (associated(self%globularClusterFormationRateSpheroids_          )) globularClusterFormationRateExtract=globularClusterFormationRateExtract+self%globularClusterFormationRateSpheroids_          %rate(node)
    return
  end function globularClusterFormationRateExtract

  function globularClusterFormationRateName(self)
    !!{
    Return the name of the globularClusterFormationRate property.
    !!}
    implicit none
    type (varying_string                                   )                :: globularClusterFormationRateName
    class(nodePropertyExtractorGlobularClusterFormationRate), intent(inout) :: self

    globularClusterFormationRateName=self%name_
    return
  end function globularClusterFormationRateName

  function globularClusterFormationRateDescription(self)
    !!{
    Return a description of the globularClusterFormationRate property.
    !!}
    implicit none
    type (varying_string                                   )                :: globularClusterFormationRateDescription
    class(nodePropertyExtractorGlobularClusterFormationRate), intent(inout) :: self

    globularClusterFormationRateDescription=self%description_
    return
  end function globularClusterFormationRateDescription

  double precision function globularClusterFormationRateUnitsInSI(self)
    !!{
    Return the units of the globularClusterFormationRate property in the SI system.
    !!}
    use :: Numerical_Constants_Astronomical, only : massSolar, gigaYear
    implicit none
    class(nodePropertyExtractorGlobularClusterFormationRate), intent(inout) :: self
    !$GLC attributes unused :: self

    globularClusterFormationRateUnitsInSI=massSolar/gigaYear
    return
  end function globularClusterFormationRateUnitsInSI
