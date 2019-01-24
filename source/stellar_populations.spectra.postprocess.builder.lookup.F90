!! Copyright 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018,
!!           2019
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

  !% Implements a stellar population spectra postprocessor builder which simply looks up postprocessors by name.
  
  !# <stellarPopulationSpectraPostprocessorBuilder name="stellarPopulationSpectraPostprocessorBuilderLookup">
  !#  <description>A stellar population spectra postprocessor builder which simply looks up postprocessors by name.</description>
  !# </stellarPopulationSpectraPostprocessorBuilder>
  type, extends(stellarPopulationSpectraPostprocessorBuilderClass) :: stellarPopulationSpectraPostprocessorBuilderLookup
     !% A stellar population spectra postprocessor builder which simply looks up postprocessors by name.
     private
     type(varying_string                           ), allocatable, dimension(:) :: names
     type(stellarPopulationSpectraPostprocessorList), allocatable, dimension(:) :: postprocessors
   contains
     final     ::          lookupDestructor
     procedure :: build => lookupBuild
  end type stellarPopulationSpectraPostprocessorBuilderLookup

  interface stellarPopulationSpectraPostprocessorBuilderLookup
     !% Constructors for the {\normalfont \ttfamily lookup} stellar population spectra postprocessor builder class.
     module procedure lookupConstructorParameters
     module procedure lookupConstructorInternal
  end interface stellarPopulationSpectraPostprocessorBuilderLookup
  
contains

  function lookupConstructorParameters(parameters) result(self)
    !% Constructor for the {\normalfont \ttfamily lookup} stellar population spectra postprocessor builder class which takes a
    !% parameter list as input.
    use Input_Parameters
    use Galacticus_Error
    implicit none
    type   (stellarPopulationSpectraPostprocessorBuilderLookup)                              :: self
    type   (inputParameters                                   ), intent(inout)               :: parameters
    type   (varying_string                                    ), allocatable  , dimension(:) :: names
    type   (stellarPopulationSpectraPostprocessorList         ), allocatable  , dimension(:) :: postprocessors
    integer                                                                                  :: countPostprocessors, countPostprocessorNames, &
         &                                                                                      i

    countPostprocessors    =max(parameters%copiesCount('stellarPopulationSpectraPostprocessorMethod',zeroIfNotPresent=.true.),1)
    countPostprocessorNames=max(parameters%      count('names'                                      ,zeroIfNotPresent=.true.),1)
    if (countPostprocessors /= countPostprocessorNames) call Galacticus_Error_Report('number of names must match number of postprocessors'//{introspection:location})
    allocate(names         (countPostprocessors))
    allocate(postprocessors(countPostprocessors))
    !# <inputParameter>
    !#   <name>names</name>
    !#   <cardinality>1</cardinality>
    !#   <defaultValue>[var_str('default')]</defaultValue>
    !#   <description>The names assigned to stellar spectra postprocessors.</description>
    !#   <source>parameters</source>
    !#   <type>string</type>
    !# </inputParameter>
    !# <objectBuilder class="stellarPopulationSpectraPostprocessor" name="postprocessors(i)%stellarPopulationSpectraPostprocessor_" source="parameters" copy="i=1,countPostprocessors"/>
    self=stellarPopulationSpectraPostprocessorBuilderLookup(names,postprocessors)
    !# <inputParametersValidate source="parameters"/>
    return
  end function lookupConstructorParameters
  
  function lookupConstructorInternal(names,postprocessors) result(self)
    !% Internal constructor for the {\normalfont \ttfamily lookup} stellar population spectra postprocessor builder.
    use Galacticus_Error
    implicit none
    type(stellarPopulationSpectraPostprocessorBuilderLookup)                              :: self
    type(varying_string                                    ), intent(in   ), dimension(:) :: names
    type(stellarPopulationSpectraPostprocessorList         ), intent(in   ), dimension(:) :: postprocessors
    !# <constructorAssign variables="names, postprocessors"/>

    if (size(names) /= size(postprocessors)) call Galacticus_Error_Report('number of names must match number of postprocessors'//{introspection:location})
    return
  end function lookupConstructorInternal

  subroutine lookupDestructor(self)
    !% Destructor for the {\normalfont \ttfamily lookup} stellar population spectra postprocessor builder.
    implicit none
    type   (stellarPopulationSpectraPostprocessorBuilderLookup), intent(inout) :: self
    integer                                                                    :: i

    do i=1,size(self%postprocessors)
       !# <objectDestructor name="self%postprocessors(i)%stellarPopulationSpectraPostprocessor_"/>
    end do
    return
  end subroutine lookupDestructor

  function lookupBuild(self,descriptor)
    !% Return a stellar population spectra postprocessor by lookup via name.
    use Galacticus_Error
    implicit none
    class  (stellarPopulationSpectraPostprocessorClass        ), pointer       :: lookupBuild
    class  (stellarPopulationSpectraPostprocessorBuilderLookup), intent(inout) :: self
    type   (varying_string                                    ), intent(in   ) :: descriptor
    integer                                                                    :: i

    lookupBuild => null()
    do i=1,size(self%names)
       if (self%names(i) == descriptor) lookupBuild => self%postprocessors(i)%stellarPopulationSpectraPostprocessor_
    end do
    if (.not.associated(lookupBuild)) call Galacticus_Error_Report('unable to located postprocessor "'//descriptor//'"'//{introspection:location})
    return
  end function lookupBuild