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

  !# <task name="taskBuildToolFSPS">
  !#  <description>A task which builds the FSPS tool.</description>
  !# </task>
  type, extends(taskClass) :: taskBuildToolFSPS
     !% Implementation of a task which builds the FSPS tool.
     private
   contains
     procedure :: perform            => buildToolFSPSPerform
     procedure :: requiresOutputFile => buildToolFSPSRequiresOutputFile
  end type taskBuildToolFSPS

  interface taskBuildToolFSPS
     !% Constructors for the {\normalfont \ttfamily buildToolFSPS} task.
     module procedure buildToolFSPSParameters
  end interface taskBuildToolFSPS

contains

  function buildToolFSPSParameters(parameters) result(self)
    !% Constructor for the {\normalfont \ttfamily buildToolFSPS} task class which takes a parameter set as input.
    use Input_Parameters
    implicit none
    type(taskBuildToolFSPS)                :: self
    type(inputParameters  ), intent(inout) :: parameters
    !GCC$ attributes unused :: parameters
    
    self=taskBuildToolFSPS()
    return
  end function buildToolFSPSParameters

  subroutine buildToolFSPSPerform(self,status)
    !% Builds the tabulation.
    use Galacticus_Error  , only : errorStatusSuccess
    use Galacticus_Display
    use Interfaces_FSPS
    implicit none
    class  (taskBuildToolFSPS), intent(inout)           :: self
    integer                   , intent(  out), optional :: status
    type   (varying_string   )                          :: fspsPath, fspsVersion
    !GCC$ attributes unused :: self

    call Galacticus_Display_Indent  ('Begin task: FSPS tool build')
    call Interface_FSPS_Initialize(fspsPath,fspsVersion,static=.true.)
    call Galacticus_DisplaY_Message('FSPS version '//fspsVersion//' successfully built in: '//fspsPath)
    if (present(status)) status=errorStatusSuccess
    call Galacticus_Display_Unindent('Done task: FSPS tool build')
    return
  end subroutine buildToolFSPSPerform

  logical function buildToolFSPSRequiresOutputFile(self)
    !% Specifies that this task does not requires the main output file.
    implicit none
    class(taskBuildToolFSPS), intent(inout) :: self    
    !GCC$ attributes unused :: self

    buildToolFSPSRequiresOutputFile=.false.
    return
  end function buildToolFSPSRequiresOutputFile