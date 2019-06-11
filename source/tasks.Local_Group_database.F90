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

  !# <task name="taskLocalGroupDatabase">
  !#  <description>A task which updates the Local Group database.</description>
  !# </task>
  type, extends(taskClass) :: taskLocalGroupDatabase
     !% Implementation of a task which updates the Local Group database.
     private
   contains
     procedure :: perform            => localGroupDatabasePerform
     procedure :: requiresOutputFile => localGroupDatabaseRequiresOutputFile
  end type taskLocalGroupDatabase

  interface taskLocalGroupDatabase
     !% Constructors for the {\normalfont \ttfamily localGroupDatabase} task.
     module procedure localGroupDatabaseParameters
  end interface taskLocalGroupDatabase

contains

  function localGroupDatabaseParameters(parameters) result(self)
    !% Constructor for the {\normalfont \ttfamily localGroupDatabase} task class which takes a parameter set as input.
    use Input_Parameters
    implicit none
    type(taskLocalGroupDatabase)                :: self
    type(inputParameters       ), intent(inout) :: parameters
    !GCC$ attributes unused :: parameters
    
    self=taskLocalGroupDatabase()
    return
  end function localGroupDatabaseParameters

  subroutine localGroupDatabasePerform(self,status)
    !% Update the database.
    use Galacticus_Display      , only : Galacticus_Display_Indent, Galacticus_Display_Unindent
    use Galacticus_Error        , only : errorStatusSuccess
    use Interface_Local_Group_DB, only : localGroupDB
    implicit none
    class  (taskLocalGroupDatabase), intent(inout)           :: self
    integer                        , intent(  out), optional :: status
    type   (localGroupDB          )                          :: database
    !GCC$ attributes unused :: self

    call Galacticus_Display_Indent  ('Begin task: localGroupDatabase')
    database=localGroupDB()
    call database%update()
    call Galacticus_Display_Unindent('Done task: localGroupDatabase' )
    if (present(status)) status=errorStatusSuccess
    return
  end subroutine localGroupDatabasePerform

  logical function localGroupDatabaseRequiresOutputFile(self)
    !% Specifies that this task does not requires the main output file.
    implicit none
    class(taskLocalGroupDatabase), intent(inout) :: self    
    !GCC$ attributes unused :: self

    localGroupDatabaseRequiresOutputFile=.false.
    return
  end function localGroupDatabaseRequiresOutputFile