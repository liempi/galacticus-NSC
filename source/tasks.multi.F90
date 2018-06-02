!! Copyright 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018
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

  type, public :: multiTaskList
     class(taskClass    ), pointer :: task_
     type (multiTaskList), pointer :: next  => null()
  end type multiTaskList
  
  !# <task name="taskMulti">
  !#  <description>A task which performs multiple other tasks.</description>
  !# </task>
  type, extends(taskClass) :: taskMulti
     !% Implementation of a task which performs multiple other tasks.
     private
     type(multiTaskList), pointer :: tasks
   contains
     final     ::                       multiDestructor
     procedure :: perform            => multiPerform
     procedure :: requiresOutputFile => multiRequiresOutputFile
  end type taskMulti

  interface taskMulti
     !% Constructors for the {\normalfont \ttfamily multi} task.
     module procedure multiConstructorParameters
     module procedure multiConstructorInternal
  end interface taskMulti

contains

  function multiConstructorParameters(parameters) result(self)
    !% Constructor for the {\normalfont \ttfamily multi} task class which takes a parameter set as input.
    use Input_Parameters
    implicit none
    type   (taskMulti      )                :: self
    type   (inputParameters), intent(inout) :: parameters
    type   (multiTaskList  ), pointer       :: task_
    integer                                 :: i

    self %tasks => null()
    task_       => null()
    do i=1,parameters%copiesCount('taskMethod',zeroIfNotPresent=.true.)
       if (associated(task_)) then
          allocate(task_%next)
          task_ => task_%next
       else
          allocate(self%tasks)
          task_ => self%tasks
       end if
       task_%task_ => task(parameters,i)
    end do
    return
  end function multiConstructorParameters
  
  function multiConstructorInternal(tasks) result(self)
    !% Internal constructor for the {\normalfont \ttfamily multi} task class.
    implicit none
    type(taskMulti    )                        :: self
    type(multiTaskList), target, intent(in   ) :: tasks
    !# <constructorAssign variables="*tasks"/>

    return
  end function multiConstructorInternal
  
  elemental subroutine multiDestructor(self)
    !% Destructor for the {\normalfont \ttfamily multi} task class.
    implicit none
    type(taskMulti    ), intent(inout) :: self
    type(multiTaskList), pointer       :: task_, taskNext

    if (associated(self%tasks)) then
       task_ => self%tasks
       do while (associated(task_))
          taskNext => task_%next
          deallocate(task_%task_)
          deallocate(task_      )
          task_ => taskNext
       end do
    end if
    return
  end subroutine multiDestructor

  subroutine multiPerform(self)
    !% Perform all tasks.
    use Galacticus_Display
    implicit none
    class(taskMulti    ), intent(inout) :: self
    type (multiTaskList), pointer       :: task_

    call Galacticus_Display_Indent('Begin multiple tasks')
    task_ => self%tasks
    do while (associated(task_))
       call task_%task_%perform()
       task_ => task_%next
    end do
    call Galacticus_Display_Unindent('Done multiple tasks')
    return
  end subroutine multiPerform

  logical function multiRequiresOutputFile(self)
    !% Returns true if any sub-task requires that the output file be open.
    implicit none
    class(taskMulti    ), intent(inout) :: self
    type (multiTaskList), pointer       :: task_

    multiRequiresOutputFile =  .false.
    task_                   => self%tasks
    do while (associated(task_).and..not.multiRequiresOutputFile)
       multiRequiresOutputFile=task_%task_%requiresOutputFile()
       task_ => task_%next
    end do
    return
  end function multiRequiresOutputFile