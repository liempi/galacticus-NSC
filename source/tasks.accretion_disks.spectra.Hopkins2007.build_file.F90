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
  
  !# <task name="taskAGNSpectraHopkins2008BuildFile">
  !#  <description>A task which evolves galaxies within a set of merger tree forests.</description>
  !# </task>
  type, extends(taskClass) :: taskAGNSpectraHopkins2008BuildFile
     !% Implementation of a task which builds a file containing a tabulation of AGN spectra from the model of \cite{hopkins_observational_2007}.
     private
   contains
     procedure :: perform            => agnSpectraHopkins2008BuildFilePerform
     procedure :: requiresOutputFile => agnSpectraHopkins2008BuildFileRequiresOutputFile
  end type taskAGNSpectraHopkins2008BuildFile

  interface taskAGNSpectraHopkins2008BuildFile
     !% Constructors for the {\normalfont \ttfamily agnSpectraHopkins2008BuildFile} task.
     module procedure agnSpectraHopkins2008BuildFileParameters
  end interface taskAGNSpectraHopkins2008BuildFile

contains

  function agnSpectraHopkins2008BuildFileParameters(parameters) result(self)
    !% Constructor for the {\normalfont \ttfamily agnSpectraHopkins2008BuildFile} task class which takes a parameter set as input.
    use Input_Parameters
    implicit none
    type(taskAGNSpectraHopkins2008BuildFile)                :: self
    type(inputParameters                   ), intent(inout) :: parameters
    !GCC$ attributes unused :: parameters
    
    self=taskAGNSpectraHopkins2008BuildFile()
    return
  end function agnSpectraHopkins2008BuildFileParameters

  subroutine agnSpectraHopkins2008BuildFilePerform(self)
    !% Builds the tabulation.
    use Galacticus_Display    
    use Accretion_Disk_Spectra
    implicit none
    class(taskAGNSpectraHopkins2008BuildFile), intent(inout) :: self
    type (accretionDiskSpectraHopkins2007   )                :: accretionDiskSpectra_
    !GCC$ attributes unused :: self

    call Galacticus_Display_Indent  ('Begin task: hopkins2007 AGN spectra file build')
    accretionDiskSpectra_=accretionDiskSpectraHopkins2007()
    call Galacticus_Display_Unindent('Done task: hopkins2007 AGN spectra file build')
    return
  end subroutine agnSpectraHopkins2008BuildFilePerform

  logical function agnSpectraHopkins2008BuildFileRequiresOutputFile(self)
    !% Specifies that this task does not requires the main output file.
    implicit none
    class(taskAGNSpectraHopkins2008BuildFile), intent(inout) :: self    
    !GCC$ attributes unused :: self

    agnSpectraHopkins2008BuildFileRequiresOutputFile=.false.
    return
  end function agnSpectraHopkins2008BuildFileRequiresOutputFile