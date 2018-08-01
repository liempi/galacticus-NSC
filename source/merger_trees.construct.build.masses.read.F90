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

  !% Implementation of a merger tree masses class which reads masses from a file.

  !# <mergerTreeBuildMasses name="mergerTreeBuildMassesRead" defaultThreadPrivate="yes" abstract="yes">
  !#  <description>A merger tree masses class which samples masses from a distribution.</description>
  !# </mergerTreeBuildMasses>
  type, abstract, extends(mergerTreeBuildMassesClass) :: mergerTreeBuildMassesRead
     !% Implementation of a merger tree masses class which reads masses from a file.
     private
     type            (varying_string) :: fileName
     double precision                 :: massIntervalFractional
   contains
     !@ <objectMethods>
     !@   <object>mergerTreeBuildMassesRead</object>
     !@   <objectMethod>
     !@     <method>read</method>
     !@     <arguments>\doubleone\ mass\argout, \doubleone\ weight\argout</arguments>
     !@     <type>\void</type>
     !@     <description>Read the halo masses, and, optionally, weights, from file..</description>
     !@   </objectMethod>
     !@ </objectMethods>
     procedure                             :: construct => readConstruct
     procedure(readReadTemplate), deferred :: read
  end type mergerTreeBuildMassesRead

  abstract interface
     subroutine readReadTemplate(self,mass,weight)
       import mergerTreeBuildMassesRead
       class           (mergerTreeBuildMassesRead), intent(inout)                            :: self
       double precision                           , intent(  out), allocatable, dimension(:) :: mass, weight
     end subroutine readReadTemplate
  end interface

contains

  subroutine readConstruct(self,time,mass,massMinimum,massMaximum,weight)
    !% Construct a set of merger tree masses by reading from a file.
    use Sort
    use Memory_Management
    implicit none
    class           (mergerTreeBuildMassesRead), intent(inout)                            :: self
    double precision                           , intent(in   )                            :: time
    double precision                           , intent(  out), allocatable, dimension(:) :: mass        , weight     , &
         &                                                                                   massMinimum , massMaximum
    integer                                                                               :: i
    !GCC$ attributes unused :: time
    
    call self%read(mass,weight)
    if (allocated(weight)) then
       call Sort_Do(mass,weight)
    else
       call Sort_Do(mass       )
       call allocateArray(massMinimum,shape(mass))
       call allocateArray(massMaximum,shape(mass))
       do i=1,size(mass)
          massMinimum(i)=+mass(i)/sqrt(self%massIntervalFractional)
          massMaximum(i)=+mass(i)*sqrt(self%massIntervalFractional)
          if (i > 1 .and. massMinimum(i) < massMaximum(i-1)) then
             massMinimum(i  )=sqrt(           &
                  &                +mass(i-1) &
                  &                *mass(i  ) &
                  &               )
             massMaximum(i-1)=sqrt(           &
                  &                +mass(i-1) &
                  &                *mass(i  ) &
                  &               )
          end if
       end do
    end if
    return
  end subroutine readConstruct