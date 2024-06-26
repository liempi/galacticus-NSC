!! Copyright 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018,
!!           2019, 2020, 2021, 2022, 2023, 2024
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
Contains a module which implements a class for computation and output of star formation histories for galaxies.
!!}

module Star_Formation_Histories
  !!{
  Implements a class for computation and output of star formation histories for galaxies.
  !!}
  use            :: Abundances_Structure      , only : abundances
  use            :: Galacticus_Nodes          , only : treeNode
  use            :: Galactic_Structure_Options, only : enumerationComponentTypeType
  use            :: Histories                 , only : history
  use, intrinsic :: ISO_C_Binding             , only : c_size_t
  use            :: Kind_Numbers              , only : kind_int8
  use            :: Locks                     , only : ompLock
  implicit none
  private

  !![
  <functionClass>
   <name>starFormationHistory</name>
   <descriptiveName>Star Formation Histories</descriptiveName>
   <description>
     Class providing recording and output of star formation histories.
   </description>
   <default>null</default>
   <method name="create" >
    <description>Create the star formation history object.</description>
    <type>void</type>
    <pass>yes</pass>
    <argument>type            (treeNode), intent(inout)           :: node</argument>
    <argument>type            (history ), intent(inout)           :: historyStarFormation</argument>
    <argument>double precision          , intent(in   )           :: timeBegin</argument>
    <argument>double precision          , intent(in   ), optional :: timeEnd</argument>
   </method>
   <method name="scales" >
    <description>Set ODE solver absolute scales for a star formation history object.</description>
    <type>void</type>
    <pass>yes</pass>
    <argument>type            (history   ), intent(inout) :: historyStarFormation</argument>
    <argument>double precision            , intent(in   ) :: massStellar</argument>
    <argument>type            (abundances), intent(in   ) :: abundancesStellar</argument>
   </method>
   <method name="rate" >
    <description>Record the rate of star formation in this history.</description>
    <type>void</type>
    <pass>yes</pass>
    <argument>type            (treeNode  ), intent(inout) :: node</argument>
    <argument>type            (history   ), intent(inout) :: historyStarFormation</argument>
    <argument>type            (abundances), intent(in   ) :: abundancesFuel</argument>
    <argument>double precision            , intent(in   ) :: rateStarFormation</argument>
   </method>
   <method name="metallicityBoundaries" >
    <description>Return a (zero-indexed) array of metallicity boundaries for this history.</description>
    <type>double precision, allocatable, dimension(:)</type>
    <pass>yes</pass>
   </method>
   <method name="times" >
    <description>Return an array of times for this history \emph{if} the tabulation in time is static per output.</description>
    <type>double precision, allocatable, dimension(:)</type>
    <pass>yes</pass>
    <argument>integer(c_size_t), intent(in   ) :: indexOutput</argument>
    <modules>Error</modules>
    <code>
      allocate(starFormationHistoryTimes(0))
      call Error_Report('times are not static'//{introspection:location})
    </code>
   </method>
   <method name="perOutputTabulationIsStatic" >
    <description>Return true if the tabulation (in time and metallicity) is static (independent of node) per output.</description>
    <type>logical</type>
    <pass>yes</pass>
    <code>
     !$GLC attributes unused :: self
     starFormationHistoryPerOutputTabulationIsStatic=.false.
    </code>
   </method>
   <method name="update">
     <description>Update the star formation history after an output time is reached.</description>
     <type>void</type>
     <pass>yes</pass>
     <argument>type   (treeNode), intent(inout), target :: node                </argument>
     <argument>integer(c_size_t), intent(in   )         :: indexOutput         </argument>
     <argument>type   (history ), intent(inout)         :: historyStarFormation</argument>
     <code>
       !$GLC attributes unused :: self, node, indexOutput, historyStarFormation
     </code>
   </method>
  </functionClass>
  !!]

end module Star_Formation_Histories
