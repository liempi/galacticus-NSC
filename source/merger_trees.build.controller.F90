!! Copyright 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018,
!!           2019, 2020, 2021
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
Contains a module which provides controller objects for building merger trees.
!!}

module Merger_Tree_Build_Controllers
  !!{
  Provides controller objects for building merger trees.
  !!}
  use :: Galacticus_Nodes   , only : treeNode
  use :: Merger_Tree_Walkers, only : mergerTreeWalkerClass
  private

  !![
  <functionClass>
   <name>mergerTreeBuildController</name>
   <descriptiveName>Merger Tree Build Controllers</descriptiveName>
   <description>Class providing merger tree build controllers.</description>
   <default>uncontrolled</default>
   <method name="control" >
    <description>Control the behavior of a tree build.</description>
    <type>logical</type>
    <pass>yes</pass>
    <argument>type (treeNode             ), intent(inout), pointer :: node       </argument>
    <argument>class(mergerTreeWalkerClass), intent(inout)          :: treeWalker_</argument>
   </method>
  </functionClass>
  !!]

end module Merger_Tree_Build_Controllers