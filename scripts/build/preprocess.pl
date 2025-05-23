#!/usr/bin/env perl
use strict;
use warnings;
use Cwd;
use lib $ENV{'GALACTICUS_EXEC_PATH'}."/perl";
use Galacticus::Build::SourceTree;
use File::Changes;
use utf8;
use open ":std", ":encoding(UTF-8)";

# Preprocess a Galacticus Fortran source file.
# Andrew Benson (17-April-2015)

# Get arguments.
die("Usage: preprocess.pl <infile> <outfile>")
    unless ( scalar(@ARGV) == 2 );
my $inputFileName  = $ARGV[0];
my $outputFileName = $ARGV[1];

# Parse the file to build a tree.
my $tree = &Galacticus::Build::SourceTree::ParseFile($inputFileName);

# Process the tree.
&Galacticus::Build::SourceTree::ProcessTree($tree);

# Analyze the preprocessed tree.
&Galacticus::Build::SourceTree::AnalyzeTree($tree)
    if ( exists($ENV{'GALACTICUS_PREPROCESSOR_ANALYZE'}) && $ENV{'GALACTICUS_PREPROCESSOR_ANALYZE'} eq "yes" );

# Get the serialized source code.
(my $codeSerialized, my $mappings) = &Galacticus::Build::SourceTree::Serialize($tree, stripMappings => 1);

# Serialize back to source code.
open(my $outputFile,">:raw",$outputFileName.".tmp");
print $outputFile $codeSerialized;
close($outputFile);
&File::Changes::Update($outputFileName,$outputFileName.".tmp", proveUpdate => "yes");

# Output line number mappings.
open(my $lineMapFile,">",$outputFileName.".lmap");
print $lineMapFile $mappings;
close($lineMapFile);

exit;
