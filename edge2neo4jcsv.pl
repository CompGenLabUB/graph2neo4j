=head1 NAME

edge2neo4jcsv.pl
graphcompare - A command-line tool to compare graph files in DOT or tabular format.


=head1 SYNOPSIS

    edge2neo4jcsv.pl   -wholegraph FILE \
                    -biogrid FILE \
                    -string FILE \
                    -ppaxe FILE \
                    -maxlvl INT \
                    -prefix PREFIX


=head1 OPTIONS

=over 8

=item B<-h>, B<-help>

Shows this help.


=item B<-w>, B<-wholegraph> <FILE>

REQUIRED. Wholegraph file is a required graph file (DOT or graphviz)
created by filter_interactions_to_graph2.pl.


=item B<-b>, B<-biogrid> <FILE>

REQUIRED. Biogrid file (tabular) contains index (col1), interactions (col2 and
col3), type of interactions (col4) and pubmed ID (col 5)from the BioGrid database.\n


=item B<-s>, B<-string> <FILE>

REQUIRED. String file (tabular) contains protein interaction (col1 and col2),
type of interaction (col3), and othe data (col4, col5, col6) from the string protein
database. Only the interactions are necessary and must be human only.


=item B<-pp>, B<-ppaxe> <FILE>

REQUIRED. Ppaxe file (tabular) is created from the ppaxe program that sifts through
pubmed articles looking for possible interactions. The file contains the Pubmed ID
(col1), the interaction (col2, col3) and the score(col4).


=item B<-m>, B<-maxlvl> <INT>

REQUIRED. Is an integer that represents the max level graph created, that is not
the wholegraph, from the filter_interactions_to_graph2.pl program. This value can
be acquiered by looking at the folder where the graphs are stored and finding the
highest number after "graph_lvl+".


=item B<-pr>, B<-prefix> <PREFIX>

This contains the prefix of the graph level files. Example: "all_graph_lvl+".


=item B<-o>, B<-output> <FILE>

REQUIRED. Output CSV file.


=back

=head1 FUNCTIONS

=cut

use warnings;
use strict;
use Getopt::Long qw(:config no_ignore_case);
use Data::Dumper;
use Pod::Usage;


# -----------------------------------------------------------------------------
=over 2

=item get_cmdline_arguments()

Reads command line arguments

=back
=cut
sub get_cmdline_arguments {
  # INPUTS:
  #     whole_graph_file
  #     interaction_prefix
  #     max_lvl
  #     interaction_type_file
  #     interaction_source_file
  # OPTIONS:
  #     \%options
  my $options = shift;

  # Print synopsis if no arguments
  pod2usage(
      -verbose => 0,
      -output => \*STDOUT
    ) if (!@ARGV);

  GetOptions(
    $options,
    'help|?',
    'verbose|?',
    'wholegraph=s',
    'maxlvl=s',
    'prefix=s',
    'biogrid=s',
    'string=s',
    'ppaxe=s',
    'output=s'
  );

  # Print help if -help option
  if (defined $options->{'help'}){
    pod2usage(
      -verbose => 1,
      -output => \*STDOUT
    );
  }
  my @required = qw( wholegraph maxlvl prefix biogrid string ppaxe output );
  foreach my $req (@required) {
    defined $options->{$req}
      or die "\n# Command line option $req is missing!\n";
  }
  return;
}


# -----------------------------------------------------------------------------
=over 2

=item init_graph()

Initializes whole graph

=back
=cut

sub init_graph {
  my $wholegraph = shift;
  my $maxlvl     = shift;
  my $verbose    = shift;

  print STDERR "\tReading graph file: ", $wholegraph, " ..." if $verbose;
  my %data       = ();
  open (my $fh, "<", $wholegraph)
    or die "Can't read $wholegraph : $!\n";

  while (my $line = <$fh>) {
    chomp($line);
    next unless $line =~ m/->/;
    my ($gene1, $gene2) = split /->/, $line;
    $gene1 = gene_cleaner($gene1);
    $gene2 = gene_cleaner($gene2);
    my $interaction1 = $gene1 . "->" . $gene2;
    my $interaction2 = $gene2 . "->" . $gene1;
    my @inter = ($interaction1, $interaction2);
    foreach my $interaction (@inter){
      $data{$interaction} = () unless exists $data{$interaction};
      $data{$interaction}->{"level"} = $maxlvl + 1;
      $data{$interaction}->{"genetic_interaction"} = 0;
      $data{$interaction}->{"physical_interaction"} = 0;
      $data{$interaction}->{"unknown_interaction"} = 0;
      $data{$interaction}->{"biogrid"} = 0;
      $data{$interaction}->{"string"}  = 0;
      $data{$interaction}->{"ppaxe"}   = 0;
      $data{$interaction}->{"ppaxe_score"}      = "NA";
      $data{$interaction}->{"biogrid_pubmedid"} = "NA";
      $data{$interaction}->{"ppaxe_pubmedid"}   = "NA";
    }
  }
  print STDERR "done.\n" if $verbose;
  return \%data;
}


# -----------------------------------------------------------------------------
=over 2

=item gene_cleaner()

Removes non-valid characters from genes by taking only the symbol between quotes.

=back
=cut
sub gene_cleaner {
  my $dirty_gene = shift;
  if ($dirty_gene =~ m/"(.*?)"/) {
    return $1;
  } else {
    die "Cannot clean $dirty_gene\n";
  }
}


# -----------------------------------------------------------------------------
=over 2

=item add_orphan_interactions()

This function initializes interactions that are misteriously
not present in the whole graph.

=back
=cut
sub add_orphan_interactions {
  my $interaction = shift;
  my $data        = shift;
  my $maxlvl      = shift;
  my $filename    = shift;
  my $verbose     = shift;
  print STDERR "\n\t\t[WARNING] Orphan found in lvl file: $filename => $interaction ";
  my ($gene1, $gene2) = split /->/, $interaction;
  my $rev_interaction = $gene2 . "->" . $gene1;
  foreach my $int ($interaction, $rev_interaction) {
    $data->{$int} = () unless exists $data->{$int};
    $data->{$int}->{"level"} = $maxlvl + 1;
    $data->{$int}->{"genetic_interaction"} = 0;
    $data->{$int}->{"physical_interaction"} = 0;
    $data->{$int}->{"unknown_interaction"} = 0;
    $data->{$int}->{"biogrid"} = 0;
    $data->{$int}->{"string"}  = 0;
    $data->{$int}->{"ppaxe"}   = 0;
    $data->{$int}->{"ppaxe_score"}      = "NA";
    $data->{$int}->{"biogrid_pubmedid"} = "NA";
    $data->{$int}->{"ppaxe_pubmedid"}   = "NA";
  }
  return;
}

# -----------------------------------------------------------------------------
sub get_level {
  # INPUTS: filename, data structure
  my $prefix  = shift;
  my $maxlvl  = shift;
  my $data    = shift;
  my $verbose = shift;
  print STDERR "\tAdding interaction level... " if $verbose;

  for my $currlvl (0..$maxlvl) {
    my $filename = $prefix . $currlvl . ".dot";
    open (my $fh, "<", $filename)
      or die "Cannot open $filename!\n";
    while (my $line = <$fh>) {
      chomp($line);
      next unless $line =~ m/->/;
      my ($gene1, $gene2) = split /->/, $line;
      $gene1 = gene_cleaner($gene1);
      $gene2 = gene_cleaner($gene2);
      my $interaction = $gene1 . "->" . $gene2;

      if (not exists $data->{$interaction}) {
        add_orphan_interactions($interaction, $data, $maxlvl, $filename, $verbose);
      }
      $data->{$interaction}{"level"} = $currlvl if ($data->{$interaction}{"level"} == ($maxlvl + 1));
    }
  }
  print STDERR "done.\n" if $verbose;
  return;
}

# -----------------------------------------------------------------------------
sub biogrid_type_converter {
  my $biogrid_type = shift;
  my %interaction_dict = (
    "Affinity Capture-Luminescence" => 1,
    "Affinity Capture-MS"      => 1,
    "Affinity Capture-RNA"     => 1,
    "Affinity Capture-Western" => 1,
    "Biochemical Activity"     => 1,
    "Co-crystal Structure"     => 1,
    "Co-fractionation" => 1,
    "Co-localization" => 1,
    "Co-purification" => 1,
    "Far Western" => 1,
    "FRET" => 1,
    "PCA" => 1,
    "Protein-peptide" => 1,
    "Protein-RNA" => 1,
    "Proximity Label-MS" => 1,
    "Reconstituted Complex" => 1,
    "Two-hybrid" => 1,
    "Dosage Growth Defect" => 0,
    "Dosage Lethality" => 0,
    "Dosage Rescue" => 0,
    "Negative" => 0,
    "Phenotypic Enhancement" => 0,
    "Phenotypic Suppression" => 0,
    "Positive" => 0,
    "Synthetic Growth Defect" => 0,
    "Synthetic Haploinsufficiency" => 0,
    "Synthetic Lethality" => 0,
    "Synthetic Rescue" => 0,
    "## UNKOWN ##"                  => 4
);

  my $tp = exists ($interaction_dict{$biogrid_type})
                 ? $interaction_dict{$biogrid_type}
                 : $interaction_dict{"## UNKOWN ##"};
  #print "$biogrid_type -> $tp\n";
  return $tp;
}

###############################################################################
##                            Source and Type                                ##
###############################################################################

# Note: Check file type before finishing this!

# -----------------------------------------------------------------------------
sub biogrid_filter {
  # INPUTS: filename, datastructure, data for biogrid type and source
  my $biogrid_source = shift;
  my $data           = shift;
  my $verbose        = shift;

  print STDERR "\tFiltering biogrid... " if $verbose;

  my $filename = $biogrid_source;
  open (my $fh, "<", $filename)
    or die "Cannot open $filename!\n";

  my $header = <$fh>;
  while (my $line = <$fh>) {
    chomp($line);
    my @columns = split /\t/, $line;
    my $gene1 = $columns[7];
    my $gene2 = $columns[8];
    my $type = $columns[12];
    my $PMID = $columns[14];

    $type = biogrid_type_converter($type);

    my $interaction = $gene1 . "->" . $gene2;
    if (exists ($data->{$interaction})) {
      $data->{$interaction}{"biogrid"} = 1;
      $data->{$interaction}{"biogrid_pubmedid"} = $PMID;
    }
    if (exists ($data->{$interaction}) and $type = 0) {
      $data->{$interaction}{"genetic_interaction"}++;
    }
    if (exists ($data->{$interaction}) and $type = 1) {
      $data->{$interaction}{"physical_interaction"}++;
    }
    if (exists ($data->{$interaction}) and $type = 4) {
      $data->{$interaction}{"unkown_interaction"}++;
    }
  }
  print STDERR "done.\n" if $verbose;
  return;
}

# -----------------------------------------------------------------------------
sub string_filter {
  # INPUTS: filename, datastructure, data for string type and source
  my $string_source = shift;
  my $data          = shift;
  my $verbose       = shift;
  print STDERR "\tFiltering string... " if $verbose;

  my $filename = $string_source;
  open (my $fh, "<", $filename)
    or die "Cannot open $filename!\n";

  while (my $line = <$fh>) {
    chomp($line);
    my ($gene1, $gene2, $type, undef, undef, undef) = split /\t/, $line;
    my $interaction = $gene1 . "->" . $gene2;
    if (exists ($data->{$interaction})) {
      $data->{$interaction}{"physical_interaction"}++;
      $data->{$interaction}{"string"} = 1;
    }
  }
  print STDERR "done.\n" if $verbose;
  return;
}

# -----------------------------------------------------------------------------
sub ppaxe_filter {
  # INPUTS: filename, datastructure, data for ppaxe source and PMID
  my $ppaxe_source = shift;
  my $data         = shift;
  my $verbose      = shift;
  print STDERR "\tFiltering ppaxe... " if $verbose;

  my $filename = $ppaxe_source;
  open (my $fh, "<", $filename)
    or die "Cannot open $filename!\n";

  while (my $line = <$fh>) {
    chomp($line);
    my ($gene1, $gene2, $PMID, $score, undef) = split /\t/, $line;
    my $interaction = $gene1 . "->" . $gene2;
    if (exists ($data->{$interaction})) {
      $data->{$interaction}{"unknown_interaction"}++;
      $data->{$interaction}{"ppaxe_pubmedid"} = $PMID;
      $data->{$interaction}{"ppaxe"} = 1;
      $data->{$interaction}{"ppaxe_score"} = $score;
    }
  } # $line from $fh

  close($fh);
  print STDERR "done.\n" if $verbose;
  return;
}


# -----------------------------------------------------------------------------
=over 2

=item print_csv()

Explain format for CSV/NEO4j

=back
=cut
sub print_csv {
  my $data    = shift;
  my $output  = shift;
  my $verbose = shift;
  print STDERR "\tPrinting csv to $output ... " if $verbose;

  open (my $fh, ">", $output)
    or die "Cannot open $output!\n";
  foreach my $interaction (keys %$data) {
    my ($gene1, $gene2) = split /->/, $interaction;
    # Skip interactions not present in biogrid OR string OR ppaxe
    if (not $data->{$interaction}{biogrid} and not $data->{$interaction}{string} and not not $data->{$interaction}{ppaxe}) {
      next;
    }
    print $fh "$gene1", ",", "$gene2";
    foreach my $attribute (qw(level genetic_interaction physical_interaction unknown_interaction biogrid biogrid_pubmedid string ppaxe ppaxe_score ppaxe_pubmedid)) {
      print $fh ",$data->{$interaction}{$attribute}";
    }
    print $fh "\n";
  }
  close($fh);

  print STDERR "done.\n" if $verbose;
  return;
}

# -----------------------------------------------------------------------------
# MAIN PROGRAM
my %options;
get_cmdline_arguments(\%options);

# START REPORT
my $start_time   = time();
my $current_time = localtime();
print STDERR "\nPROGRAM STARTED\n",
            "\tProgram         edge2neo4jcsv.pl
            "\tVersion         v0.1.0\n",
            "\tStart time      $current_time\n\n" if $options{'verbose'};


my $data = init_graph($options{'wholegraph'}, $options{'maxlvl'}, $options{'verbose'});
get_level($options{'prefix'}, $options{'maxlvl'}, $data, $options{'verbose'});
biogrid_filter($options{'biogrid'}, $data, $options{'verbose'});
string_filter($options{'string'}, $data, $options{'verbose'});
ppaxe_filter($options{'ppaxe'}, $data, $options{'verbose'});
print_csv($data, $options{'output'}, $options{'verbose'});


# END REPORT
if ($options{'verbose'}) {
  my $end_time  = time();
  $current_time = localtime();
  my $sec = $end_time - $start_time;

  my $hours = ($sec/3600) % 24;
  my $minutes = ($sec/60) % 60;
  my $seconds = $sec % 60;
  print STDERR "\nPROGRAM FINISHED\n",
               "\tEnd time \t$current_time\n\n",
               "\tJob took ~ $hours hours, $minutes minutes and $seconds seconds\n\n";
}
