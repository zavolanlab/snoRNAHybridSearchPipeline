use strict;
use warnings;

#usage
# scriptname mapped_sequences genome_mappings

my %copies = ();
my %anno   = ();
open( F, "$ARGV[0]" );
while (<F>) {
  chomp;
  next if ( $_ =~ m/^id/ );
  my @F = split(/\t/);
  next if ( $F[$#F] == 1 );
  next if ( $F[ $#F - 1 ] eq 'bacterial' );
  next if ( $F[ $#F - 1 ] eq 'fungus' );
  next if ( $F[ $#F - 1 ] eq 'vector' );
  next if ( $F[ $#F - 1 ] eq 'viral' );

  # only take unique-mappers
  next if ( $F[11] eq 'NULL' );
  next if ( $F[11] != 1 );
  $copies{ $F[0] } = $F[2];
  $anno{ $F[0] }   = $F[ $#F - 1 ];
}
close(F);

open( G, "$ARGV[1]" );
while (<G>) {
  chomp;
  next if ( $_ =~ m/^id/ );
  my @F = split(/\t/);
  next if ( not defined $copies{ $F[1] } );
  print "seq$F[1]\t$F[2]\t$F[5]\t$F[6]\t$F[7]\t$copies{$F[1]}\t$anno{$F[1]}\n";
}
close(G);

