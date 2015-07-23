#!/import/bc2/soft/bin/perl5/perl -w
use lib "/import/bc2/home/zavolan/clipz/clipz/development.build/perl";
#!/usr/bin/perl -w

use strict;
use Getopt::Long;
$| = 1;

my $loci_file = '';
my $str_dir = '';
my $windowleft = 50000;
my $windowright = 50000;
my ($currentchr, $currentchrseq, $ext);
my %data;

$ext = "fa";

my $results = GetOptions ( 'loci:s'   => \$loci_file,
			   'dir:s'    => \$str_dir,
			   'windowleft:s' => \$windowleft,
			   'windowright:s' => \$windowright,
			   'chrfileext:s' => \$ext
			   );

unless ( $loci_file ne '' && $str_dir ne '' && $windowleft ne '' && $windowright ne '' ) {
    usage();
    exit;
}
open ( LOCI, "<$loci_file" ) || die "Cannot open $loci_file $!";
while (<LOCI>) {
    chomp;

    if($_ =~ /^(\S+)\s+\((\d+)\s+nc\)\s+(\d+)\D\D(\d+)\s+(\S+)\s+(\d+)\D\D(\d+)/) {
	my ($clone_id, $clone_len, $clone_begin, $clone_end, $chr_id, $chr_start, $chr_end) =
	    ($1, $2, $3, $4, $5, $6, $7);
	my ($begin, $end) = ($chr_start, $chr_end);
	my $strand;
	if($clone_begin > $clone_end) {
	    $strand = "-";
	    $chr_start -= $windowright;
	    $chr_end += $windowleft;
	}
	else {
	    $strand = "+";
	    $chr_start -= $windowleft;
	    $chr_end += $windowright;

	}

	$chr_start = $chr_start < 1 ? 1 : $chr_start;
 
	my $chr_len = $chr_end - $chr_start + 1;
	my $rec = "$clone_id|$clone_len|$chr_id|$strand|$begin|$end|$windowleft|$windowright $chr_start $chr_len $strand";
	push @{$data{$chr_id}}, $rec;
    }
    elsif($_ =~ /^(\S+)\s+(\S+)\s+([\+\-])\s+(\d+)\s+(\d+)/) {
	my ($clone_id, $chr_id, $strand, $chr_start, $chr_end) = ($1, $2, $3, $4, $5);
	my ($begin, $end) = ($chr_start, $chr_end);
	if($strand =~ /\-/) {
	    $chr_start -= $windowright;
	    $chr_end += $windowleft;
	}
	else {
	    $chr_start -= $windowleft;
	    $chr_end += $windowright;

	}

	$chr_start = $chr_start < 1 ? 1 : $chr_start;
	my $chr_len = $chr_end - $chr_start + 1;
	my $rec = "$clone_id|$chr_id|$strand|$begin|$end|$windowleft|$windowright $chr_start $chr_len $strand";
	push @{$data{$chr_id}}, $rec;
    }
}
#now process one chromosome at a time
foreach $currentchr (keys(%data)) {
    $currentchrseq = "";
    my @recs = @{$data{$currentchr}};
    for(my $i = 0; $i < @recs; $i++) {
	my ($name, $chr_start, $chr_len, $strand) = split(/\s+/, $recs[$i]);
	my $chr_seq = extract_sequence( $str_dir, $currentchr, $chr_start, $chr_len );
	print ">$name\n";
	if($strand =~ /\+/) {
	    print "$chr_seq\n";
	}
	else {
	    my $seq = $chr_seq;
	    $seq =~ tr/[ACGTNRYWSMKBDHVacgtnrywsmkbdhv]/[TGCANYRWSKMVDHBtgcanyrwskmvdhb]/;                        
	    my @s = split(//, $seq);
	    @s = reverse(@s);
	    $seq = join("", @s);
	    print "$seq\n";
	}
    }
}
close(LOCI);

#####
# End MAIN
#####

####################
# Subs
####################

sub usage {
    print "Usage: $0 <-l|--loci loci_file> <-d|--dir dir> <-w|--window window_size>\n";
}

sub extract_sequence {
    my $db_location = shift;
    my $chromosome = shift;
    my $start = shift;
    my $length = shift;
    my $chr_file = $db_location . $chromosome;
    if($ext) {
	$chr_file .=".$ext";
    }
    print STDERR "Extracting $chr_file\n";
    if(!($currentchrseq)) {
	if(open (SEQ, "<$chr_file")) {
	    while(<SEQ>) {
		if($_ !~ /^\>/) {
		    $_ =~ s/\s//g;
		    $currentchrseq .= $_;
		}
	    }
	    close(SEQ);
	}
	elsif(open (SEQ, "<$chr_file.fa")) {
	    while(<SEQ>) {
		if($_ !~ /^\>/) {
		    $_ =~ s/\s//g;
		    $currentchrseq .= $_;
		}
	    }
	    close(SEQ);
	}
	else {
	    die "Cannot open $chr_file: $!";
	}
    }
    my $l = length($currentchrseq);
    if($length > $l-$start) {
	$length = $l-$start;
    } 
    my $seq = substr($currentchrseq, $start-1, $length);

    return $seq;
}
