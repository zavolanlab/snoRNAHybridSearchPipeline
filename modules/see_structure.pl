#! /usr/bin/env perl
#
# Short description for see_structure.pl
#
# Author Rafal Gumienny <r.gumienny@unibas.ch>
# Version 0.1
# Copyright (C) 2015 Rafal Gumienny <r.gumienny@unibas.ch>
# Modified On 2015-09-02 08:57
# Created  2015-09-02 08:57
#
use strict;
use warnings;



sub parse_structure
{
  my $structure = shift;
    my $gap_left;
    my $gap_right;
    my $U_gap;
    my $U_near_gap_l;
    my $U_near_gap_r;
    my $b_i_gap;
    my $i_t_gap;
    my $t_i_gap;
    my $i_b_gap;
    my $stem_ass;
    my $stem_max_ass;
    my $stem_length;
    my $count=0;
  ###################
  $structure=~/^([^\|]+)\|/;
  # print "1. $structure\n";
  my $temp=$1;
  # print "2. $temp\n";
  while($temp=~/<\.</gi){
    $count++;
  }
  # print "3.gap_left_count $count\n";
  $gap_left=$count;
  $count=0;
   ####################
  if($structure=~/<\.\|\.</){
    $U_gap=3;
  }else{
    $U_gap=2;
  }
  # print "U_gap $U_gap, \n";
  # print $structure=~/<\.\|\.</;
  #################
  $structure=~/\|\.([^\&]+)/;
  $temp=$1;
  while($temp=~/<\.</gi){
    $count++;
  }
  $gap_right=$count;
  #################
  if($structure=~/<\.<\.{0,1}\|/){
    $U_near_gap_l =1 ;
  }
  else{
    $U_near_gap_l =0 ;
  }
  #################
  if($structure=~/\|\.<\.</){
    $U_near_gap_r =1 ;
  }
  else{
    $U_near_gap_r =0 ;
  }
  ##################
  $structure=~/&([^>]+>)/;
  $temp=$1;
  if(defined($1)){
    if($temp=~/(\(\.*\>)/){
      $b_i_gap = length($1)-2;
    }
    else{
      #    print $structure," ",$temp,"\n";exit(0);
      $b_i_gap="100";
    }
  }
  else{
    $b_i_gap="100";
  }
  ######################
  $structure=~/>(\.*)\(/g;
  $temp=$1;
  print "tmp: $temp\n";
  $i_t_gap = length($temp);
  #########################
  $structure=~/>\.*(\(.+\))\.*>/;
  $temp=$1;
  $stem_length=length($temp);
  if($stem_length==0){
    print $structure," ",$temp,"\n";
    exit(0);
  }
  my @temp_array;
  my $temp_pos;
  my $position = -2;

  while(($position)){
    $position=index($temp,"(",$position);
    $position++;
    if($position>0){
      push @{$temp_pos},$position;
    }

  }
  $position=-2;
  while(($position)){
    $position=index($temp,")",$position);
    $position++;
    if($position>0){
      push @{$temp_pos},$position;
    }
  }
  $stem_ass=0;
  $stem_max_ass=0;
  for(my $i=0; $i<$#$temp_pos/2; $i++){
    $stem_max_ass = $stem_max_ass > abs( $temp_pos->[$i+1]-$temp_pos->[$i] - ($temp_pos->[$#$temp_pos-$i] - $temp_pos->[$#$temp_pos-$i-1])) ? $stem_max_ass : abs( $temp_pos->[$i+1]-$temp_pos->[$i] - ($temp_pos->[$#$temp_pos-$i] - $temp_pos->[$#$temp_pos-$i-1]));
    $stem_ass += $temp_pos->[$i+1]-$temp_pos->[$i] - ($temp_pos->[$#$temp_pos-$i] - $temp_pos->[$#$temp_pos-$i-1]);
  }
  ############################
  $structure =~/(\)\.*\>)/;
  $temp=$1;
  $t_i_gap = length($temp)-2;
  ############################
  $structure =~/(\>[^>]*)$/;
  $temp=$1;
  if($temp=~/\>([^\)]*)\)/){
    $temp=$1;
    $i_b_gap=length($temp);
  }
  else{
    $i_b_gap="-1";
  }

  my @a = ($gap_left, $gap_right, $U_gap, $U_near_gap_l, $U_near_gap_r, $b_i_gap, $i_t_gap, $stem_length, $stem_max_ass, $stem_ass, $t_i_gap, $i_b_gap);
  my @b = ($t_i_gap, $U_gap, $i_b_gap, $i_t_gap, $gap_right, $stem_length, $stem_ass);

  return join(", ", @b);
}


my @structures = ("<<<|.<<<<<<<&..................(((((>>>>>>>((((((((((..........))).)))))........)).>>>...)))))",
				  "<<.|.<<<&>>>.(((..((((.(((((((((((....)))..))))..)))).))))...)))>>.........",
			      "<.<.|.<<<<&(>>>>..(((....(((((((.(((((.......)))))))))))))))..>.>.)......",
			      "<<<|.<<<<<<<&..................(((((>>>>>>>((((((((((..........))).)))))........)).>>>...)))))",
			      "<<.|.<<<&>>>.(((..((((.(((((((((((....)))..))))..)))).))))...)))>>.........",
			      "<|.<<<<<<<<<<&>>>>>>>>>>.((.(..(((..(((((((((((........))))))))))).))).)))>...........",
			      "<<<.<<<<.<.|.<<&.>>.(((..((((.(((((((((((....)))..))))..)))).))))...)))>.>>>>.>>>.",
			      "<<<.<<<|.<.<<<&...(..>>>.>.(((..(((..(((((((((((........))))))))))).)))..)))>>>.>>>..).",
			      "<<<|.<<<<&((..>>>>..((..((((((.(((.(((......))).))).))))))..))>>>.......)).",
			      "<<.<<|.<<<<&>>>>..(((.(((...........((((.(((((.((....)).)))))...))))))).)))..>>.>>......",
			      "<<.<|.<<<<<<&((((((...>>>>>>(((((.(((.(((......))).))).))))).>.>>)))))).......",
			      "<<<<<<|.<&(.........>..((((....((.((((.(((((.((....)).)))))...)))).)).)).))>>>>>>)....",
			      "<<.<<.|.<<<<<&>>>>>.(((((.((((.(((..............))).)))).))).)).>>.>>............",
			      "<<<<|.<<&(.((.....>>.((.(((.....(((((..(((((.......)))))..)))))..))).))>>>>)).)....",
			      "<<<<.<.|.<.<<<<<&..((.............>>>>>.>.(.(((((((..(((.....)))))))))).)>.>>>>..)).",
			      "<<.<<<<|.<<<<&(((.....>>>>..(((......(((((..(((((.......)))))..)))))......)))>>>>.>>.)))",
			      "<<<<<<<<|.<<<<.<&...(((>.>>>>..(((((.((.........))...)))))>>>>>>>>)))",
			      "<<<.<|.<<<&.....((..........>>>.(((((((......(((((((((((((.((....)).))))))))))))))))))))>.>>>.))...",
			      "<<<<<|.<.<<<<<<<<<<<.<<&(.>>.>>>>>>>>>>>.>..(((..((((.((....)).)))).)))>>>>>......)....",
			      "<<<<<<.<|.<<<<&.((>>>>.(((...........((((((((((.........))))))..))))))).>.>>>>>>.))......",
			      "<<<.|.<&(.>..((..(((...((((((((....))))))))...)))....))>>>........)....",
			      "<<<|.<<<<&.((.....>>>>.((.((.((.((((((((((.........))))))..)))).))...))))>>>))......");
# my $data = parse_structure($initial_struct);
# print $data;
open F, ">test_from_perl.tab";
# print F join("\t", ("gap_left", "gap_right", "U_gap", "U_near_gap_l", "U_near_gap_r", "b_i_gap", "i_t_gap", "stem_length", "stem_max_ass", "stem_ass", "t_i_gap", "i_b_gap\n"));
print F "(";
foreach (@structures){
	my $tmp = parse_structure($_);
	print F "('$_', [$tmp]),\n";
}
print F ")";
close F;
# print join(", ", @data);
# my $second_struct = "<<.|.<<<&>>>.(((..((((.(((((((((((....)))..))))..)))).))))...)))>>.........";
# my @data2 = parse_structure($second_struct);
