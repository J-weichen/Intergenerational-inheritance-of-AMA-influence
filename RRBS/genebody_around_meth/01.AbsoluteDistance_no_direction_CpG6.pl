#!/usr/bin/perl -w
use strict;

my %promoter;
open IN,"$ARGV[0]" or die $!; # mm9.longest_Intragenic.xls or hg19.Intragenic.xls 
while(<IN>){
	chomp;
	my @a = split /\s+/,$_;
	my $s5 = $a[1];
	my $s3 = $a[2];
#	if ($a[4] eq "-"){
#		$s3 = $a[1];
#		$s5 = $a[2];
#	}
	$promoter{$a[0]}{$a[3]}{$s5} = $s3;
}
close IN;

my %CpG;
open IN,"$ARGV[1]" or die $!; # ../8cell.SingleCmet
while(<IN>){
	chomp;
	next if (/#/);
	my @a = split /\s+/,$_;
	next if ($a[5] ne "CG");
	my $total = int($a[3])+int($a[4]);
	next if ($total<6);

	$CpG{$a[0]}{$a[1]} = int($a[3])/$total;
}
close IN;


foreach my $chr (sort keys %promoter){
	foreach my $gene (sort keys %{$promoter{$chr}}){
		foreach my $s5 (sort keys %{$promoter{$chr}{$gene}}){
			my $s3 = $promoter{$chr}{$gene}{$s5};
			my $length = abs($s5-$s3);
			my $bin = int($length/100);
			if ($s5 > $s3){
				print "$gene";
				for(my $pos=$s5+15000;$pos-99>=$s5;$pos-=100){
                                        my ($sum,$n) = (0,0);
                                        foreach my $i (($pos-99)..$pos){
                                                if (exists $CpG{$chr}{$i}){
                                                        $sum += $CpG{$chr}{$i};
                                                        $n ++;
                                                }
                                        }
                                        my $mean = ($n != 0) ? $sum/$n : "NA";
                                        print "\t$mean";
                                }
				my $pos = $s5;
                        	for(my $i=1;$i<=100;$i++){
                                	my ($sum,$n) = (0,0);
                                	foreach my $j (($pos-$bin)..$pos){
                                        	if (exists $CpG{$chr}{$j}){
                                                	$sum += $CpG{$chr}{$j};
                                                	$n ++;
                                        	}
					}
                                	my $mean = ($n != 0) ? $sum/$n : "NA";
                                	print "\t$mean";
					$pos -= $bin;
                        	}
				for(my $pos=$s3;$pos-99>=$s3-15000;$pos-=100){
                                        my ($sum,$n) = (0,0);
                                        foreach my $i (($pos-99)..$pos){
                                                if (exists $CpG{$chr}{$i}){
                                                        $sum += $CpG{$chr}{$i};
                                                        $n ++;
                                                }
                                        }
                                        my $mean = ($n != 0) ? $sum/$n : "NA";
                                        print "\t$mean";
                                }
                        	print "\n";
			}
			else{
				print "$gene";
				for(my $pos=$s5-15000;$pos+99<=$s5;$pos+=100){
                                        my ($sum,$n) = (0,0);
                                        foreach my $i ($pos..($pos+99)){
                                                if (exists $CpG{$chr}{$i}){
                                                        $sum += $CpG{$chr}{$i};
                                                        $n ++;
                                                }
                                        }
                                        my $mean = ($n != 0) ? $sum/$n : "NA";
                                        print "\t$mean";
                                }
				my $pos = $s5;
                                for(my $i=1;$i<=100;$i++){
                                        my ($sum,$n) = (0,0);
                                        foreach my $j ($pos..($pos+$bin)){
                                                if (exists $CpG{$chr}{$j}){
                                                        $sum += $CpG{$chr}{$j};
                                                        $n ++;
                                                }
                                        }
                                        my $mean = ($n != 0) ? $sum/$n : "NA";
                                        print "\t$mean";
					$pos += $bin;
                                }
				for(my $pos=$s3;$pos+99<=$s3+15000;$pos+=100){
                                        my ($sum,$n) = (0,0);
                                        foreach my $i ($pos..($pos+99)){
                                                if (exists $CpG{$chr}{$i}){
                                                        $sum += $CpG{$chr}{$i};
                                                        $n ++;
                                                }
                                        }
                                        my $mean = ($n != 0) ? $sum/$n : "NA";
                                        print "\t$mean";
                                }
				print "\n";
			}
		}
	}
}
