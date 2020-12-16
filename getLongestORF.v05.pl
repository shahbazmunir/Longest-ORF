#!/usr/bin/perl
use strict;
use warnings;
use List::Util qw(sum);
use Term::ANSIColor;

=begin
	This script accepts the multifasta and gff file. It generates all possible reading frames and select the longest frame. It also converts the ORF coordinates w.r.t to gene coordinates. The gene coordinates was obtained from gff file.

	Section A: Preprocess the multifasta file and GFF file
		

	Section B: Find Largest ORF
		It process the multifasta file and generates all possible reading frames. It returns the ORF coordinates w.r.t to gene coordinates.

=cut

my ($Fastafile, $GFFfile) = @ARGV;
if(scalar @ARGV < 1){
	print "\nUsage:\n\t./getLongestORF.pl dna.fasta dna.fasta.gff\n\n";
	exit;
}

if(not defined $GFFfile){
	$GFFfile    = "$Fastafile.gff";
	generateGFF($Fastafile);

}

my $tmpfile    = "tmp.fa";
my $OutputFile = "$Fastafile.Annotation.tab";
my $NoORFfound = "$Fastafile.NoORF.txt";

preprocess($Fastafile, $GFFfile, $tmpfile);
EachGeneORF($tmpfile, $OutputFile);
system("rm $tmpfile");


#----------------Find Largest ORF------------------------------------
sub EachGeneORF{
	my $tmpfile = $_[0];
	my $OutputFile = $_[1];
	my $header = "FALSE";
	my $DNA = "";
	my $GeneAcc = "";
	my @arr;
	my @modifiedHeader;
	my @previousHeader;
	my $GeneStart;
	my $GeneEnd;
	my $LongestORF;

	#Write Longest ORF for each gene to a separate file
        open(my $OUTPUT, ">$OutputFile") or die "Could not open file '$OutputFile' $!";

	#Genes with no ORF	
	open(my $NOORF, ">$NoORFfound") or die  "Unable to open file '$NoORFfound' $!";

	print $OUTPUT "#GeneAccession\tFeature\tFrame\tFrom\tTo\tCDS_Length\tDNA_Seq\tCDS\n";
	open(DATA, "<$tmpfile");
	while(<DATA>){	#Read Multifasta file
		if($_ =~ m/^>/g){	#identify is it header or DNA seq?
			#Extract gene start and gene end
			@modifiedHeader = split("\t", $_);
			chomp @modifiedHeader;
			#print "$modifiedHeader[1]";

			#Extract Gene Accession number
			$_ =~ tr/\|/ /;
			@arr = split(" ", $_);
			$header = "TRUE";
		}else{
			$header = "FALSE";
		}
	
		if($header eq "FALSE"){	#Prepare DNA sequence
			chomp $_;
			$GeneAcc = $arr[0];
			$GeneAcc =~ s/^>//;
                        $GeneStart = $modifiedHeader[1];
                        $GeneEnd = $modifiedHeader[2];
			@previousHeader = @modifiedHeader;
			$DNA = $DNA.$_;
			$DNA = uc($DNA);
		}
		if($header eq "TRUE" && $DNA ne ""){
			#Find all six reading frames and return the longest Reading frame
			$LongestORF = LongestORF(ReadingFrame($DNA, $GeneAcc, "+", $GeneStart, $GeneEnd), ReadingFrame(ReverseComplement($DNA), $GeneAcc, "-", $GeneStart, $GeneEnd));

			# split the contig into 5'UTR and 3'UTR on the basis of CDS coordinates
			displayORF($LongestORF, $GeneStart, $GeneEnd, $OUTPUT, $NOORF, @previousHeader);
		        $DNA = "";
			#$GeneStart = "";
			#$GeneEnd = "";

		}#end of last if condition
	
	}#end of while loop

	#For the last sequence
	#Find all six reading frames and return the longest Reading frame
	$LongestORF = LongestORF(ReadingFrame($DNA, $GeneAcc, "+", $GeneStart, $GeneEnd), ReadingFrame(ReverseComplement($DNA), $GeneAcc, "-", $GeneStart, $GeneEnd));
	displayORF($LongestORF, $GeneStart, $GeneEnd, $OUTPUT, $NOORF, @previousHeader);

	close DATA;
	close $NOORF;
	close $OUTPUT;

}


sub displayORF{
	chomp $_[0];
	my @LORF = split("\t", $_[0]);
	
	# if no open reading frame was found
	if((scalar @LORF) != 0){# if open reading frame found then in which direction(forward/reverse)
		my $geneStart = $_[1];
		my $geneEnd = $_[2];
		my $OUTPUT_WRITE = $_[3];

		my $stringA;
		my $stringB;

		# if the open reading frame is in reverse direction
		if($LORF[2] =~ m/\+/g){
			$stringA = "5'UTR";
			$stringB = "3'UTR";	
		}else{
			$stringA = "3'UTR";
			$stringB = "5'UTR";
		}
	
		# split the contig into 5'UTR and 3'UTR on the basis of CDS coordinates
		print $OUTPUT_WRITE "$LORF[0]\t$stringA\t.\t$geneStart\t".($LORF[3]-1)."\t.\t.\t.\n";
		print $OUTPUT_WRITE "$_[0]\n";
		# if there is partial CDS/ORF, then add 0 0 as gene 3'UTR start and end coordinate
		if($LORF[4] == $geneEnd){
			print $OUTPUT_WRITE "$LORF[0]\t$stringB\t.\t0\t0\t.\t.\t.\n";
		}else{
			print $OUTPUT_WRITE "$LORF[0]\t$stringB\t.\t".($LORF[4]+1)."\t$geneEnd\t.\t.\t.\n";
		}
	}else{
		# if no open reading frame was found
		my $NOORF_WRITE = $_[4];
		my @previousheader = $_[5];
		print $NOORF_WRITE $previousheader[0]."\n";
		
	}

}
#______________________________*********************Subroutines**************______________________________________
#Generate Reading frames of a DNA sequences
sub ReadingFrame{
        my $DNA = uc($_[0]);
	my $GeneAcc = $_[1];
	my $Strand = $_[2];
 	my $GeneStart = $_[3];
	my $GeneEnd = $_[4];
        my $Triplet = "";
        my $ORF = "";
	my $check = "FALSE";
	my $AminoAcid = "";
	my $ORFStart;
	my $ORFEnd;
	my $seq = "";
	my $Output = "";

	#print "GeneAcc\tFeature\tFrame\tFrom\tTo\tLength\tDNA_Seq\tCDS\n";
        for(my $j=0;$j<3;$j++){	#Generate three reading frames either the sequenc is on forward or reverse strand
                #my $readingFrame = "";
		my $MCount = 0;
		$ORF = "";
		$seq = "";

                for(my $i=$j;$i<length($DNA);$i=$i+3){	#Generate Triplets
                        $Triplet = substr($DNA, $i ,3);
			$AminoAcid = Translate($Triplet);	#Find the approptiate amino acid
			if($AminoAcid eq 'M'){	#Start of the ORF
				$check = "TRUE";
				if(++$MCount == 1){
					$ORFStart = $i;
				}
			}
			if($check eq "TRUE"){
				$ORF .= $AminoAcid;
				$seq .= $Triplet;
			}
			if($AminoAcid eq '*' && $check eq "TRUE" && $MCount > 0){	#End of the ORF
				$check = "FALSE";
				$ORFEnd = $i + 2;
				$MCount = 0;
			}elsif($i+3 >= (length $DNA) && $check eq "TRUE" && $MCount > 0){#if no stop codon found
				$check = "FALSE";
				$ORFEnd = $i + length($Triplet)-1;
				$MCount = 0;
			}
			if($ORF ne "" && $check eq "FALSE"){	#Output this ORF
				#$Output is a multiline variable which stores all ORFs for current DNA sequence
                                $Output .= "$GeneAcc\tCDS\t$Strand".sum($j+1)."\t".sum($ORFStart+$GeneStart)."\t".sum($ORFEnd+$GeneStart)."\t".(length $ORF)."\t$seq\t$ORF\n"; #$seq\t$ORF
				#$Output .= "$GeneAcc\t$Strand".sum($j+1)."\t".sum($ORFStart+1+$GeneStart)."\t".sum($ORFEnd+1+3+$GeneStart)."\t".sum(($ORFEnd+1+3)-($ORFStart+1))."\t.\t.\n";
				$ORF = "";
				$seq = "";
			}

                }
        }#end of the outer loop

	return $Output;
}

#Longest Reading Frame
sub LongestORF{
	#
	my $ORFs = $_[0].$_[1];
	my $LongestORF;
	
	#print "$ORFs\n\n";
	$LongestORF = `echo \"$ORFs\" | sort -nrk6 | head -1`;
	#print $LongestORF;

	return $LongestORF;
}


sub ReverseComplement{
        my $DNA = $_[0];
        my @arr = split("", $DNA);
        my $CompDNA;
        for(my $k=0;$k<length($DNA);$k++){
                #$CompDNA = $CompDNA.@arr[$k];

                if($arr[$k] eq "A"){
                        $CompDNA =$CompDNA."T";
                }elsif($arr[$k] eq "T"){
                        $CompDNA =$CompDNA."A";
                }elsif($arr[$k] eq "G"){
                        $CompDNA =$CompDNA."C";
                }elsif($arr[$k] eq "C"){
                        $CompDNA =$CompDNA."G";
                }elsif($arr[$k] eq "N"){
                        $CompDNA =$CompDNA."N";
                }
        }#end of for loop
	#print  scalar reverse($CompDNA);
        return scalar reverse($CompDNA);
}


#Translate nucleotides into amino acids
sub Translate{
        my $Triplet = $_[0];
        my $aa;
        if($Triplet =~ m/GC[ATGCN]/g){  #Ala / A
                return "A";
        }elsif($Triplet =~ m/TG[CT]/g){ #Cys / C
                return "C";
        }elsif($Triplet =~ m/ATG/g){ #Met / M
                return "M";
        }elsif($Triplet =~ m/GA[TC]/g){ #Asp / D
               return "D";
        }elsif($Triplet =~ m/GA[AG]/g){ #Glu / E
                return "E";
        }elsif($Triplet =~ m/TT[TC]/g){ #Phe / F 
                return "F";
        }elsif($Triplet =~ m/GG[ATGCN]/g){ #Gly / G
                return "G";
        }elsif($Triplet =~ m/CA[TC]/g){ #His / H
                return "H";
        }elsif($Triplet =~ m/AT[ATC]/g){ #Ile / I
                return "I";
        }elsif($Triplet =~ m/AA[AG]/g){ #Lys / K
                return "K";
        }elsif($Triplet =~ m/CT[ATGCN]|TT[AG]/g){ #Leu / L 
                return "L";
        }elsif($Triplet =~ m/AA[TC]/g){ #Asn / N
                return "N";
        }elsif($Triplet =~ m/CC[ATGCN]/g){ #Pro / P
                return "P";
        }elsif($Triplet =~ m/CA[AG]/g){ #Gln / Q
                return "Q";
        }elsif($Triplet =~ m/CG[ATGCN]|AG[AG]/g){ # Arg / R
                return "R";
        }elsif($Triplet =~ m/TC[ATGCN]|AG[CT]/g){ #Ser / S
                return "S";
        }elsif($Triplet =~ m/AC[ATGCN]/g){ #Thr / T
                return "T";
        }elsif($Triplet =~ m/GT[AGCTN]/g){ #Val / V
                return "V";
        }elsif($Triplet =~ m/TGG/g){ #Trp / W
                return "W";
        }elsif($Triplet =~ m/TA[TC]/g){ #Tyr / Y
                return "Y";
        }elsif($Triplet =~ m/TA[GA]|TGA/g){ #*
                return "*";
        }

}


#-------------------------Preprocess the fasta and GFF file-------------

sub preprocess{
	my $Fastafile = $_[0];
	my $GFFfile   = $_[1];
	my $WriteFile = $_[2];
	
	my %GeneList;

	print "Preprocessing the fasta and GFF file...\n";

	# Process GFF file and prepare the GeneList
	open(DATA, "<$GFFfile");
	while(<DATA>){
		my $line = $_;
		chomp;
		if($_ =~ m/^#/g){
			next;
		}
		$line =~ tr/\t;/ /;
		my @geneInfo = split(" ", $line);
		$geneInfo[8] =~ s/ID=//;
		#print "$geneInfo[8]\n";
		# HashMap of GeneList
		$GeneList{$geneInfo[8]} = "$geneInfo[3]\t$geneInfo[4]";

	}

	# Write preprocessed data to a temporary file
	open(my $write, ">$WriteFile") or die "Could not open file '$WriteFile' $!";



	# Process the fasta file and concatenate the GeneStart/GeneEnd in the header
	my $GeneAbsent = "FALSE"; 
	open(DATA2, "<$Fastafile");
	while(<DATA2>){
		if($_ =~ m/^>/g){ # Header
			my $line2 = $_;
		        chomp;
		        #Extract Gene Accession number
		        $line2 =~ s/^>//;
			#print "$line2";
		        my @header = split(" ", $line2);
			#$header[0] =~ s/gnl\|UG\|//g;
			#print "$header[0]\n";
		        if(exists $GeneList{$header[0]}){
				my($start, $end) = split("\t", $GeneList{$header[0]});		
				#print $write "$_\t/GeneStart=$start /GeneEnd=$end\n";
				$_ =~ s/gnl\|UG\|//;
				print $write "$_\t$start\t$end\n";
				$GeneAbsent = "FALSE";
			}else{
				$GeneAbsent = "TRUE";
				print color("red"), "\nFollowing sequence was absent in GFF file\n", color("reset");
				print $_."\n";
				#print $write "$_\t/GeneStart=-1 /GeneEnd=-1\n";
				#print $write "$_\t-1\n";
			}
		}elsif($GeneAbsent eq "TRUE"){
			print $_;
		}else{# print sequences
			print $write $_;
		}

	}

	close $write;
	close DATA;
	close DATA2;
	print "Done\nGenerating ORF...\n";
}


#------------------Generate gff file--------------------------

sub generateGFF{
	my $fastaFile = $_[0];
	my $outFile = "$fastaFile.gff";

	open(WRITE, ">$outFile") or die "unbale to write data to '$outFile' $!";
	print WRITE "##gff-version 3\n";

	my @arr;
	my $DNA_len;
	my $DNA = "";
	my $header;

	open(DATA3, "<$fastaFile") or die "Unable to open file '$fastaFile' $!";
	while(<DATA3>){
		chomp $_;
		# split the heder into gff format
		if($_ =~ m/^>/g){
			$_ =~ s/^>//g;
			if($DNA ne ""){
				@arr = split(" ", $header);
				$DNA_len = length $DNA;
				print WRITE "$arr[0]\t.\tgene\t1\t$DNA_len\t.\t.\t.\t$header\n";
				$DNA = "";
			}
			$header = $_;	
			next;
		}
		$DNA.= $_;
	}
	@arr = split(" ", $header);
	$DNA_len = length $DNA;
	print WRITE "$arr[0]\t.\tgene\t1\t$DNA_len\t.\t.\t.\t$header\n";

	close DATA3;
	close WRITE;
}










