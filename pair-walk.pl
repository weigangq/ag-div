#!/usr/bin/env perl
use strict;
use warnings;
use Bio::AlignIO;
use Bio::SeqIO;
use Inline::Files;
use Algorithm::Numerical::Sample  qw /sample/;
use Data::Dumper;
use Bio::LocatableSeq;
use Getopt::Long;
#use Statistics::Basic qw(:all);

die "Usage: $0 -s 5 (with mutation step from 5 [default 1]) -c 100 (repeated runs, default 100) lk.txt [this is likelihood table]\n" unless @ARGV >= 1;

$|++; # stop bufferred printing
my %opts;
GetOptions(\%opts,
#	   "all|a=s",
#     "mut_site|m=i",
     "mut_rep|c=i",
     "mut_step|s=i",
);


# 1. Parse likelihood table
my %lk_table;
my %seen_pos;
while(<>) {
    chomp;
    next if /^allele/;
    my ($allele, $pos, $pnear, $pfar, $pmed) = split;
    $lk_table{$allele}{$pos} = { 'pnear' => $pnear, 'pfar' => $pfar, 'pmed' => $pmed};
    $seen_pos{$pos}++;
}
#print Dumper(\%lk_table); exit;

# 2. Parse reference alignment
my $aln = Bio::AlignIO->new(-fh=>\*ALN)->next_aln();
my $seq_len = $aln->length();

my @seq_aas;
my @ne_ids;
foreach ($aln->each_seq){
  my @aa = split //, $_->seq();
  push @seq_aas, \@aa;
  push @ne_ids, $_->id() if $_->id =~ /^\S$/;
}
@ne_ids = sort @ne_ids;

my %ref_aas;
for (my $i=1; $i<=$seq_len; $i++){
  my %aas;
  foreach (@seq_aas){
    my @aa = @{$_};
    $aas{$aa[$i-1]}++ unless (($i<=30 ||  $i >= 207) && $aa[$i-1] eq '-'); # two sites with N-terminus signals
  }
  my @a = keys %aas;
  $ref_aas{$i} = \@a unless scalar @a ==1
}
#print Dumper(\%ref_aas); exit;
#print scalar keys %ref_aas, "\n"; exit;

my $count = $opts{'mut_rep'} || 100;
my $mutRate;
my ($al_start, $al_end, $id_start, $id_end, %diff_pos, %lk_seed);
for (my $i = 0; $i <= $#ne_ids; $i++) {
    $al_start = $aln->get_seq_by_id($ne_ids[$i]);
    $id_start = $al_start->id();
    for (my $j = 0; $j <= $#ne_ids; $j++) {
	next if $i == $j;
	$al_end = $aln->get_seq_by_id($ne_ids[$j]);
	$id_end = $al_end->id();
	%diff_pos = &diff_pair($al_start, $al_end);
#print Dumper(\@diff_pos); exit;
	%lk_seed = %{$lk_table{$id_end}};
	my $total_diff = scalar keys %diff_pos;

	for ($mutRate = 1; $mutRate <= $total_diff; $mutRate += $opts{mut_step} || 1){
	    my @fit;
	    for (my $i=0; $i<$count; $i++){
		my $gamete = &mutateGamete($al_start);
		push @fit, sprintf "%.6f", &obtain_fitness($gamete, $al_end);
	    }
	    print join "\t", ($id_start, $id_end, $mutRate, $total_diff, @fit);
	    print "\n";
#	    printf "%s\t%s\t%d\t%.4f\t%.4f\n", $id_start, $id_end, $mutRate, mean(@fit), stddev(@fit);
	}
#printf "%s\t%s\t%d\t%.4f\n", $id_end, $id_end, 0, &obtain_fitness($al_end, $al_end);
    }
}
exit;

sub diff_pair {
    my ($s1, $s2) = @_;
    my %dpos;
#    for (my $pos = 1; $pos <= $seq_len; $pos++) {
    foreach my $pos (sort {$a <=> $b} keys %seen_pos ) {
	my $query_aa = $s1->subseq($pos, $pos);
	my $ref_aa = $s2->subseq($pos, $pos);
	$dpos{$pos}{$query_aa} = $ref_aa unless $query_aa eq $ref_aa;
    }
    return %dpos;
}

sub mutateGamete {
  my $seq = shift;
  my $query = Bio::LocatableSeq->new(-id=>'query', -seq=>$seq->seq()); 
  my @siteToMutate = sample(-set => [ keys %diff_pos ], -sample_size => $mutRate );
  for (my $i = 0; $i < $#siteToMutate; $i++) {
      my $site = $siteToMutate[$i];
      my $a = $seq->subseq($site, $site);
      my $aa = $diff_pos{$site}{$a};
      my $q_seq = $query->seq();
      substr $q_seq, $site-1, 1, $aa; 
      $query->seq($q_seq);
  }
  return $query;
}

sub obtain_fitness {
  my $query = shift;
  my $target = shift;
  my @bayes_numr = (0, 0, 0);
  foreach my $pos (keys %lk_seed) {
    my $prob = $lk_seed{$pos};
    my $query_aa = $query->subseq($pos, $pos);
    my $ref_aa = $target->subseq($pos, $pos);
    if ($query_aa eq $ref_aa) {
        $bayes_numr[0] += log($prob->{'pnear'});
        $bayes_numr[1] += log($prob->{'pfar'});
        $bayes_numr[2] += log($prob->{'pmed'});
    } else {
        $bayes_numr[0] += log(1-$prob->{'pnear'});
        $bayes_numr[1] += log(1-$prob->{'pfar'});
        $bayes_numr[2] += log(1-$prob->{'pmed'});
    }
  }
  return exp($bayes_numr[0])/(exp($bayes_numr[0]) + exp($bayes_numr[1]) + exp($bayes_numr[2]))
}


__ALN__
CLUSTAL W (1.81) multiple sequence alignment


G                      MKKNTLSAILMTLFLFISCNNSGKDGNASTNSADESVKGPNLAEISKKITESNAVVLAVK
A3                     ------------------------------NSADESVKGPNLTEISKKITESNAVVLAVK
E3                     ------------------------------NSADESVKGPNLTEISKKITDSNAVVLAVK
E                      MKKNTLSAILMTLFLFISCNNSGKDGNASANSADESVKGPNLTEISKKITESNAVVLAVK
N                      MKKNTLSAILMTLFLFISCNNSGKDGNASTNSADESVKGPNLTEISKKITESNAVVLAVK
D3                     ------------------------------NPADESVKGPNLTEISKKITDSNAIVMAVK
D                      MKKNTLSAILMTLFLFISCNNSGKDGNTSANSADESVKGPNLTEISKKITDSNAVLLAVK
F3                     ------------------------------NSADESVKGPNLTEISKKITDSNAVVLAVK
C3                     ------------------------------NSADESVKGPNLAEISKKITESNAVVLAVK
H3                     ------------------------------NSADESVKGPNLTEISKKITNSNAVVLAVK
M                      MKKNTLSAILMTLFLFISCNNSGKDGNTSANSADESVKGPNLTEISKKITESNAVVLAVK
K                      MKKNTLSAILMTLFLFISCNNSGKDGNTSANSADESVKGPNLTEISKKITESNAVVLAVK
U                      MKKNTLSAILMTLFLFISCNNSGKDGNASANSADESVKGPNLAEISKKITESNAVVLAVK
A                      MKKNTLSAILMTLFLFISCNNSGKDGNTSANSADESVKGPNLTEISKKITDSNAVLLAVK
B                      MKKNTLSAILMTLFLFISCNNSGKDGNTSANSADESVKGPNLTEISKKITDSNAVLLAVK
I                      MKKNTLSAILMTLFLFISCNNSGKDGNTSANSADESVKGPNLTEISKKITESNAVVLAVK
C                      MKKNTLSAILMTLFLFISCNNSGKDGNASANSADESVKGPNLTEISKKITESNAVVLAVK
T                      MKKNTLSAILMTLFLFISCNNSGKDGNASVNSADESVKGPNLTEISKKITESNAVVLAVK
F                      MKKNTLSAILMTLFLFISCNNSGKDGNTSANSADESVKGPNLTEISKKITESNAVVLAVK
I3                     ------------------------------NSADESVKGPNLTEISKKITESNAVVLAVK
L                      MKKNTLSAILMTLFLFISCNNSGKDGNASVNSADESVKGPNLVEISKKITDSNAVVIAVK
H                      MKKNTLSAILMTLFLFISCNNSGKDGNASANSADESVKGPNLTEISKKITESNAVVLAVK
J                      MKKNTLSAILMTLFLFISCNNSGKDGNTSANSADESVKGPNLTEISKKITESNAVVLAVK
                                                     *.**********.*******:***:::***


G                      EVAALLSSIDELA-KAIGKKIEQN-GLGADANHNTSLLAGAHEISTLIKQK---LDGLKN
A3                     EVETLLASIDELA-KAIGQKIGAN-GLTAQAAQNGSLLAGAYAISSLITEK---LNKFKD
E3                     EIETLLSSIDELATKAIGKKIQQNSGLGVEANRNESLLAGACVISTLITEK---LSKLSG
E                      EVETLLASIDELATKAIGKKIGNN-GLEANQSKNTSLLSGAYAISDLIAEK---LNVLKN
N                      EVAALLSSIDELA-KAIGKKINNN-GLDDVQNFNASLLAGAHTISKLVTEK---LSKLKN
D3                     EVETLLTSIDELA-KAIGKKINNN-GLDALQNFNASLLGGAHTISKLITEK---LSKLNG
D                      EVEVLLSSIDELAKKAIGKKIDQNNALGTLDNHNGSLLAGAYAISALITEK---LSSIKD
F3                     EVETLLASIDEIAAKAIGKKIDQNNALGTLDNHNGSLLAGAYAISALITKK---LSSIKD
C3                     EVETLISSIDEIAAKAIGKKIKNDGSLDNDANHNGSLLAGAYAISTLITQK---LGGLKN
H3                     EVETLLASIDELANKAIGKQITAN-GLNAQAGQNGTLLAGAYAISVLIAEK---LNILKN
M                      EVETLLASIDEVAKKAIGNLIAQN-GLNAGANQNGSLLAGAYVISTLIAEK---LDGLKN
K                      EIETLLASIDELATKAIGKKIQQNGGLAVEAGHNGTLLAGAYTISKLITQK---LDGLKN
U                      EVEALLASIDEIGSKAIGKRIQAN-GLQDLQGQNGSLLAGAYAISNLITQKINVLNGLKN
A                      EVEALLSSIDEIAAKAIGKKIHQNNGLDTENNHNGSLLAGAYAISTLIKQK---LDGLKN
B                      EVEALLSSIDELA-KAIGKKIKNDGSLDNEANRNESLLAGAYTISTLITQK---LSKLNG
I                      EVETLLTSIDELA-KAIGKKIKNDVSLDNEADHNGSLISGAYLISTLITKK---ISAIKD
C                      EVETLLASIDELA-KAIGKKIKNDVSLDNEADNNGSLIAGAYTISTLITQK---LSKLNG
T                      EVETLLASIDELANKAIGQKIDQNNGLSVDAGHNGPLLAGAYAISALITEK---LNGLTI
F                      EIETLLSSIDELATKAIGQKIDAN-GLGVQANQNGSLLAGAYAISTLITQK---LSAL-N
I3                     EIETLLSSIDELATKAIGQKIDAN-GLGVQADQNGSLLAGAYAISTLITQK---LSAL-N
L                      EVETLLVSIDELA-KAIGKKIEAGGTLGSDGAHNGSLLAGAYKIATEITAN---LSKLKA
H                      EVETLLASINQLA-KAIGKKIDQNGTLGDDGGQNGSLLAGAYAISTVIIEK---LSTLKN
J                      EIETLLASIDELATKAIGKKIDNNAGLGAEVGQNGSLLAGAYAISTVIIEK---LSTLKN
                       *: .*: **:::. ****: *  .  *      * .*:.**  *:  :  :   :. :  


G                      -EGLNKEIEAAKKCSAAFTKKLADSNADLGVAAGNATDDNAKRAILKTHGH-EDKGGKEL
A3                     LEGLKEKINDAKKCSEAFTNKLKGEQATLGQQN--ASDDDAKKAILKTHGT-PDKGAKEL
E3                     SEGLKEKIEAAKNCSEQFTKKLKDSHGELGVQN--VTDDNAKKAILKTNND-KTKGADEL
E                      -EELKEKIDTAKQCSTEFTNKLKSEHAVLGLDN--LTDDNAQRAILKKHAN-KDKGAAEL
N                      SEGLKEKIEDAKKCSDDFTKKLQSSHAQLGVAGGATTDEEAKKAILRTNAI-KDKGADEL
D3                     SEGLKEKIEDAKKCSDDFTKKLQSSHAVLGVAGGATTDEEAKKAILRTNSV-KDKGAEEL
D                      SGELKAEIEKAKKCSESFTKKLSDNQAELGIEN--ATDDNAKKAILKTHNA-KDKGAEEL
F3                     SGELKAEIEKAKKCSEAFTAKLKGEHTDLGKED--VTNAHAQEAILKTNGA-NDKGAKEL
C3                     SEELKEKIEEAKKCNKAFTEKLKSSHAELGKQD--AQDDDAKKAILRTHNT-KDKGAEEL
H3                     -DELKEKIEDAKKCTKAFTDKLKSSHGNLGQA---ADDDNAKKAILKTNQAGKDKGAEEL
M                      SEELKEKIEDAKKCNKAFTDKLKSSHAELGIANGAATDANAKAAILKTNGT-KDKGAQEL
K                      SEKLKEKIENAKKCSEDFTKKLEGEHAQLGIEN--VTDENAKKAILITDAA-KDKGAAEL
U                      SEELKEKINEAKGCSEKFTKKLSESHADIGIQA--ATDANAKDAILKTNPT-KTKGAEEL
A                      -EGLKEKIDAAKKCSETFTNKLKEKHTDLGKEG--VTDADAKEAILKTNGT-KTKGAEEL
B                      SEGLKEKIAAAKKCSEEFSTKLKDNHAQLGIQG--VTDENAKKAILKANAAGKDKGVEEL
I                      SGELKAEIEKAKKCSEEFTAKLKGEHTDLGKEG--VTDDNAKKAILKTNND-KTKGADEL
C                      SEGLKEKIAEAKKCSEEFTKKLKEKHTDLGKKD--ATDVHAKEAILKTNGT-KDKGAAEL
T                      SEELKEKIEKAKKCSTGFTNKLKSGHAELGPVGGNATDENAKQAILKTHGN-VTKGAKEL
F                      SEDLKEKVAKVKKCSEDFTNKLKNGNAQLGLAA--ATDDNAKAAILKTNGT-NDKGAKEL
I3                     SENLKEKVAKVKKCSEDFTNKLKNGNAQLGLAA--ATDADAKEAILKTNGT-KTKGAEEL
L                      SEDLKEKITKAKECSEKFTDKLKSENVALGKQD--ASDDDAKKAILKTHND-ITKGAKEL
H                      VEELKEKITKAKDCSEKFAGKLKNEHASLGKKD--ATDDDAKKAILKTHGN-TDKGAKEL
J                      VEELKEKITKAKDCSEKFTKKLKDSHAELGKKD--ASDDDAKKAILKTNQA-NDKGAKEL
                          *: ::  .* *.  *: **   :  :*       : .*: ***  .     **  **


G                      KELSEAVKSLLKAAQAALANSVQELTSPVVAETPKKP
A3                     KDLSESVAGLLKVAKEMLANSVKELTSPVVAET----
E3                     EKLFKAVEVLSKAAKEMLTNAVKELTSPVVAET----
E                      EKLFKAVENLSKAAQDTLKNAVKELTSPIVAESPKKP
N                      EKLFKSVESLAKAAQDALANSVNELTGPVVAETPKKP
D3                     GKLFESVESLSKAAKDKLNNSVKELTSPVVAET----
D                      VKLSESVAGLLKAAQAILANSVKELTSPVVAESPKKP
F3                     KDLSDSVESLVKAAKEMLTNSVKELTSPVVAES----
C3                     DKLFKAVENLSKAAKEMLSNSVKELTSPVVAET----
H3                     EKLFESVRALSKAAQDALTNAVKELTSPVVAES----
M                      EKLFESVKNLSKAAQETLNNSVKELTSPVVAENPKKP
K                      EKLFKAVENLAKAAKEMLANSVKELTSPIVAESPKNP
U                      DKLFKAVENLSKAAKEMLANSVKELTSPVVAESPKKP
A                      GKLFESVEVLSKAAKEMLANSVKELTSPVVAESPKKP
B                      EKLSGSLESLSKAAKEMLANSVKELTSPVVVESPKKP
I                      EKLFESVKNLSKAAKEMLTNSVKELTSPVVAESPKKP
C                      EKLFESVENLAKAAKEMLSNSVKELTSPVVAESPKNP
T                      KDLSESVEALAKAAQAMLTNSVKELTSPVVAETPKKP
F                      KDLSDSVESLVKAAQVMLTNSVKELTSPVVAESPKKP
I3                     GKLFESVEVLSKAAKEMLANSVKELTSPVVAES----
L                      KELSESVETLLKAAKEMLANSVKELTSPVVAESPKKP
H                      KDLSDSVESLVKAAKEMLTNSVKELTSPVVAESPKKP
J                      KELFEAVESLSKAAKEMLNKSVKELTSPIVAESPKNP
                        .*  ::  * *.*:  * ::*:***.*:*.*.    


