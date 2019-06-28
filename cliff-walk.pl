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

$|++; # stop bufferred printing
my %opts;
GetOptions(\%opts,
     "all|a=s", # allele
     "mut_num|m=i",
     "mut_rep|c=i",
     "mut_step|s=i",
     "help|h",
);

die "Usage: $0 -a J (walk from allele J) -s 1 (with mutation step from 1) -m (num mutations) lk.txt [this is likelihood table]\n" if $opts{'help'};


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
foreach ($aln->each_seq){
  my @aa = split //, $_->seq();
  push @seq_aas, \@aa;
}

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
#print Dumper([keys %seen_pos]); exit;
#print scalar keys %ref_aas, "\n"; exit;

my $id_seed = $opts{all};
my $seq_seed = $aln->get_seq_by_id($id_seed);
my %lk_seed = %{$lk_table{$id_seed}};
#print  &obtain_fitness($seq_seed), "\n";

# 3. mutate
my $mutRate;
#my $mutRate = $opts{'mut_site'} || 1;
my $count = $opts{'mut_rep'} || 100;
if ($opts{mut_num}) {
    $mutRate = $opts{mut_num};
    for (my $i=0; $i<$count; $i++){
	my $gamete = &mutateGamete($seq_seed);
	printf "%s\t%d\t%.4f\t%s\n", $id_seed, $mutRate, &obtain_fitness($gamete), $gamete->id();
    }
} else {
    for ($mutRate=1; $mutRate<= scalar keys %seen_pos; $mutRate += $opts{mut_step} || 5){
	my @fit;
	for (my $i=0; $i<$count; $i++){
	    my $gamete = &mutateGamete($seq_seed);
	    push @fit, sprintf "%.4f", &obtain_fitness($gamete);
	}
	print join "\t", ($id_seed, $mutRate, @fit);
	print "\n";
#	printf "%s\t%d\t%.4f\t%.4f\n", $id_seed, $mutRate, mean(@fit), stddev(@fit);
    }
}
exit;


sub mutateGamete {
  my $seq = shift;
  my $query = Bio::LocatableSeq->new(-id=>'query', -seq=>$seq->seq()); 
  my @siteToMutate = sample(-set => [keys %seen_pos], -sample_size => $mutRate );
  my $qid = join "_", sort {$a <=> $b} @siteToMutate;
  $query->id("q|" . $qid);
  for (my $i = 0; $i < $#siteToMutate; $i++) {
      my $site = $siteToMutate[$i];
      my @aaa = @{$ref_aas{$site}};
      my $a = $seq->subseq($site, $site);
      my @aa_new;
      foreach (@aaa){
        push @aa_new, $_ unless $_ eq $a
      }
      my @sample = sample(-set=>\@aa_new);
      my $aa = shift @sample;
      my $q_seq = $query->seq();
      substr $q_seq, $site-1, 1, $aa; 
      $query->seq($q_seq);
  }
  return $query;
}

sub obtain_fitness {
  my $query = shift;
  my @bayes_numr = (0, 0, 0);
  foreach my $pos (keys %lk_seed) {
    my $prob = $lk_seed{$pos};
    my $query_aa = $query->subseq($pos, $pos);
    my $ref_aa = $seq_seed->subseq($pos, $pos);
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


