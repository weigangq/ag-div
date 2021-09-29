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
use Bio::SimpleAlign;

$|++; # stop bufferred printing
my %opts;
GetOptions(\%opts,
#     "num-rep|r=i",
     "seed-file|s=s", # at least 10
     "mut|m=i",
     "num-elite|e=i",
     "num-gen|n=i",
     "help|h"
);

if ($opts{'help'}) {
  print "Usage: $0 [--seed-file <FASTA> --mut (default 10) --num-elite (default 10) --num-gen (default 100)]\n";
  exit;
}

# 1. Parse reference alignment
my $aln = Bio::AlignIO->new(-fh=>\*ALN)->next_aln();
my $seq_len = $aln->length();

my (@seq_aas, @ref_seqs);
foreach ($aln->each_seq){
  push @ref_seqs, $_;
  my @aa = split //, $_->seq();
  push @seq_aas, \@aa;
}
my $nRef = scalar @ref_seqs;

my @site_aa;
my %varSite;
for (my $i=1; $i<=$seq_len; $i++){
  my @aa;
  my %aas;
  foreach (@seq_aas){
    my @seq = @{$_};
    push @aa, $seq[$i-1];
    $aas{$seq[$i-1]}++;
  }
  push @site_aa, \@aa;
  my @a = keys %aas;
  $varSite{$i} = \@a unless scalar @a ==1
}

my $numElites = $opts{'num-elite'} || 10; # select top 10 gametes
my @seed;
if ($opts{'seed-file'}) {
  my $seed_fas = Bio::SeqIO->new(-file=>$opts{'seed-file'}, -format=>'fasta');
  while (my $seq = $seed_fas->next_seq()){
    my $s = $seq->seq();
    push @seed, {seq=>$s, fit=>&dis2ref($s)}
  }
} else {
#  my $nRdm = $opts{'num-rep'} || 10;
#  my @seeds;
  for (my $i=0; $i<$numElites; $i++){
    my $s = &random_seq();
    push @seed, {seq=>$s, fit=>&dis2ref($s)}
  }
  @seed = sort {$a->{fit} cmp $b->{fit}} @seed;
#  @seed = @seeds[0..$numElites-1];
}
#print Dumper(\@seed); exit;

my $mutRate = $opts{'mut'} || 10; # introduce 10 mutations per sequence
my $gameteSize = 10; # 10 offspring per parent
my $n = $opts{'num-gen'} || 100; # num of generations

for (my $j=0; $j<$n; $j++){
  my @gametePool = @seed;
  foreach my $p (@seed){
    my $s = $p->{seq};

    for (my $i=0; $i<$gameteSize; $i++){
      my $mut = &mutate($s);
      push @gametePool, {seq=>$mut, fit=>&dis2ref($mut)};
    }
  }
  my @pop = sort {$a->{fit} cmp $b->{fit}} @gametePool;
  if ($pop[$numElites-1]->{fit} eq $seed[-1]->{fit}){
    print STDERR "Generation ", ($j+1), "\n"
  } else {
    @seed = @pop[0..$numElites-1];
    print STDERR "\tGeneration ", ($j+1), "\n";
    print "Generation ", ($j+1), "\n";
    print $_->{fit}, "\t", $_->{seq}, "\n" foreach @seed;
   }
}
exit;

sub random_seq{
  my $s;
  foreach (@site_aa){
    my @ss = sample(-set=>$_);
    $s .= shift @ss
  }
  return $s;  
}

sub mutate {
  my $s = shift;
  my @sitesToMutate = sample(-set => [keys %varSite], -sample_size => $mutRate );
  foreach my $site (@sitesToMutate) {
    my @aas = sample(-set => $varSite{$site});
    substr $s, $site-1, 1, shift @aas;
  }
  return $s
}

sub dis2ref {
  my $s = shift;
  my $sim = Bio::LocatableSeq->new(-id=>'sim', -seq=>$s);
  my $sum_d;
  my $max=0;
  foreach (@ref_seqs){
    my $pair = Bio::SimpleAlign->new();
    $pair->add_seq($sim);
    $pair->add_seq($_);
    my $d = 1 - $pair->percentage_identity()/100;
#    printf STDERR "%s(%0.3f)", ($_->id(), $d);
    $max = $d if $d>$max;
    $sum_d += $d
  }
  return $max . "\t" . $sum_d/$nRef
}

__ALN__
CLUSTAL W (1.81) multiple sequence alignment


L                      CNNSGKDGNASVNSADESVKGPNLVEISKKITDSNAVVIAVKEVETLLVSIDELA-KAIG
N                      CNNSGKDGNASTNSADESVKGPNLTEISKKITESNAVVLAVKEVAALLSSIDELA-KAIG
E                      CNNSGKDGNASANSADESVKGPNLTEISKKITESNAVVLAVKEVETLLASIDELATKAIG
G                      CNNSGKDGNASTNSADESVKGPNLAEISKKITESNAVVLAVKEVAALLSSIDELA-KAIG
B                      CNNSGKDGNTSANSADESVKGPNLTEISKKITDSNAVLLAVKEVEALLSSIDELA-KAIG
U                      CNNSGKDGNASANSADESVKGPNLAEISKKITESNAVVLAVKEVEALLASIDEIGSKAIG
O                      CNNSGKDGNTSANSADESVKGPNLTEISKKITESNAVVLAVKEVEALLSSIDELA-KAIG
D                      CNNSGKDGNTSANSADESVKGPNLTEISKKITDSNAVLLAVKEVEVLLSSIDELAKKAIG
H                      CNNSGKDGNASANSADESVKGPNLTEISKKITESNAVVLAVKEVETLLASINQLA-KAIG
J                      CNNSGKDGNTSANSADESVKGPNLTEISKKITESNAVVLAVKEIETLLASIDELATKAIG
I                      CNNSGKDGNTSANSADESVKGPNLTEISKKITESNAVVLAVKEVETLLTSIDELA-KAIG
A                      CNNSGKDGNTSANSADESVKGPNLTEISKKITDSNAVLLAVKEVEALLSSIDEIAAKAIG
K                      CNNSGKDGNTSANSADESVKGPNLTEISKKITESNAVVLAVKEIETLLASIDELATKAIG
C                      CNNSGKDGNASANSADESVKGPNLTEISKKITESNAVVLAVKEVETLLASIDELA-KAIG
M                      CNNSGKDGNTSANSADESVKGPNLTEISKKITESNAVVLAVKEVETLLASIDEVAKKAIG
T                      CNNSGKDGNASVNSADESVKGPNLTEISKKITESNAVVLAVKEVETLLASIDELANKAIG
F                      CNNSGKDGNTSANSADESVKGPNLTEISKKITESNAVVLAVKEIETLLSSIDELATKAIG
                       *********:*.************.*******:****::****: .** **:::. ****


L                      KKIEAGGTLGSDGAHNGSLLAGAYKIATEITAN---LSKLKASEDLKEKITKAKECSEKF
N                      KKINNN-GLDDVQNFNASLLAGAHTISKLVTEK---LSKLKNSEGLKEKIEDAKKCSDDF
E                      KKIGNN-GLEANQSKNTSLLSGAYAISDLIAEK---LNVLKN-EELKEKIDTAKQCSTEF
G                      KKIEQN-GLGADANHNTSLLAGAHEISTLIKQK---LDGLKN-EGLNKEIEAAKKCSAAF
B                      KKIKNDGSLDNEANRNESLLAGAYTISTLITQK---LSKLNGSEGLKEKIAAAKKCSEEF
U                      KRIQAN-GLQDLQGQNGSLLAGAYAISNLITQKINVLNGLKNSEELKEKINEAKGCSEKF
O                      KEIGAN-GLVNQANHNVSLLAGAYEISTLITEK---LSKLNGSEGLKEKIAAAKKCSEEF
D                      KKIDQNNALGTLDNHNGSLLAGAYAISALITEK---LSSIKDSGELKAEIEKAKKCSESF
H                      KKIDQNGTLGDDGGQNGSLLAGAYAISTVIIEK---LSTLKNVEELKEKITKAKDCSEKF
J                      KKIDNNAGLGAEVGQNGSLLAGAYAISTVIIEK---LSTLKNVEELKEKITKAKDCSEKF
I                      KKIKNDVSLDNEADHNGSLISGAYLISTLITKK---ISAIKDSGELKAEIEKAKKCSEEF
A                      KKIHQNNGLDTENNHNGSLLAGAYAISTLIKQK---LDGLKN-EGLKEKIDAAKKCSETF
K                      KKIQQNGGLAVEAGHNGTLLAGAYTISKLITQK---LDGLKNSEKLKEKIENAKKCSEDF
C                      KKIKNDVSLDNEADNNGSLIAGAYTISTLITQK---LSKLNGSEGLKEKIAEAKKCSEEF
M                      NLIAQN-GLNAGANQNGSLLAGAYVISTLIAEK---LDGLKNSEELKEKIEDAKKCNKAF
T                      QKIDQNNGLSVDAGHNGPLLAGAYAISALITEK---LNGLTISEELKEKIEKAKKCSTGF
F                      QKIDAN-GLGVQANQNGSLLAGAYAISTLITQK---LSAL-NSEDLKEKVAKVKKCSEDF
                       : *  .  *      * .*::**: *:  :  :   :. :     *: ::  .* *.  *


L                      TDKLKSENVALGKQD--ASDDDAKKAILKTHND-ITKGAKELKELSESVETLLKAAKEML
N                      TKKLQSSHAQLGVAGGATTDEEAKKAILRTNAI-KDKGADELEKLFKSVESLAKAAQDAL
E                      TNKLKSEHAVLGLDN--LTDDNAQRAILKKHAN-KDKGAAELEKLFKAVENLSKAAQDTL
G                      TKKLADSNADLGVAAGNATDDNAKRAILKTHGH-EDKGGKELKELSEAVKSLLKAAQAAL
B                      STKLKDNHAQLGIQG--VTDENAKKAILKANAAGKDKGVEELEKLSGSLESLSKAAKEML
U                      TKKLSESHADIGIQA--ATDANAKDAILKTNPT-KTKGAEELDKLFKAVENLSKAAKEML
O                      STKLKSSNAQLNQAN--ANDANAKAAILKTHNT-KDKGAEELVKLAESVAGLLKVAQEML
D                      TKKLSDNQAELGIEN--ATDDNAKKAILKTHNA-KDKGAEELVKLSESVAGLLKAAQAIL
H                      AGKLKNEHASLGKKD--ATDDDAKKAILKTHGN-TDKGAKELKDLSDSVESLVKAAKEML
J                      TKKLKDSHAELGKKD--ASDDDAKKAILKTNQA-NDKGAKELKELFEAVESLSKAAKEML
I                      TAKLKGEHTDLGKEG--VTDDNAKKAILKTNND-KTKGADELEKLFESVKNLSKAAKEML
A                      TNKLKEKHTDLGKEG--VTDADAKEAILKTNGT-KTKGAEELGKLFESVEVLSKAAKEML
K                      TKKLEGEHAQLGIEN--VTDENAKKAILITDAA-KDKGAAELEKLFKAVENLAKAAKEML
C                      TKKLKEKHTDLGKKD--ATDVHAKEAILKTNGT-KDKGAAELEKLFESVENLAKAAKEML
M                      TDKLKSSHAELGIANGAATDANAKAAILKTNGT-KDKGAQELEKLFESVKNLSKAAQETL
T                      TNKLKSGHAELGPVGGNATDENAKQAILKTHGN-VTKGAKELKDLSESVEALAKAAQAML
F                      TNKLKNGNAQLGLAA--ATDDNAKAAILKTNGT-NDKGAKELKDLSDSVESLVKAAQVML
                       : **   :. :.      .* .*: ***  .     **  ** .*  ::  * *.*:  *


L                      ANSVKELTSPVVAESPKKPLE
N                      ANSVNELTGPVVAETPKKPLE
E                      KNAVKELTSPIVAESPKKPLE
G                      ANSVQELTSPVVAETPKKPLE
B                      ANSVKELTSPVVVESPKKPLE
U                      ANSVKELTSPVVAESPKKPLE
O                      NNSVKELTSPVVAENPKKPLE
D                      ANSVKELTSPVVAESPKKPLE
H                      TNSVKELTSPVVAESPKKPLE
J                      NKSVKELTSPIVAESPKNPLE
I                      TNSVKELTSPVVAESPKKPLE
A                      ANSVKELTSPVVAESPKKPLE
K                      ANSVKELTSPIVAESPKNPLE
C                      SNSVKELTSPVVAESPKNPLE
M                      NNSVKELTSPVVAENPKKPLE
T                      TNSVKELTSPVVAETPKKPLE
F                      TNSVKELTSPVVAESPKKPLE
                        ::*:***.*:*.*.**:***

