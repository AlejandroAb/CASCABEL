#!/usr/bin/perl
use strict;

#Este programa calcula valores de N, tamanio promedio de contig, contig mas largo y numero de contigs contenidos en dicho N a partir de secuencias ensambladas con base en el tañanio de cada contig

#El archivo de entrada es un archivo de texto en formato FASTA descomprimido donde cada identificador corresponde a un contig

my @longitudes ;
my @long ;
my @ordenado ; 
my @parametros ;

my $elemento ;
my $elementoS ;
my $elementoN ;
my $linea ; 
my $N ;
my $totalBases ;
my $valorMedio ;
my $sumaElem ;
my $contadorN ;
my $posicionN ;
my $tamanoProm ;
my $contcontig ;
my $Totalong ;
my $longlinea ;
my $ultVal ;

scalar(@ARGV) == 1 || die "usage: $0 <assembly.fasta>\n";

######  ABRIENDO EL ARCHIVO Y CALCULANDO LAS LONGITUDES DE CADA CONTIG ENSAMBLADO  ###################3
my $file = $ARGV[0] ;
# '/home/ales/Documents/4_semestre/manejoSeq/SeqArch/scaftigEFD_K47' ;
open (ARCH, $file) or die ("No puedo abrir el archivo $file") ;
  while (<ARCH>) {#abre while ARCH
  chomp ;
  my $linea = $_ ;
	if ($linea =~ /^>/) {#abre if
	$contcontig++ ;
	push (@longitudes, $Totalong) ;
	$Totalong = 0 ;
	   }else{
	   $longlinea = length($linea) ;
	   $Totalong += $longlinea ;
	   }#cierra else
  $ultVal = $Totalong ;
  }#cierra while ARCH
push (@longitudes, $ultVal) ;

##################### ORDENANDO LOS ELEMENTOS DEL ARREGLO QUE CONTIENE LAS LONGITUDES  #############
@ordenado = sort { $b <=> $a } (@longitudes) ;
  foreach $elemento (@ordenado) {
  $totalBases += $elemento ;
  }#cierra foreach
$tamanoProm = int $totalBases/$contcontig ;
print "Numero total de contigs: $contcontig\n" ;
print "Tamanio promedio de los contigs: $tamanoProm\n" ;
print "El contig mas largo mide: $ordenado[0]\n" ;
print "El contig mas corto mide: $ordenado[-2]\n" ;

#####################  ESPECIFICANDO LOS VALORES DE N ##############################################
@parametros = (10,20,30,40,50,60,70,80,90,100) ;
  foreach $elementoS (@parametros) {
  $valorMedio = 0 ;
  $sumaElem = 0 ;
  $contadorN = 0 ;
  &calcN($elementoS) ;
  }# cierra foreach

exit;

####################  SUBRUTINA PARA CALCULAR LOS N ####################

sub calcN {
$N = $_[0] ;
$valorMedio = ($N/100)*$totalBases ; 
   foreach $elementoN (@ordenado) {
	unless ($sumaElem >= $valorMedio) {
	$sumaElem += $elementoN ;
	$contadorN ++ ;
	}#cierra unless
   }#cierra foreach
$posicionN = $contadorN -1 ;
return print "N$N = $ordenado[$posicionN]	contiene $contadorN contigs	total = $sumaElem bases\n" ;
}# cierra sub

