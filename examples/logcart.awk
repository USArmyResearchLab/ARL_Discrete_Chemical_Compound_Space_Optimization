#
#  This awk script will take the cartesian coordinates from the last standard orientation in a g03 log file and make a pdb file  
#
#  This script was created by Lisa M. Perez at the Laboratory for Molecular Simulation and is provided as is.  This script
#     has not been fully tested.  Use at your own risk.
#
#  Please e-mail Lisa M. Perez at mouse@mail.chem.tamu.edu with problems.
#
#  To use this script you should make an alias in your .tcshrc, .cshrc or .login file:  
#
#  alias logtopdb 'awk -f /home/mouse/log_to_pdb.awk \!*.log > \!*.pdb'
#    
#  This alias implies that you put log_to_pdb.awk in your home directory and that your Gaussian output files end in .log  
#     You will type logtopdb filename to execute.  The script will create a file named filename.pdb
#
#  last updated 03/05/04
#

BEGIN{i=0} 


# Get the coordinates from the output file

/Standard/,/Rotational/ {
if (/       0     /) {i=$1;an[i]=$2;x[i]=$4;y[i]=$5;z[i]=$6
} 
}

END {
printf "   %d\n\n",i

for (k=1;k<=i;k++) {

if (an[k]=="1") {as[k]="H"};if (an[k]=="2") {as[k]="He"};if (an[k]=="3") {as[k]="Li"};if (an[k]=="4") {as[k]="Be"};if (an[k]=="5") {as[k]="B"}
if (an[k]=="6") {as[k]="C"};if (an[k]=="7") {as[k]="N"};if (an[k]=="8") {as[k]="O"};if (an[k]=="9")  {as[k]="F"};if (an[k]=="10") {as[k]="Ne"}
if (an[k]=="11") {as[k]="Na"};if (an[k]=="12") {as[k]="Mg"};if (an[k]=="13") {as[k]="Al"};if (an[k]=="14") {as[k]="Si"};if (an[k]=="15") {as[k]="P"}
if (an[k]=="16") {as[k]="S"};if (an[k]=="17") {as[k]="Cl"};if (an[k]=="18") {as[k]="Ar"};if (an[k]=="19") {as[k]="K"};if (an[k]=="20") {as[k]="Ca"}
if (an[k]=="21") {as[k]="Sc"};if (an[k]=="22") {as[k]="Ti"};if (an[k]=="23") {as[k]="V"};if (an[k]=="24") {as[k]="Cr"};if (an[k]=="25") {as[k]="Mn"}
if (an[k]=="26") {as[k]="Fe"};if (an[k]=="27") {as[k]="Co"};if (an[k]=="28") {as[k]="Ni"};if (an[k]=="29") {as[k]="Cu"};if (an[k]=="30") {as[k]="Zn"}
if (an[k]=="31") {as[k]="Ga"};if (an[k]=="32") {as[k]="Ge"};if (an[k]=="33") {as[k]="As"};if (an[k]=="34") {as[k]="Se"};if (an[k]=="35") {as[k]="Br"}
if (an[k]=="36") {as[k]="Kr"};if (an[k]=="37") {as[k]="Rb"};if (an[k]=="38") {as[k]="Sr"};if (an[k]=="39") {as[k]="Y"};if (an[k]=="40") {as[k]="Zr"}
if (an[k]=="41") {as[k]="Nb"};if (an[k]=="42") {as[k]="Mo"};if (an[k]=="43") {as[k]="Tc"};if (an[k]=="44") {as[k]="Ru"};if (an[k]=="45") {as[k]="Rh"}
if (an[k]=="46") {as[k]="Pd"};if (an[k]=="47") {as[k]="Ag"};if (an[k]=="48") {as[k]="Cd"};if (an[k]=="49") {as[k]="In"};if (an[k]=="50") {as[k]="Sn"}
if (an[k]=="51") {as[k]="Sb"};if (an[k]=="52") {as[k]="Te"};if (an[k]=="53") {as[k]="I"};if (an[k]=="54") {as[k]="Xe"};if (an[k]=="55") {as[k]="Cs"}
if (an[k]=="56") {as[k]="Ba"};if (an[k]=="57") {as[k]="La"};if (an[k]=="58") {as[k]="Ce"};if (an[k]=="59") {as[k]="Pr"};if (an[k]=="60") {as[k]="Nd"}
if (an[k]=="61") {as[k]="Pm"};if (an[k]=="62") {as[k]="Sm"};if (an[k]=="63") {as[k]="Eu"};if (an[k]=="64") {as[k]="Gd"};if (an[k]=="65") {as[k]="Tb"}
if (an[k]=="66") {as[k]="Dy"};if (an[k]=="67") {as[k]="Ho"};if (an[k]=="68") {as[k]="Er"};if (an[k]=="69") {as[k]="Tm"};if (an[k]=="70") {as[k]="Yb"}
if (an[k]=="71") {as[k]="Lu"};if (an[k]=="72") {as[k]="Hf"};if (an[k]=="73") {as[k]="Ta"};if (an[k]=="74") {as[k]="W"};if (an[k]=="75") {as[k]="Re"}
if (an[k]=="76") {as[k]="Os"};if (an[k]=="77") {as[k]="Ir"};if (an[k]=="78") {as[k]="Pt"};if (an[k]=="79") {as[k]="Au"};if (an[k]=="80") {as[k]="Hg"}
if (an[k]=="81") {as[k]="Tl"};if (an[k]=="82") {as[k]="Pb"};if (an[k]=="83") {as[k]="Bi"};if (an[k]=="84") {as[k]="Po"};if (an[k]=="85") {as[k]="At"}
if (an[k]=="86") {as[k]="Rn"};if (an[k]=="87") {as[k]="Fr"};if (an[k]=="88") {as[k]="Ra"};if (an[k]=="89") {as[k]="Ac"};if (an[k]=="90") {as[k]="Th"}
if (an[k]=="91") {as[k]="Pa"};if (an[k]=="92") {as[k]="U"};if (an[k]=="93") {as[k]="Np"};if (an[k]=="94") {as[k]="Pu"};if (an[k]=="95") {as[k]="Am"}
if (an[k]=="96") {as[k]="Cm"};if (an[k]=="97") {as[k]="Bk"};if (an[k]=="98") {as[k]="Cf"};if (an[k]=="99") {as[k]="Es"};if (an[k]=="100") {as[k]="Fm"}
if (an[k]=="101") {as[k]="Md"};if (an[k]=="102") {as[k]="No"};if (an[k]=="103") {as[k]="Lr"}

printf "%5s %8.3f%8.3f%8.3f \n",as[k],x[k],y[k],z[k]
}

}
