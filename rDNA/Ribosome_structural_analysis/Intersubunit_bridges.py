####  Variants with VF >= 5% and within 5A of intersubunit bridges  ####
# Original: Jan 2021; modified: Feb 2022                               #
# Daniel Sultanov                                                      #
# Prerequisites:                                                       #
# - PyMOL with pyhon v.3 and Pandas package                            #
# - file with rRNA variant annotation (rRNA_annotation.txt)            #
########################################################################
#Note: can be further optimized to run more efficiently if needed

#You can invoke the script directly in PyMOL through "@" (using .txt extension is preferable by .py also worked)
#to convert, simply change .py to .txt

###PyMOL part
#load ribosome structure
fetch 4v88-assembly1

#set background
bg_color white
set ray_opaque_background, off

#make OsHex and Mg selections (note: there is also Zn in the stucture)
select OsHex, r. OHX
select Mg, r. MG
hide spheres

#make selections for ribosome constituents
select SSU, c. A2
select LSU, c. A1 or c. A4 or c. A3

split_chains 4v88-assembly1, RP
#rename chains
set_name RP0001, 25S 
set_name RP0002, 18S
set_name RP0003, 5S
set_name RP0004, 5.8S

color yelloworange, 18S
color aquamarine, 25S
color slate, 5.8S
color deepteal, 5S

show_as surface
set transparency, 0.85

#set cartoon representation
cartoon tube
set cartoon_transparency, 0
hide cartoon


#Group ribosomal proteins
group RPs, RP*

#Ribosomal protein part of intersubunit bridges
#note: PyMOL might spit out "invalid character in identifier". To avoid, re-copy "set_name RP0024, S19" etc from PyMOL back to the script (open in text_editor)
set_name RP0024, S19
set_name RP0055, L19e
set_name RP0011, S6e
set_name RP0060, L24e
set_name RP0006, S1e
set_name RP0020, S15
set_name RP0039, L2
set_name RP0077, L41e
set_name RP0013, S8e
set_name RP0079, L43e

#rainbow color model
util.cbc(selection='RPs',first_color=1,quiet=1,legacy=0,_self=cmd)

##Add rRNA-rRNA intersubunit bridge contacts
hide spheres
#LSU
select LSU_b_RNA_RNA, resi 2257 or resi 2258 or resi 2259 or resi 2255 or resi 2262 or resi 2263 or resi 2264 or resi 2272 or resi 2195 or resi 2196 or resi 2275 or resi 2290 or resi 2291 or resi 2303 or resi 2302 or resi 2301 or resi 2292 or resi 2125 or resi 2124 or resi 2126 or resi 2305 or resi 2294 or resi 846 or resi 847 or resi 1935 or resi 2239 or resi 2240 in 25S
show_as spheres, LSU_b_RNA_RNA
set sphere_transparency, 0.6, LSU_b_RNA_RNA

#SSU
select SSU_b_RNA_RNA, resi 1646 or resi 1645 or resi 1644 or resi 1780 or resi 1004 or resi 1757 or resi 1758 or resi 1643 or resi 1759 or resi 996 or resi 994 or resi 995 or resi 1781 or resi 983 or resi 1746 or resi 1655 or resi 1747 or resi 1656 or resi 1657 or resi 1749 or resi 1659 or resi 971 or resi 972 or resi 973 or resi 630 or resi 628 or resi 1667 or resi 909 or resi 910 in 18S
show_as spheres, SSU_b_RNA_RNA
set sphere_transparency, 0.6, SSU_b_RNA_RNA

#RNA part of RNA-RP intersubunit bridges
#LSU
select LSU_b_RNA_RP, resi 1025 or resi 847 or resi 2536 or resi 2537 or resi 3354 or resi 3355 or resi 3353 or resi 3345 or resi 2107 in 25S
show_as spheres, LSU_b_RNA_RP
set sphere_transparency, 0.6, LSU_b_RNA_RP

#SSU
select SSU_b_RNA_RP, resi 1734 or resi 411 or resi 1670 or resi 1724 or resi 987 or resi 986 or resi 1013 or resi 1012 or resi 983 or resi 419 or resi 412 or resi 815 or resi 813 or resi 850 or resi 851 or resi 814 or resi 852 or resi 855 or resi 853 or resi 854 or resi 1641 or resi 1783 or resi 1773 or resi 1774 or resi 1784 or resi 1785 or resi 1112 or resi 1777 or resi 1778 or resi 1782 or resi 1642 or resi 1114 or resi 1775 or resi 1127 or resi 1779 or resi 1115 or resi 1126 or resi 1116 or resi 1125 or resi 1654 or resi 1118 or resi 1117 or resi 1655 in 18S
show_as spheres, SSU_b_RNA_RP
set sphere_transparency, 0.6, SSU_b_RNA_RP

#variant positions within 5A of intersubunit bridges
select milieu, byres all within 5 of (LSU_b_RNA_RNA or SSU_b_RNA_RNA or LSU_b_RNA_RP or SSU_b_RNA_RP or S19 or L19e or S6e or L24e or S1e or S15 or L2 or L41e or S8e or L43e)

##Extract all residues (check the presence of variant positions in those)
iterate milieu, print (model + " " + resi)

###Python part
python

from pymol import cmd
import pandas as pd

#load data
#specify path to rRNA_annotation.txt (i.e. ~/Desktop/work/database/rRNA_annotation.txt)
dfRNA = pd.read_table("~/path/to/rRNA_annotation/.txt", sep="\t", comment="#")
#only get variants with VF >= 5%
isVF=dfRNA["VF"]>=0.05
Var_positions = dfRNA[isVF] #indices are from initial df; VFP = variant frequency polymorphisms

##Color variable nucleotides and represent them as spheres
#Note: do for LSU and SSU separately
#for LSU
LSU = Var_positions[Var_positions["SUBUNIT"] == "LSU"]
for i in range(0,len(LSU.index)):
    res = str(LSU.iat[i,5])
    sele = "resi " + res + " in " + LSU.iat[i,1] + " and milieu"
    cmd.color("purpleblue", sele)
    cmd.show_as("spheres", sele)

#for SSU
SSU = Var_positions[Var_positions["SUBUNIT"] == "SSU"]
for i in range(0,len(SSU.index)):
    res = str(SSU.iat[i,5])
    sele = "resi " + res + " in " + SSU.iat[i,1] + " and milieu"
    cmd.color("purpleblue", sele)
    cmd.show_as("spheres", sele)

python end

##PyMOL
show cartoon, S19
show cartoon, L19e
show cartoon, S6e
show cartoon, L24e
show cartoon, S1e
show cartoon, S15
show cartoon, L2
show cartoon, L41e
show cartoon, S8e
show cartoon, L43e
