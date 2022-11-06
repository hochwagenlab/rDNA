####              Map polymorphic sites with VF >= 5%               ####
# Original: Dec 2021; modified: Feb 2022                               #
# Daniel Sultanov                                                      #
# Prerequisites:                                                       #
# - PyMOL with pyhon v.3 and Pandas package                            #
# - file with rRNA variant annotation (rRNA_annotation.txt)            #
########################################################################

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

#set cartoon representation
cartoon tube, RP*
set cartoon_transparency, 0.5
set cartoon_transparency, 0.8, RP*

#Group ribosomal proteins
group RPs, RP*

#rainbow color model
util.cbc(selection='RPs',first_color=1,quiet=1,legacy=0,_self=cmd)

###Python part
python

from pymol import cmd
import pandas as pd

#load data
#specify path to rRNA_annotation.txt (i.e. ~/Desktop/work/database/rRNA_annotation.txt)
dfRNA = pd.read_table("~/Desktop/work/database/rRNA_annotation.txt", sep="\t", comment="#")
#only get variants with VF >= 5%
isVF=dfRNA["VF"]>=0.05
Var_positions = dfRNA[isVF] #indices are from initial df; VFP = variant frequency polymorphisms

##Color variable nucleotides and represent them as spheres
#Blue if nucleotide does not interact with RP
#Red if nucleotide interacts with RP
#Note: do for LSU and SSU separately
#for LSU
LSU = Var_positions[Var_positions["SUBUNIT"] == "LSU"]
for i in range(0,len(LSU.index)):
    if pd.isna(LSU.iat[i,11]): #pd.isna returns TRUE if cell is NaN (no RP interaction)
        res = str(LSU.iat[i,5])
        #has to reference rRNA because COORDINATE values overlap!
        sele = "resi " + res + " in " + LSU.iat[i,1] #value must correspond to selection made above with /select/ in PyMOL
        cmd.color("blue", sele)
        cmd.show_as("spheres", sele)
    else:
        res = str(LSU.iat[i,5])
        sele = "resi " + res + " in " + LSU.iat[i,1] 
        cmd.color("red", sele)
        cmd.show_as("spheres", sele)

#for SSU
SSU = Var_positions[Var_positions["SUBUNIT"] == "SSU"]
for i in range(0,len(SSU.index)):
    if pd.isna(SSU.iat[i,11]):
        res = str(SSU.iat[i,5])
        sele = "resi " + res + " in " + SSU.iat[i,1]
        cmd.color("blue", sele)
        cmd.show_as("spheres", sele)
    else:
        res = str(SSU.iat[i,5])
        sele = "resi " + res + " in " + SSU.iat[i,1]
        cmd.color("red", sele)
        cmd.show_as("spheres", sele)

python end

