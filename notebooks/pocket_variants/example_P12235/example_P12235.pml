# load AF-P12235-F1-filled_v4.cif
# /select/
# save AF-P12235-F1-filled_v4_ligand.pdb, (sele)

delete all

bg_color white
set_color tab_blue, [31,119,180]
set_color tab_orange, [255,127,14]
set_color tab_green, [44,160,44]
set_color tab_red, [214,39,40]
set_color tab_purple, [148,103,189]
set_color tab_brown, [140,86,75]
set_color tab_pink, [227,119,194]
set_color tab_gray, [127,127,127]
set_color tab_olive, [188,189,34]
set_color tab_cyan, [23,190,207]

delete all
load ~/euler-home/project/22.12_pocketomes/results/23.04_bindfunc/af2/P1/22/35/P12235-F1.pdb
color grey50, P12235-F1
load ~/euler-home/project/22.12_pocketomes/results/23.04_bindfunc/af2.obabel_hxr.autosite/P1/22/35/P12235-F1/P12235-F1_cl_001.pdb
hide everything, P12235-F1_cl_001
show surface, P12235-F1_cl_001
set transparency, 0.5, P12235-F1_cl_001
color tab_blue, P12235-F1_cl_001
color tab_blue, P12235-F1 & resid 4+9+17+20+24+25+27+32+43+56+63+81+84+89+92+108+113+115+117+118+123+125+139+141+145+148+149+157+164+167+172+182+189+192+193+209+216+217+218+221+227+229+237+247+252+264+283+287+288+289+294
color tab_red, P12235-F1 & resid 33+80+278
color tab_red, P12235-F1 & resid 123
show sticks, P12235-F1 & resid 33+80+278
show sticks, P12235-F1 & resid 123

load ~/euler-home/project/af2genomics/notebooks/pocket_variants/example_P12235/AF-P12235-F1-filled_v4_ligand.pdb
color tab_pink, AF-P12235-F1-filled_v4_ligand

# >get_view
### cut below here and paste into script ###
set_view (\
    -0.650603294,   -0.389876902,    0.651697516,\
    -0.708973646,    0.619374156,   -0.337240577,\
    -0.272161335,   -0.681444287,   -0.679379165,\
     0.000000000,    0.000000000, -191.857254028,\
    -0.766267776,    5.128052711,   -3.638110161,\
    73.120101929,  310.594329834,  -20.000000000 )
### cut above here and paste into script ###

set ray_opaque_background, 0
png ~/euler-home/project/af2genomics/notebooks/pocket_variants/example_P12235/example_P12235_.png, width=3cm, height=3cm, dpi=600, ray=1