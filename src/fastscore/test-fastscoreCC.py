import fastscoreCC
inter=fastscoreCC.Interaction()
inter.init_complete_score()
print inter.score_complete("AAAAAAAAAAAAA","AAAAAAAAAAAAA",0) 
print inter.score_complete("DEIQALEEENAQLEQENAALEEEIAQLEYG","DKIAQLKEKNAALKEKNQQLKEKIQALKYG",0) 
print inter.score_complete("AAA","AAA",0) 
print inter.score_complete("HHHHHHHHHHHHHHHH","HHHHHHHHHHHHHHHH",0) 