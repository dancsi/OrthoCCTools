import fastscoreCC
inter=fastscoreCC.Interaction()
inter.init_complete_score()
print inter.score_complete("AAAAAAAAAAAAA","AAAAAAAAAAAAA") 
print inter.score_complete("DEIQALEEENAQLEQENAALEEEIAQLEYG","DKIAQLKEKNAALKEKNQQLKEKIQALKYG") 
print inter.score_complete("AAA","AAA") 
print inter.score_complete("HHHHHHHHHHHHHHHH","HHHHHHHHHHHHHHHH") 