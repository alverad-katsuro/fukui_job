import os

orbitais = os.popen("awk '/Summary of Natural Population Analysis:/{f=1;next} /=======================================================================/{f=0} f'  3_stage_rank_0.log")