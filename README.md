# SETUP

python3 -m venv .venv
source .venv/bin/activate
pip install -r requirements.txt

python rqc.py -q 0 -m 5 plot "36C1 Pfal" /Users/joshualevendis/Downloads/36C1_to_pfal.sorted.bam "36C1 Yeast" /Users/joshualevendis/Downloads/36C1_to_yeast_sorted.bam "36C1 Human" /Users/joshualevendis/Downloads/36C1_to_humans.sorted.bam

deactivate