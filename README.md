## RQC (RNA Quality Control)

RQC helps you analyse the quality of your BAM files. Input BAM files must be sorted and indexed. RQC can do the following:

* Summarise read PHRED scores for one or multiple BAM files


## SETUP
```
# set up and enter virtual environment
python3 -m venv .venv
source .venv/bin/activate
pip install -r requirements.txt

# run rqc
python rqc.py -q 7 -m 5 plot "36C1 Pfal" /Users/joshualevendis/Downloads/36C1_to_pfal.sorted.bam

# exit the virtual environment
deactivate
```
If you have a file that looks like this:
samples.txt
```
36C1_Pfal /Users/joshualevendis/Downloads/36C1_to_pfal.sorted.bam 
36C1_Yeast /Users/joshualevendis/Downloads/36C1_to_yeast_sorted.bam 
36C1_Human /Users/joshualevendis/Downloads/36C1_to_humans.sorted.bam
```
You can run something like this:
```
python rqc.py -q 7 -m 5 plot $(cat samples.txt)
```