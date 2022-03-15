SRCFOLDER:=/import/freenas-m-03-geodynamics/jhayek/TEAR/Results/T3/Runs
SRCFILES:=/P1_P1_100x100/SlipAtReceiver.txt /P2_P1_050x050/SlipAtReceiver.txt /P3_P1_025x025/SlipAtReceiver.txt /P4_P3_025x025/SlipAtReceiver.txt  /P5_P3_050x050/SlipAtReceiver.txt  /P6_P3_100x100/SlipAtReceiver.txt
#SRCFILES:=/P1_P1_100x100/SlipAtReceiver.txt /P2_P1_100x100/SlipAtReceiver.txt


CODEFILES := $(addprefix $(SRCFOLDER),$(SRCFILES))  
OUTFILES := $(addprefix .,$(SRCFILES))
#SPECIFICFILE := $(word 4,$(CODEFILES))
SPECIFICFILE := /import/freenas-m-03-geodynamics/jhayek/TEAR/Results/PM/Runs/P1_P1_100x100/SlipAtReceiver.txt


all:
	scp -r jhayek@wolfel:${SPECIFICFILE} .


