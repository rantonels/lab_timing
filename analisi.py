#!/usr/bin/python

import glob
import os
import sys
import subprocess

def precheck():
#	os.chdir(os.path.realpath(__file__))
	if not os.path.isdir("data"):
		print >> sys.stderr, "analisi.py ERRORE: cartella dati non trovata"
		sys.exit()


def istogrammi():
	print "Calcolo istogrammi."
	print
	count = 0
	for file in os.listdir("data"):
	    if (file.endswith(".root")) and (file != "1830HV40CG10min.root"):   #e' importante escludere questo file corrotto
			print "*",file
			calls = ["bin/root2csv", "data/"+file, "tmp/"+os.path.splitext(file)[0]] 
			print "$"," ".join(calls)
			ret = subprocess.call(" ".join(calls), shell=True)
			if (ret > 0):
				print >> sys.stderr, "analisi.py ERRORE: root2csv e' stato terminato con errore."
				sys.exit()
			count +=1
	
	print "calcolati",count,"x4 istogrammi."
			
	open("tmp/donehisto", 'a').close()

CLIST = [#lista dei file del quale fare i profili compton
	"tmp/giovedi11_1_h2",
	"tmp/giovedi11_1_h3",

        "1430HV100GC10min_h3",
        "1530HV100GC10min_h3",
        "1630HV100GC10min_h3",
        "1730HV100GC25min_h3",
        "1930HV40GC10min_h3",
        "2R1420HV100GC10min_h3",
        "2R1520HV100GC10min_h3",
        "2R1620HV100GC10min_h3",
        "2R1720HV100GC25min_h3",
        "2R1830HV40GC10min_h3",
        "2R1930HV40CG10min_h3"
]

def comptonfits():
	print "fit profili compton (molto lungo)"
	print
	fcount = 1
	for file in CLIST:
		print "* FIT",fcount,"su",len(CLIST),"-",file
		
		if not os.path.exists(file):
			print >> sys.stderr, "analisi.py ERRORE: il file "+file+"necessario per il fit compton non esiste"
			sys.exit()
		
		call = "bin/comptonfit analisi "+file
		print call
		ret = subprocess.call(call,shell=True)
		if (ret>0):
			print >> sys.stderr, "analisi.py ERRORE: comptonfit e' stato terminato con errore"
			sys.exit()
		fcount +=1
	


def analisi():
	print
	print "inizio analisi completa."
	print
	if not os.path.exists("tmp/donehisto"):
		istogrammi()
	else:
		print "istogrammi gia' fatti. (Pulire tmp/ per ripetere il calcolo degli istogrammi)"
		
	print
	
	if not os.path.exists("tmp/donefits"):
		comptonfits()
	else:
		print "fit gia' fatti. (Pulire tmp/ per ripetere il calcolo dei fit)"
		
	print
		
	print "FINE ANALISI, senza errori :D"
	

def aiuto():
	print "analisi()	esegue analisi completa"
	print "istogrammi()	calcola gli istogrammi"


print "Programma di analisi per l'esperienza \"timing\" "

precheck()

if __name__ == "__main__":
	print "nota: per eseguire solo parte delle operazioni di analisi, importare questo script come modulo in un interpreter python e seguire le istruzioni."
	analisi()
else:
	print "analisi() esegue l'analisi completa. aiuto() mostra i comandi disponibili"
