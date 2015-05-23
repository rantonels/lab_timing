#!/usr/bin/python

import glob
import os
import sys
import subprocess
import re
import pickle

HVpattern = re.compile(r'(2R)?(\d{4})HV(\d+)GC\d+min\.root')

CLIST = []

PLIST = []


##costanti per le sigma temporali
t_binsize = 33
t_binnum = 16
TLIST = [
            "anticipo1ns",
            "ritardo0ns",
            "ritardo1ns",
            "ritardo2ns"
        ]

ritardi = [
    5.4,
    6.8,
    8.2,
    9.7
           ]



def precheck():
#	os.chdir(os.path.realpath(__file__))
	if not os.path.isdir("data"):
		print >> sys.stderr, "analisi.py ERRORE: cartella dati non trovata"
		sys.exit()

class HV():
    pass


def istogrammi():
        global PLIST,CLIST
        print "Calcolo istogrammi."
	print
	count = 0
        #CLIST = []
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


                        mo = HVpattern.match(file)
                        if mo:
                            nhv = HV()
                            nhv.fname = file
                            nhv.voltaggio   = int(mo.group(2))
                            nhv.gain        = int(mo.group(3))
                            if mo.group(1) == '2R':
                                nhv.ch = 3
                            else:
                                nhv.ch = 2
                            CLIST.append(nhv)
                            print "Trovato file HV con fname %s, voltaggio %d, gain %d, canale %d"%(nhv.fname,nhv.voltaggio,nhv.gain,nhv.ch)
        
        print "trovati %d HV"%len(CLIST)

	
	print "calcolati",count,"x4 istogrammi."
			
	open("tmp/donehisto", 'a').close()


#for file in os.listdir("data"):


#CLIST = [#lista dei file del quale fare i profili compton
#	"tmp/giovedi11_1_h2",
#	"tmp/giovedi11_1_h3",
#        "tmp/1430HV100GC10min_h3",
#        "tmp/1530HV100GC10min_h3",
#        "tmp/1630HV100GC10min_h3",
#        "tmp/1730HV100GC25min_h3",
#        "tmp/1930HV40GC10min_h3",
#        "tmp/2R1420HV100GC10min_h3",
#        "tmp/2R1520HV100GC10min_h3",
#        "tmp/2R1620HV100GC10min_h3",
#        "tmp/2R1720HV100GC25min_h3",
#        "tmp/2R1830HV40GC10min_h3",
#        "tmp/2R1930HV40CG10min_h3"
#]

def comptonfits():
        global PLIST,CLIST
	print "fit profili compton (molto lungo)"
	print

	fcount = 1
	for hvfile in CLIST:
                file = "tmp/" + os.path.splitext(hvfile.fname)[0] + "_h" + str(hvfile.ch)
                print file
                print "* FIT",fcount,"su",len(CLIST),"-",file
		
		if not os.path.exists(file):
			print >> sys.stderr, "analisi.py ERRORE: il file "+file+"necessario per il fit compton non esiste"
			sys.exit()
		
                gainstr = str( float(hvfile.gain)/100.0        )

                call = "bin/comptonfit analisi "+file # + " " +gainstr
		print call
		ret = subprocess.call(call,shell=True)
		if (ret>0):
			print >> sys.stderr, "analisi.py ERRORE: comptonfit e' stato terminato con errore"
			sys.exit()
		
                PLIST.append( (file, file + ".ccurve") )
                
                fcount +=1




	
	open("tmp/donefits", 'a').close()

def generate_gnuplot_script():
    global PLIST,CLIST
    print "generazione script gnuplot"
    fl = open('tmp/fitplots','w')

    fl.write('''set term png\n\n''')

    fl.write('''set xrange [0:1000]\n''')

    for f,fc in PLIST:
        fl.write('''set output "%s"\n\n'''%(f+".png"))
        strinka = '''plot "%s" u 0:1 with lines, "%s" u 1:2 with lines\n''' % (f,fc)
        fl.write(strinka)

    fl.close()


def sigmas():
    global PLIST,CLIST
    print "calcolo delle sigma di risoluzione"

    for h in CLIST:
        fl = open("tmp/"+os.path.splitext(h.fname)[0]+"_h"+str(h.ch) + ".cfit")
        com,e1,e2,k,y,sigma,S = fl.readlines()
        fl.close()
        if com[0] != '#':
            print "ERROR WRONG HEADER"
            sys.exit(1)

        h.sigma_E = float(sigma)

        print h.voltaggio,h.ch,h.sigma_E

def timesig():
    global TLIST

    print
    print "calcolo delle sigma temporali"

    for i in [1,2]:

        print "CANALE %d"%i

        call = "bin/timegaussians data/giovedi11_1.root %d %d %d tmp/giovedi11_1.root.time.R%d"%(t_binnum,t_binsize,i,i)
        print call
        ret = subprocess.call(call,shell=True)
        if (ret>0):
            print >> sys.stderr, "ERRORE giovedi"
            sys.exit()

        for f in TLIST:
            print f

            call = "bin/timegaussians data/"+f+".root %d %d %d"%(1,1000,i)+" "+"tmp/"+f+".time.R%d"%i
            print call
            ret = subprocess.call(call,shell=True)
            if (ret>0):
                print >> sys.stderr, "analisi.py ERRORE: timegaussians e' stato terminato con errore"
                sys.exit()


def bundletimesigs():
    global TLIST
    

    print
    print "raccoglimento sigma temporali..."


    #j = open("tmp/timesigmas",'w')
    
    #for i in range(len(TLIST)):
    #    f = TLIST[i]
    #    rit = ritardi[i]
    #    counter = 0
    #    for l in open("tmp/"+f+".time",'r').readlines():
    #        v,verr = [float(x) for x in l.split(',')]

    #        if counter != 0:
    #            j.write("%f\t%d\t%f\t%f\n"%(rit,counter,v,verr))

    #        counter += 1

    #j.close()

    for ch in [1,2]:
        j = open("tmp/timesigmasR%d"%ch,'w')
        for i in range(len(TLIST)):
                f = TLIST[i]
                rit = ritardi[i]
                counter = 0
                for l in open("tmp/"+f+".time.R%d"%ch,'r').readlines():
                    v,verr = [float(x) for x in l.split(',')]

                    j.write("%f\t%f\t%f\n"%(rit,v,verr))

                    counter+=1
        j.close()


def analisi():
        global PLIST,CLIST
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
                pickle.dump(PLIST,open("tmp/PLIST",'w'))
                pickle.dump(CLIST,open("tmp/CLIST",'w'))
	else:
		print "fit gia' fatti. (Pulire tmp/ per ripetere il calcolo dei fit)"
	


	print
	
        PLIST = pickle.load(open("tmp/PLIST",'r'))
        CLIST = pickle.load(open("tmp/CLIST",'r'))
        generate_gnuplot_script()

        sigmas()

        timesig()

        bundletimesigs()

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
