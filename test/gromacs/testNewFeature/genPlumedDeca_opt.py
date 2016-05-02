import sys
def myRestraint(num_atoms, giv_id, num_sim, max_steps, force_const, stride, rmsdThres):
    myRes = "myRestraint: MyForce_opt ARG="
    myRes += "p1.x,p1.y,p1.z"
    for i in range(2, num_atoms+1):
        myRes += ",p" + str(i)+ ".x,p" + str(i) + ".y,p" + str(i) + ".z"
    #myRes += " ID=" + str(giv_id) + " NUM_SIM="+str(num_sim) + " FORCE_CONST="+str(force_const) +" NUM_STRIDE="+str(stride)+ " NUM_TURNS="+str(turns)+ " INIT_TURNS="+ str(init_turns) + "\n"
    myRes += " ID=" + str(giv_id) + " MAX_STEPS="+str(max_steps) +" FORCE_CONST="+str(force_const) +" NUM_STRIDE="+str(stride)+ " RMSD_THR="+ str(rmsdThres)+" WRITE_RMSD=0\n"
    return myRes

def genPlumed(num_atoms, num_sim, max_steps ,force_const=0.001, stride=1, rmsdThres=0.01) :
    atom_pos = ""
    finp = open("./filesImplicit2/selectedAtoms.txt",'r');
    inp = str(finp.readline())
    inpList = inp.split()
    finp.close()
    for i in range(1,len(inpList)+1):
        atom_pos += "p"+str(i)+": POSITION ATOM="+str(inpList[i-1]) + "\n"
    for i in range(num_sim):
        f = open("plumed.dat."+str(i), "w")
        f.write(atom_pos)
        f.write(myRestraint(len(inpList), i, num_sim, max_steps, force_const, stride, rmsdThres))
        #f.write("END2END: DISTANCE ATOMS=1,95\n")
        f.write("END2END: DISTANCE ATOMS="+inpList[0]+","+inpList[-1]+"\n")
        #f.write("ENERGY LABEL=PE\n")
        #f.write("PRINT STRIDE="+str(stride)+ " ARG=myRestraint.force2,myRestraint.rmsdAll,myRestraint.bias,END2END, FILE=COLVAR_opt\n")
        f.write("PRINT STRIDE="+ "100" + " ARG=myRestraint.force2,myRestraint.rmsdAll,myRestraint.bias,END2END, FILE=COLVAR_opt\n")
        f.close()

if len(sys.argv) == 6 :
    genPlumed(int(sys.argv[1]), int(sys.argv[2]), int(sys.argv[3]), float(sys.argv[4]), rmsdThres=float(sys.argv[5]))
elif len(sys.argv) == 5 :
    genPlumed(int(sys.argv[1]), int(sys.argv[2]), int(sys.argv[3]), float(sys.argv[4]))
else :
    genPlumed(int(sys.argv[1]), int(sys.argv[2]), int(sys.argv[3]))
