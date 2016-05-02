f1 = open("./2new0.gro",'r')
f2 = open("./new0.gro",'r')
f3 = open("selectedAtoms.txt",'w')

f1.readline()
f2.readline()

n1 = int( f1.readline())
n2 = int( f2.readline())

ans = []
s1 = str(f1.readline())
sSplit1 = s1.split()

for i in range (1,n2+1):    
    s2 = str(f2.readline())
    sSplit2 = s2.split()
    if(sSplit2[1] == sSplit1[1] and sSplit2[3] == sSplit1[3] ):
        ans.append(i)
        f3.write(str(i))
        if(len(ans) == n1):
            f3.write("\n")
            break            
        f3.write(" ")
        s1 = str(f1.readline())
        sSplit1 = s1.split()

print (ans)

f1.close();f2.close();f3.close()
