#!/usr/bin/env python
#
# Script to convert Abaqus mesh to Gmsh (.inp ---> .msh)
# Bernardo M. Rocha, 2014
#

import sys, time, os
import numpy as np

def abaqus2Gmsh (inpMesh, outputMesh):  
    """
    Given a mesh in ABAQUS .inp file format, convert it to Gmsh's format.
    It creates only nodes and elements sections.
    """
    
    ptslist, elelist =  list(), list()
    etagdict  = {} # holds { (elementId, surfId), ... }
    etagcount = 0
    eltype = []

    with open(inpMesh, "r") as f:
        searchlines = f.readlines()

    for i, line in enumerate(searchlines):

        if ("*NODE" in line):           
            for l in searchlines[i+1:]:
                if ("**" in l): break
                lm = map(float,l.split(',')[1:])
                ptslist.append(lm)                        

        if ("*ELEMENT" in line):
            elinfo = line.split(",")
            eltype.append( elinfo[1].strip()[5:] )
            for l in searchlines[i+1:]:
                if ("**" in l): break
                if ("*ELEMENT" in l): break
                enodes = map(int,l.split(',')[1:])
                elelist.append( enodes ) 
        
        if ("ELSET, ELSET=" in line):
            print("Converting element set")
            etagcount = etagcount + 1
            elsetinfo = line.split("=")[1]
            print(searchlines[i])
            print(searchlines[i+1])
            for l in searchlines[i+1:]:
                if ("**" in l): break
                if ("*SURFACE" in l or "SS" in l):
                    pass
                else:
                    l = l.strip()
                    l = l.rstrip(',')
                    etag = filter(int,l.split(','))
                    etag = map(int,l.split(','))
                    for i in range(len(etag)):
                        #etaglist[etagcount-1].append(etag[i])            
                        elemId = etag[i]
                        etagdict[elemId] = etagcount

    f.close()
    print 'Parsing OK'

    # split element list in case we have surface elements
    nvElem = 0
    nsElem = 0

    nl0 = len(elelist[0])
    idxsplit = -1
    for i in range(len(elelist)):
        nl = len(elelist[i])
        if(nl != nl0):
            idxsplit = i
            break

    if(idxsplit != -1):
        slst = elelist[:idxsplit]
        vlst = elelist[idxsplit:]
    else:
        slst = []
        vlst = elelist

    # create numpy arrays
    print 'Creating arrays'  
    vpts = np.array(ptslist)    

    numElem = []
    if ( len(slst) > 0 ): 
        selem = np.array(slst, dtype=np.int32)
        numElem.append( len(slst) )
    else:
        selem = None
    
    velem = np.array(vlst,dtype=np.int32)
    numElem.append( len(vlst) )                

    # write output file
    # basic header
    outputFile = open(outputMesh, 'w')
    outputFile.write('$MeshFormat\n2.0 0 8\n$EndMeshFormat\n')

    # nodes section
    numPts = np.shape(vpts)[0]
    print(" Number of nodes = %d" % numPts)
    outputFile.write('$Nodes\n')
    outputFile.write('%d\n' % numPts)
    for i in range(numPts):
        if( np.shape(vpts)[1] == 2):
            outputFile.write('%d %f %f 0.0\n' % (i+1,vpts[i,0],vpts[i,1]))
        elif ( np.shape(vpts)[1] == 3):
            outputFile.write('%d %f %f %f\n' % (i+1,vpts[i,0],vpts[i,1],vpts[i,2]))
    outputFile.write('$EndNodes\n')

    # elements section
    # first surface elements, then volume elements
    totalElem  = numElem[0] 
    if(len(numElem)>1): 
        totalElem += numElem[1]
    outputFile.write('$Elements\n')
    outputFile.write('%d\n' % totalElem)
    print " Number of surface elements %d" % len(slst)
    print " Number of volume elements %d" % len(vlst)
    print " Total number of elements %d" % totalElem


    # pre-process eltype list to remove duplicates
    s = set(eltype)
    eltype = list(s)

    # go
    if(len(eltype) > 0):
        for k in range(len(eltype)):
            elem = eltype[k]
            if  (elem=="STRI3"): numNodesElem=3
            elif(elem=="S4R"):   numNodesElem=4
            elif(elem=="C3D8R"): numNodesElem=8
            else:                numNodesElem = int(eltype[k][3])

            print(" Element type = %s" % elem)
            print("  Number of elements = %d" % numElem[k])
            print("  Number of nodes per element = %d" % numNodesElem)

            for i in range(numElem[k]):  
                surfindex = 99
                if (i+1 in etagdict.keys()):            
                    surfindex = etagdict[i+1]*10

                e = i
                if(k==1): e = i + numElem[0]

                if (numNodesElem == 8):
                    outputFile.write('%d 5 2 %d 2 ' % (e+1,surfindex))
                elif (numNodesElem == 4 and elem=='S4R'):
                    outputFile.write('%d 3 2 %d 2 ' % (e+1,surfindex))
                elif (numNodesElem == 4):
                    outputFile.write('%d 4 2 %d 2 ' % (e+1,surfindex))
                elif (numNodesElem == 3):
                    outputFile.write('%d 2 2 %d 2 ' % (e+1,surfindex))
                
                if (len(slst) > 1):
                    outputFile.write('%s\n' % ' '.join(map(str,slst[i])))
                else:
                    outputFile.write('%s\n' % ' '.join(map(str,vlst[i])))

        # end elemtype loop

        outputFile.write('$EndElements\n')
        outputFile.close()  
        print("Done")


if __name__ == "__main__":

    if (len(sys.argv) < 3):
        print("\n Usage: abaqus2gmsh <mesh.inp> <mesh.msh>\n")
        sys.exit(-1)

    # parse and check input
    inpfile = sys.argv[1]

    if (not os.path.isfile(inpfile)):
        print("\n Error: the input %s does not exist.\n" % (inpfile))
        sys.exit(-1)

    # convert
    abaqus2Gmsh(inpfile, sys.argv[2]) 

# end of main
