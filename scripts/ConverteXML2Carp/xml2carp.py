
from xml.dom import minidom
#
def getCarpName(folder, extension):
    return folder+"."+extension;

    
def writecarp(fileout, data, prefix='', sufix=''):
    size = len(data);    
    if(size==0):
        return;    
    f = open(fileout, 'w');
    f.write("%d\n" % (size))
    for item in data:
        cols=len(item);
        if(prefix!=''):
            f.write(prefix);    
            f.write(' ');
        cont=0;
        for i in item:            
            f.write(i);
            if(cont<cols-1):
                f.write(' ');
            cont=cont + 1
        if(sufix!=''):
            f.write(' ');
            f.write(sufix);
        f.write('\n');
    f.close()
#
def xmlReadTags(root, tag, aAttr):
    aRet = [];    
    items = root[0].getElementsByTagName(tag)
    for item in items:
        aI = [];
        for i in aAttr:
            aI.append(item.attributes[i].value)
        aRet.append(aI)
    return aRet;

def readFibers(root, tag):
    aRet = [];    
    items = root[0].getElementsByTagName(tag)
    for item in items:        
        linha=[];
        eFiber  = item.getElementsByTagName("fiber")[0]
        eSheet  = item.getElementsByTagName("sheet")[0]
        eNormal = item.getElementsByTagName("normal")[0]
        #
        aFiber = eFiber.firstChild.nodeValue.split(',')
        aNormal= eNormal.firstChild.nodeValue.split(',')
        aSheet = eSheet.firstChild.nodeValue.split(',')
        #
        linha.append(aFiber[0])
        linha.append(aFiber[1])
        linha.append(aFiber[2])
        linha.append(aSheet[0])
        linha.append(aSheet[1])
        linha.append(aSheet[2])
        linha.append(aNormal[0])
        linha.append(aNormal[1])
        linha.append(aNormal[2])
        aRet.append(linha)
    return aRet;
###
#filein = 'in/.xml';
#nameout = 'out/bench3';
filein = 'in/mesh_patient_cardiax.xml';
nameout = 'out/mesh_patient_cardiax';

# parse an xml file by name
xml = minidom.parse(filein);
#
bounds = xml.getElementsByTagName('prescribed_displacement');
aBounds = xmlReadTags(bounds, 'node', ['id', 'direction', 'value']);
writecarp(getCarpName(nameout,"bound"), aBounds)
#
mesh = xml.getElementsByTagName('mesh');
mshNodes = mesh[0].getElementsByTagName('nodes');
aNodes = xmlReadTags(mshNodes, 'node', ['x', 'y', 'z']);
writecarp(getCarpName(nameout,"pts"), aNodes)
#
elements = mesh[0].getElementsByTagName('elements')
aElements = xmlReadTags(elements, 'element', ['v0', 'v1', 'v2', 'v3']);
writecarp(getCarpName(nameout,"elem"), aElements, 'Tt', '1')
#
press = mesh[0].getElementsByTagName('boundary')
aPress = xmlReadTags(press, 'element', ['v0', 'v1', 'v2']);
writecarp(getCarpName(nameout,"press"), aPress)
#
fibers = mesh[0].getElementsByTagName('element_data')
aFibers = readFibers(fibers, 'element');
writecarp(getCarpName(nameout,"fib"), aFibers)
