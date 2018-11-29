import yasara as y
import os

def dock(macrotarget):
    y.ApplyMacro(os.path.join(y.info.dir,"mcr/dock_run.mcr"), \
                 targets=macrotarget, \
                 remove="FromUnderscore", \
                 newextension="")

#y.info.mode='txt'
target = "a_receptor.pdb"
print("Docking "+target)
dock(target)
print("Finished "+target) 
