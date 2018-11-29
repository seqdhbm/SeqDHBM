import yasara as y

y.info.mode='txt'

y.LoadPDB("1crn")
# List all arginine residues
y.ListRes("Arg")
# Exit YASARA
y.Exit()
