sp = "sp1  sp2  sp3  sp4  sp5  sp6  sp7  sp8  sp9 sp10 sp11 sp12 sp14 sp15 sp16 sp17 sp18 sp19 sp20 sp21 sp22 sp23"
sp_names = sp.replace("sp", "\",\"sp")
sp_names += "\""
sp_names = sp_names.replace(" ", "")

print(sp_names)


states = "0    0    0    1    0    0    1    0    0    0    0    0    0    1    1    0   0    0    1    1    0    0"
states = states.split()
states = [int(s) + 1 for s in states]

sp = sp.split()

state_str = "\""
for i, j in zip(sp, states):
    state_str += i + "=" + str(j) + ","

state_str += "\""
print(state_str)


