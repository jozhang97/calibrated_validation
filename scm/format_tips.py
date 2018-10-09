sp = " sp2  sp3  sp4  sp5  sp6  sp7  sp8  sp9 sp10 sp12 sp13 sp15 sp17 sp18 sp19 sp20 sp21 sp22 sp23 sp24 sp25 sp26"
sp_names = sp.replace("sp", "\",\"sp")
sp_names += "\""
sp_names = sp_names.replace(" ", "")

print(sp_names)


states = "   0    1    1    1    1    1    1    1    1    0    0    0    0    0    0    0     0    0    1    1    1    1"
states = states.split()
states = [int(s) + 1 for s in states]

sp = sp.split()

state_str = "\""
for i, j in zip(sp, states):
    state_str += i + "=" + str(j) + ","

state_str += "\""
print(state_str)


nd = "nd1  nd2  nd6 nd22  nd7  nd9 nd11  nd3  nd4  nd8 nd10 nd12 nd15 nd16 nd17 nd13 nd21  nd5 nd14 nd18 nd20"
print(nd.split())
states = "0    1    1    0    1    1    1    0    0    0    0    0    0    0    1    0    0    0    0    0    0"
states = states.split()
states = [int(s) + 1 for s in states]
print(states)
