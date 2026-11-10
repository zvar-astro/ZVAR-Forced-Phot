import os

filename = "rsync.log"

messages = []
fields = []
ccds = []
quads = []
filters = []
with open(filename) as file:

    for line in file.readlines():
        if "WARNING" in line:
            msg = 'WARNING'
            quadID = line.strip("\n").split(" ")[4]
        elif "ERROR" in line:
            msg = 'ERROR'
            quadID = line.strip("\n").split(" ")[4]
        else:
            continue
        
        fields.append(int(quadID.split("_")[1]))
        ccds.append(int(quadID.split("_")[2]))
        quads.append(int(quadID.split("_")[3]))
        filters.append(quadID.split("_")[4][1])
        messages.append(msg)

Nfailed = len(fields)

savename = "./rsync_warnings.txt"
if os.path.isfile(savename):
    os.remove(savename)

for i in range(Nfailed):
    line_text = f"{fields[i]} {ccds[i]} {quads[i]} {filters[i]}"
    print(messages[i],line_text)
    with open(savename,'a') as file:
        file.write(f"{line_text}\n")
