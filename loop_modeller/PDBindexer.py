class PDBindexer:
    def __init__(self, filename):
        self.fn = filename
        self.resis = None
        self.parse()

    def parse(self):
        with open(self.fn) as f:
            lines = f.readlines()
        self.resis = []
        n = len(lines)
        resi = []
        for i in range(n):
            l = lines[i]
            if (l[0:4] == 'ATOM' or l[0:4] == 'HETA'):
                if len(resi) == 0:
                    resi.append(l)
                else:
                    resnum = int(l[22:26])
                    curresnum = int(resi[0][22:26])
                    if resnum == curresnum:
                        resi.append(l)
                    else:
                        self.resis.append(resi)
                        resi = [l]
        if len(resi) > 0:
            self.resis.append(resi)

    def renumber(self):

        def renumber_helper(l, resnum, atomi):
            return (l[0:6] +
                    "{0: 5d}".format(atomi + 1) +
                    l[11:22] +
                    "{0: 4d}".format(resnum + 1) +
                    l[26:])

        atomcount = 0
        for i in range(len(self.resis)):
            resi = self.resis[i]
            for j in range(len(resi)):
                resi[j] = renumber_helper(resi[j], i, atomcount)
                atomcount += 1

    def write(self, filename):
        with open(filename, "w") as f:
            for resi in self.resis:
                for line in resi:
                    f.write(line)
