from operator import itemgetter

# Input: A group: a list that contains lines belonging to the same gene
class Group:
    # Modify "type" column and "attributes" colunm, initialize id, gene, source, stream
    def __init__(self, group):
        self.group = group
        self.id = str(group[0][0])
        self.source = str(group[0][1])
        self.stream = str(group[0][6])
        self.gene = str(group[0][8])
        for x in range(0, len(group)):
            self.group[x][2] = "CDS"
            self.group[x][8] = "Parent=mRNA_" + self.gene
            self.group[x][3] = int(self.group[x][3])
            self.group[x][4] = int(self.group[x][4])

    # Order the group elements accoriding to Stream, +: ascanding order, -: descanding order   
    def order(self):
        self.num = len(self.group)
        if self.stream == "+":
            self.group = sorted(self.group, key=itemgetter(3))
            self.min_item = self.group[0][3]
            self.max_item = self.group[self.num-1][4]
        elif self.stream == "-":
            self.group = sorted(self.group, key=itemgetter(3), reverse=True)
            self.min_item = self.group[self.num-1][3]
            self.max_item = self.group[0][4]
        else:
            print("Stream in invalid!\n")
    
    def phaseCalculator(self, i, donor = 0):
        if i >= self.num:
            pass
        else:
            self.type = self.group[i][2]
            self.size = self.group[i][4] - self.group[i][3] + 1
        if self.num == 1:
            if self.type == "Eterm":
                self.group[i][7] = str(self.size % 3)
            else:
                self.group[i][7] = "0"
        elif self.num > 1 and i < self.num:
            accept = (3 - donor) % 3
            self.group[i][7] = str(accept)
            donor = (self.size - accept) % 3
            i = i + 1
            self.phaseCalculator(i, donor)
            
    
    def writer(self, gff3):
        self.order()
        self.phaseCalculator(0)
        gff3.write(self.id + "\t" + self.source + "\tgene\t" + str(self.min_item) + "\t" + str(self.max_item) + "\t.\t" + self.stream + "\t.\t" + "ID=" + self.gene + "\n")
        gff3.write(self.id + "\t" + self.source + "\tmRNA\t" + str(self.min_item) + "\t" + str(self.max_item) + "\t.\t" + self.stream + "\t.\t" + "ID=mRNA_" + self.gene + ";Parent=" + self.gene + "\n")
        for x in range(0, len(self.group)):
            self.group[x][3] = str(self.group[x][3])
            self.group[x][4] = str(self.group[x][4])
            gff3.write("\t".join(self.group[x]) + "\n")