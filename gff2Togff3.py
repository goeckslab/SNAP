import argparse
import sys
import fileinput
from Group import Group

def main():
        parser = argparse.ArgumentParser(description='Get a gff file and the output gff3 file')
        parser.add_argument('--input', help='input gff file')
        parser.add_argument('--output', help='output gff3 file', required=True)
        args = parser.parse_args()
        input = args.input
        output = args.output
        if not sys.stdin.isatty():
            c = Convertor(sys.stdin, output)
        else: 
            c = Convertor(input, output)
        c.convert()
        
class Convertor:    
    def __init__(self, input, output):
        if type(input) is str:
            with open(input) as self.f:
                self.li = [line.rstrip().split("\t") for line in self.f]
        else:
            self.li = [line.rstrip().split("\t") for line in input]
        self.gff3 = open(output, "w")
        self.gff3.write("##gff-version 3\n")

    def convert(self):
        index = 0
        while index in range(0, len(self.li)):
            index = self.groupAsgene(index)
        self.gff3.close()
                
                    
    def groupAsgene(self, start = 0):
        gene = self.li[start][8]
        index = len(self.li)
        for i in range(start+1, len(self.li)):
            line = self.li[i]
            if gene != line[8]:
                index = i
                break
        if index >= len(self.li):
            group = self.li[start:len(self.li)]
        else:
            group = self.li[start:index]
        g = Group(group)
        g.writer(self.gff3)
        return index

   
        

if __name__ == "__main__":
    main()


    