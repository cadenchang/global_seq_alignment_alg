import sys

class GlobalSeqAlignment:
    def __init__(self, wordx, wordy, matchScore, mismatchScore, gapScore):
        self.x = wordx
        self.y = wordy
        self.M = matchScore
        self.m = mismatchScore
        self.g = gapScore
        #initializing score matrix and pointer matrix
        rows = len(self.x) + 1
        cols = len(self.y) + 1
        self.matrix = [[0 for i in range(cols)] for j in range(rows)]
        self.pointers = [['x' for i in range(cols)] for j in range(rows)]

    # member function that runs the main algorithm (populates score and 
    # pointer matrix)
    def get_optimal_alignment(self):
        for i in range(len(self.matrix)):
            for j in range(len(self.matrix[0])):
                #calculate edge values
                if i == 0 and j == 0:
                    self.matrix[i][j] = 0
                elif i == 0 and j >= 1:
                    self.matrix[i][j] = self.matrix[i][j-1] + self.g
                    self.pointers[i][j] = 'l' 
                elif j == 0 and i >= 1:
                    self.matrix[i][j] = self.matrix[i-1][j] + self.g
                    self.pointers[i][j] = 'u'
                else:
                    xi = self.x[i-1]
                    yj = self.y[j-1]
                    #calculate center values
                    if xi == yj:
                        self.matrix[i][j] = self.matrix[i-1][j-1] + self.M
                        self.pointers[i][j] = 'd' 
                    else:
                        mismatch = self.matrix[i-1][j-1] + self.m
                        upgap = self.matrix[i-1][j] + self.g
                        leftgap = self.matrix[i][j-1] + self.g
                        if (mismatch >= upgap and mismatch >= leftgap):
                            self.pointers[i][j] = 'd'
                            self.matrix[i][j] = mismatch
                        elif (upgap >= mismatch and upgap >= leftgap):
                            self.pointers[i][j] = 'u'
                            self.matrix[i][j] = upgap
                        else:
                            self.pointers[i][j] = 'l'
                            self.matrix[i][j] = leftgap
        # return sequence alignment solution                    
        return self.trace_back()
        
    # assuming a populated score and pointer matrix, returns the two aligned
    # sequences
    def trace_back(self):
        xtracker = len(self.matrix) - 1
        ytracker = len(self.matrix[0]) - 1
        x_prime = ""
        y_prime = ""
        while xtracker > 0 or ytracker > 0:
            direction = self.pointers[xtracker][ytracker] 
            if direction == 'd':
                x_prime = self.x[xtracker-1] + x_prime
                y_prime = self.y[ytracker-1] + y_prime
                xtracker = xtracker - 1
                ytracker = ytracker - 1
            elif direction == 'u':
                y_prime = "-" + y_prime
                x_prime = self.x[xtracker-1] + x_prime
                xtracker = xtracker - 1
            else:
                x_prime = "-" + x_prime
                y_prime = self.y[ytracker-1] + y_prime
                ytracker = ytracker - 1

        return x_prime, y_prime

# prints a 2d array for debugging purposes
def print_matrix(arr):
        for i in range(len(arr)):
            for j in range(len(arr[0])):
                print(arr[i][j], end=" ")
            print("\n")

# given an input stream, reads in a file in FASTA format and returns the two
# starting sequences
def parse_FASTA_file(input):
    sequences = []
    current = ""
    for line in input:
        line = line.strip()
        if line.startswith(">"):
            if current != "":
                sequences.append(current)
            current = ""
        else:
            current += line
    sequences.append(current)
    return tuple(sequences)

 
if __name__ == "__main__":
    x, y = parse_FASTA_file(sys.stdin)
    align = GlobalSeqAlignment(x, y, 4, -2, -2)
    x_prime, y_prime = align.get_optimal_alignment()
    print(x_prime)
    print(y_prime)