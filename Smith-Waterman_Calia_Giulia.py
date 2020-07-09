import argparse
import sys
import os 

MATCH = +3
MISMATCH = -3
GAP = -2

def params():
    obj_parse_args = argparse.ArgumentParser(description="""Given two sequences
        or a file containing sequences return the alignments.
        The sequences in the input file have to be tab separated""")
    obj_parse_args.add_argument("--seq1",
                                type = str,
                                help = "The first sequence to align")
    obj_parse_args.add_argument("--seq2",
                                type = str,
                                help = "The second sequence to align")
    obj_parse_args.add_argument("-i", "--input",
                                type = str,
                                help = "The file given in input")
    obj_parse_args.add_argument("-m", "--match",
                                default = MATCH,
                                type = float,
                                required = False,
                                help = "Score of the match (default = {})".format(MATCH))
    obj_parse_args.add_argument("-s", "--mismatch",
                                default = MISMATCH,
                                type = float,
                                required = False,
                                help = "Score of the mismatch (default = {})".format(MISMATCH))
    obj_parse_args.add_argument("-g", "--gap",
                                default = GAP,
                                type = float,
                                required = False,
                                help = "Score of the gap (default = {})".format(GAP))
    obj_parse_args.add_argument("--minscore",
                                type = float,
                                help = "The minimum score")
    obj_parse_args.add_argument("--minlength",
                                type = float,
                                help = "The minimum length")
    obj_parse_args.add_argument("--numresult",
                                type = int,
                                help = "The number of alignments to return")

    args = obj_parse_args.parse_args()
    return args                              
                           
def check_params(parameters):
    parameters.seqs = []
    input_f = False
    if parameters.input is not None:
        if os.path.exists(parameters.input):
            if os.path.isfile(parameters.input):
                size = len(parameters.seqs)
                parameters.seqs += [tuple(line.strip().split("\t")) for line
                    in open(parameters.input) if len(line.strip().split("\t"))==2]
                input_f = True
                if len(parameters.seqs) <= size:
                    print("""the file could be empty or wrongly formatted: \nseq1\tseq2
                        \nseq3\tseq4\n..\n(tab separated)""")
            else:
                print("input is not a file")
        else:
            print("input file doesn't exist!")

    if ((parameters.seq1 is not None) and (parameters.seq2 is not None)):
        parameters.seqs.append((parameters.seq1, parameters.seq2))

    if ((parameters.seq1 is None) or (parameters.seq2 is None)) and input_f is False:
        print("You have to give in input at least two sequences or a file of sequences!")
        if not parameters.seqs:
            print("You did not give me anything...")
            sys.exit(1)

    if parameters.match <= 0:
        print("provided value for matches is not correctm default will be used instead")
        parameters.match = MATCH

    if parameters.mismatch > 0:
        print("provided value for mismatches is not correct, default will be used instead!")
        parameters.mismatch = MISMATCH

    if parameters.gap > 0:
        print("provided value for gap is not correct, default will be used instead")
        parameters.gap = GAP

    if parameters.minscore is not None:
        if parameters.minscore < 0:
            print("provided value for minimum score is negative, no minimum score will be used")
            parameters.minscore = None

    if parameters.minlength is not None:
        if parameters.minlength < 0:
            print("provided value for minimum length is negative, no minimum length will be used")
            parameters.minlength = None

    return parameters

def filling_table(seqs, m=MATCH, mm=MISMATCH, g=GAP):
    """Fills a score matrix taking in input values/penalties for 
        match, mismatch and gap. Moreover it fills a direction matrix
        that take in account the direction of the chosen score"""
    lseq1 = len(seqs[0])+1
    lseq2 = len(seqs[1])+1
    s1 = seqs[0]
    s2 = seqs[1]
    
    score_matrix = [[0 for i in range(lseq2)] for j in range(lseq1)]
    direction_matrix = [[[0,0,0] for i in range(lseq2)] for j in range(lseq1)]

    for i in range(1, lseq1):
        for j in range(1, lseq2):
            diag_score_m = score_matrix[i-1][j-1] + m
            diag_score_mm = score_matrix[i-1][j-1] + mm
            diag_score = diag_score_m if s1[i-1] == s2[j-1] else diag_score_mm

            left_score = score_matrix[i][j-1] + g
            up_score = score_matrix[i-1][j] + g

            score_matrix[i][j] = max(diag_score, left_score, up_score, 0)
            if score_matrix[i][j] == diag_score:
                direction_matrix[i][j][0] = 1
            elif score_matrix[i][j] == left_score:
                direction_matrix[i][j][1] = 1
            elif score_matrix[i][j] == up_score:
                direction_matrix[i][j][2] = 1
            else:
                direction_matrix[i][j] = [0,0,0]
    
    return score_matrix, direction_matrix

def connecting_tables(seqs, matrix):
    """Provides a dictionary with the score and the relative position(s)
        in the matrices in addiction to the maximum score of the matrix
        and separately its position"""
    lseq1 = len(seqs[0]) + 1
    lseq2 = len(seqs[1]) + 1
    score_dict = {}

    for i in range(lseq1):
        for j in range(lseq2):
            if matrix[i][j] not in score_dict:
                score_dict[matrix[i][j]] = [[i,j]]
            else:
                score_dict[matrix[i][j]].append([i,j])

    scores = sorted(score_dict.keys(), reverse = True)
    score_max = max(scores)
    pos_max = score_dict[score_max]

    return score_dict, scores, score_max, pos_max

def traceback(seqs, matrix, direction_matrix, pos_max):
    """Allows to build up the local alignment proceeding backwards till 
        a 0 score starting from the maximum score of the matrix"""
    s1 = seqs[0]
    s2 = seqs[1]
    seq1 = ""
    seq2 = ""
    mid_seq = ""

    for start_pos in pos_max:
        score = matrix[start_pos[0]][start_pos[1]]

        while score != 0:
            up_pos = [start_pos[0]-1, start_pos[1]]
            left_pos = [start_pos[0], start_pos[1]-1]
            diag_pos = [start_pos[0]-1, start_pos[1]-1]

            up_score = matrix[up_pos[0]][up_pos[1]]
            left_score = matrix[left_pos[0]][left_pos[1]]
            diag_score = matrix[diag_pos[0]][diag_pos[1]]

            current_direction = direction_matrix[start_pos[0]][start_pos[1]]

            if current_direction[0] == 1: #== [1, 0, 0]:
                score = diag_score
                start_pos = diag_pos
                seq1 = s1[diag_pos[0]] + seq1 
                seq2 = s2[diag_pos[1]] + seq2
                if s1[diag_pos[0]] == s2[diag_pos[1]]:
                    mid_seq = "|" + mid_seq
                else:
                    mid_seq = ":" + mid_seq            

            elif current_direction[1] == 1: #== [0, 1, 0]:
                score = left_score
                start_pos = left_pos
                seq1 = "-" + seq1
                seq2 = s2[left_pos[1]] + seq2
                mid_seq = " " + mid_seq

            elif current_direction[2] == 1: #== [0, 0, 1]:
                score = up_score
                start_pos = up_pos
                seq1 = s1[up_pos[0]] + seq1
                seq2 = "-" + seq2
                mid_seq = " " + mid_seq
    
        return seq1, mid_seq, seq2


if __name__ == "__main__":
    args = params()
    args = check_params(args)
    for i in args.seqs:
        scoreMatrix, directionMatrix = filling_table(i, MATCH, MISMATCH, GAP)
        connecting = connecting_tables(i, scoreMatrix)
        scoreDict = connecting[0]
        scores = connecting[1]
        scoreMax = connecting[2]
        posMax = connecting[3]
    
        if args.numresult == None:
            if args.minscore == None or args.minscore <= scoreMax:
                seq1, mid, seq2 = traceback(i, scoreMatrix, directionMatrix, posMax)
                if args.minlength == None or args.minlength <= len(seq1):
                    print("MAX SCORE: {}\n{}\n{}\n{}\n".format(scoreMax, seq1, mid, seq2))
                else: 
                    print("ATTENTION!\nThe length of the alignment is too short: {} (minlength: {})".format(len(seq1), args.minlength))
            else:
                print("ATTENTION!\nThe minimum score is too small: {} (minscore: {})".format(scoreMax, args.minscore))
        else:
            count = args.numresult
            result_positions = []
            n = 0 
            j = 0 

            while count != 0:
                result_positions.append([scoreDict[scores[n]][j]])
                count -= 1
                if len(scoreDict[scores[n]]) > (j+1):
                    j += 1
                else:
                    if len(scoreDict) > (n+1):
                        n += 1
                        j = 0

            for position in result_positions:
                scoreMax = scoreMatrix[position[0][0]][position[0][1]]
                
                if args.minscore == None or args.minscore <= scoreMax:
                    seq1, mid, seq2 = traceback(i, scoreMatrix, directionMatrix, position)
                    
                    if args.minlength == None or args.minlength <= len(seq1):
                        print("ACTUAL MAX SCORE: {}\n{}\n{}\n{}\n".format(scoreMax, seq1, mid, seq2))
                    else:
                        print("ATTENTION!\nThe length of the alignment is too short: {} (minlength: {})".format(len(seq1), args.minlength))
                else:
                    print("ATTENTION!\nThe minimum score is too small: {} (minscore: {})".format(scoreMax, args.minscore))