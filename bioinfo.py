# Author: Leyla Cufurovic leylacuf@uoregon.edu

# Check out some Python module resources:
#   - https://docs.python.org/3/tutorial/modules.html
#   - https://python101.pythonlibrary.org/chapter36_creating_modules_and_packages.html
#   - and many more: https://www.google.com/search?q=how+to+write+a+python+module

'''This module is a collection of useful bioinformatics functions
written during the Bioinformatics and Genomics Program coursework.
You should update this docstring to reflect what you would like it to say'''

__version__ = "Python 3.12.2"         

DNA_bases = set("ACTGNactgn")
RNA_bases = set("AUGCNaugcn")
IUPAC_bases = {
    "A":"[Aa]",
    "C":"[Cc]",
    "G":"[Gg]",
    "T":"[Tt]",
    "U":"[UuTt]",
    "W":"[AaTtUu]",
    "S":"[CcGg]",
    "M":"[AaCc]",
    "K":"[GgTtUu]",
    "R":"[AaGg]",
    "Y":"[CcTtUu]",
    "B":"[CcGgTt]",
    "D":"[AaGgTtUu]",
    "H":"[AaCcTtUu]",
    "V":"[AaCcGg]",
    "N":"[AaCcGgTtUu]"
}


def init_list(lst: list, value: float=0.0) -> list:
    '''This function takes an empty list and will populate it with
    the value passed in "value". If no value is passed, initializes list
    with 101 values of 0.0.'''
    i = 0
    for x in range(101):
        lst.append(value)
    return lst

def populate_list(file: str) -> tuple[list, int]:
    '''Creates and empty list using init_list(), opens Fastq file, loops through every individual sequence record
    and converts the phred score to a number, then adds the qual scores to an cumulative sum in the list, keeps track of
    total number of lines in the file, and returns the list and number of lines as a tuple.'''
    num_lines = 0 #Initialize line counter
    my_list: list = []
    my_list = init_list(my_list)
        
    with open(file, "r") as fastq: #open fastq file
        for line in fastq: #iterate over each line in the file
            num_lines += 1 #increment the line counter
            
            if num_lines%4==0: #use only the 4th line of the record that contains quality scores
                phred_scores = line.strip() #strip whitespaces,  new lines, etc.
                for i, phred in enumerate(phred_scores): #iterate over each quality score character
                    my_list[i] += bioinfo.convert_phred(phred) #convert quality scores to phred scores and add the quality scores to my_list
                    
    return my_list, num_lines # return list of quality scores (my_list) and line count (num_lines)

def convert_phred(letter: str) -> int:
    '''Converts a single character into a phred score'''
    return ord(letter)-33

def qual_score(phred_score: str) -> float:
    """This function calculates the quality score by taking the unmodified phred score string as a parameter and then calculates the average quality score of the entire string. ****REQUIRES convert_phred****"""
    total_score = 0
    count = 0
    for letter in phred_score:
        indv_score = convert_phred(letter)
        total_score += indv_score
        count += 1
    avg_score = total_score / count if count > 0 else 0
    return avg_score 

def calc_median(lst):
    """Calculates and returns the median of a sorted one-dimensional list."""
    length = len(lst)
    if length % 2 == 0:
        # If the list has an even length, calculate the average of the middle two elements
        mid1 = length // 2 - 1
        mid2 = length // 2
        median = (lst[mid1] + lst[mid2]) / 2
    else:
        # If the list has an odd length, the median is the middle element
        mid = length // 2
        median = lst[mid]
    return median


def validate_base_seq(seq, RNAflag=False):
    '''This function takes a string. Returns True if string is composed
    of only As, Ts (or Us if RNAflag), Gs, Cs. False otherwise. Case insensitive.'''
    return len(seq) == seq.count("A") + seq.count("U" if RNAflag else "T") + seq.count("G") + seq.count("C")

def gc_content(seq):
    '''Returns GC content of a DNA or RNA sequence as a decimal between 0 and 1.'''
    seq = seq.upper()
    return (seq.count('G') + seq.count('C'))/ len(seq)

def oneline_fasta(input_fasta):
    output_filename = "singleline_" + input_fasta

    with open(input_fasta, "r") as fh1, open(output_filename, 'w') as fh2:
        header = fh1.readline().strip()
        sequence = ""
        for line in fh1:
            if line.startswith(">"):
                fh2.write(header + "\n")
                fh2.write(sequence + "\n")
                header = line.strip()
                sequence = ""
            else:
                sequence += line.strip()
        fh2.write(header + "\n")
        fh2.write(sequence + "\n")
    return output_filename

def rev_complement(sequence):
    '''This function takes a DNA nucleotide sequence and returns the reverse complement'''
    complement = {'A':'T', 'C':'G', 'T':'A', 'G':'C', 'N':'N'}
    return ''.join(complement[base] for base in sequence[::-1])


    

if __name__ == "__main__":
    # write tests for functions above, Leslie has already populated some tests for convert_phred
    assert convert_phred("I") == 40, "wrong phred score for 'I'"
    assert convert_phred("C") == 34, "wrong phred score for 'C'"
    assert convert_phred("2") == 17, "wrong phred score for '2'"
    assert convert_phred("@") == 31, "wrong phred score for '@'"
    assert convert_phred("$") == 3, "wrong phred score for '$'"
    print("Your convert_phred function is working! Nice job")

    # test for gc_content
    assert gc_content('CGAT') == 0.5, "try again..."
    assert gc_content('CGAT') != 1, "I would love to stay and chitchat, but I gotta go (try again)"
    assert gc_content('CCGCGGTCGCCCACCCCTGATATATAGGC') == 0.6551724137931034, "pls why is this not working"
    print("wow amazing, wonderful, good job <3")
    
    # tests for qual_score
    assert qual_score("E!J-E")== 25,  "nope"
    assert qual_score("JKL") == 42.0,  "absolutely not"
    print("✯¸.•´*¨`*•✿ 𝓌𝒽𝒶𝓉 𝒶 𝓂𝒶𝑔𝓃𝒾𝒻𝒾𝒸𝑒𝓃𝓉 𝓅𝑒𝓇𝓈𝑜𝓃 𝓎𝑜𝓊 𝒶𝓇𝑒, 𝑔𝑜𝑜𝒹 𝒿𝑜𝒷 ✿•*`¨*`•.¸✯")
    
    # tests for validate_base_seq
    assert validate_base_seq("Eeeeeep!!!") == False, "Fails to recognize nonDNA"
    assert validate_base_seq("UUUUAAACCG", True) == True, "Validate base seq does not work on RNA"
    assert validate_base_seq("AGTCG") == True, "Fails to recognize nonDNA"
    assert validate_base_seq("Oh...", True) == False, "Validate base seq fails to recognize nonDNA"
    print("☄. *. ⋆Passed DNA and RNA tests☄. *. ⋆")
    
    # tests for calc_median
    assert calc_median([10, 20, 30]) == 20, "function not working :C"
    assert calc_median([2, 4, 6, 8]) == 5, "function not working :C"
    assert calc_median([798, 58, 12589]) == 58, "function not working :C"
    print("*ੈ✩‧₊˚Yay, you did math!*ੈ✩‧₊˚")
    
    # tests for populate_list
    assert type(populate_list(file)) == tuple
    assert type(populate_list(file)[0]) == list
    assert type(populate_list(file)[1]) == int
    assert len(populate_list(file)[0]) == 101
    print("(っ◔◡◔)っ ♥ you get a ✩𝕤𝕥𝕒𝕣✩ ♥")
    
    