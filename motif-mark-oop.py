#!/usr/bin/env python

'''
Author: Leyla Cufurovic
Python Version 3.12.2
PyCairo Version 1.26.0
'''

import argparse
import re
import cairo
import os
import random


from bioinfo import IUPAC_bases, oneline_fasta

def get_args():
    '''
    Initiate argparse. 
    Required inputs: fasta file, text file with motifs, and output directory
    '''
    parser = argparse.ArgumentParser()
    parser.add_argument("-f","--fasta_input", help ="Input fasta file in which to find motifs.", type = str, required = True)
    parser.add_argument("-m","--motif_input",help = "Input file containing motifs.",type = str, required = True)
    parser.add_argument("-o","--out_directory",help = "Path to directory to save output.",type = str, default="./")

    return parser.parse_args()
args = get_args()

fasta_file = args.fasta_input
motif_file = args.motif_input

# classes
class Motif:
    '''
    The Motif class identifies all motif patterns based on the given input motif text file using IUPAC bases.
    Patterns are stored in a list.
    '''
    def __init__(self, motif):
        self.motif = motif
        self.caps_motif = motif.upper()
        self.length = len(motif)
        self.pattern = self.find_pattern(self.caps_motif)
        
    def print_motif(self,motif:str) -> str:
        return motif

    def find_pattern(self,motif:str) -> list: #iterate over motifs to store a list of all the possible patterns using IUPAC
        ''' 
        input: UYRM 
        output: ['[TU][CTU][AG][AC]']
        '''
        pattern = [''.join([IUPAC_bases[base] for base in motif])]
        return pattern
    


class Gene:
    '''
The Gene class represents a gene from a FASTA file, extracting information such as gene name, chromosome location, exons, and motif matches. It provides methods to store motif information and find matches for motifs in the gene sequence.
    '''
    def __init__(self, header, sequence:str):
        self.gene = self.gene_name(header) #first line of seq
        self.chrom = self.extract_chrom(header) #chrom location
        self.sequence = sequence
        self.exons = self.extract_exons(sequence)
        self.matches = {}
        self.length = len(self.sequence)

    def gene_name(self,header:str) -> str:
        '''
        Splits the ID line using space as a delimiter.
        Takes the first part (before the first space).
        '''
        gene = header.split(' ')[0][1:]
        return gene

    def extract_chrom(self,header:str) -> str:
        chrom = header.split(' ')[1]
        return chrom

    def extract_exons(self,sequence:str) -> list:
        exon = [(match.start(0)+1,match.end(0)+1) for match in re.finditer(pattern='[A-Z]+',string=sequence)]
        return exon

    def store_motif(self,motif:Motif) -> None:
        self.matches[motif] = []
        return None
    
    def find_matches(self,motif:Motif,sequence:str) -> None:
        self.store_motif(motif)
        for patterns in motif.pattern:
            matches = []
            patterns = '(?={0})'.format(patterns)
            # plus one accounts for zero indexing
            pattern_matches = [m.start(0)+1 for m in re.finditer(patterns,sequence,re.IGNORECASE)]
            if len(pattern_matches) > 0:
                for pattern_match in pattern_matches:
                    self.matches[motif].append(pattern_match)
        return None

class DrawGene():
    '''
    The DrawGene class manages a collection of genes, keeps track of the longest gene length, and provides methods to add new genes, update the longest gene length, append genes to the list, and process a FASTA file to extract gene information and motif matches.
    '''
    def __init__(self):
        self.genes = []
        self.longest_gene = 0
        self.num_genes = 0

    def add_new_gene(self,gene:Gene) -> None:
        self.genes, self.num_genes = self.append_genes(gene,self.genes)
        self.longest_gene = self.new_gene_len(self.longest_gene,gene.length)

    def new_gene_len(self,longest_gene:int,gene_length:int) -> int:
        new_gene_len = max(longest_gene,gene_length)
        return new_gene_len

    def append_genes(self,gene:Gene,genes:list) -> tuple[list,int]:
        genes.append(gene)
        num_genes = len(genes)
        return genes, num_genes

    def make_oneline_fasta(self, fasta_file:str,motifs:list) -> None:
        open_fasta = open(fasta_file)
        while True:
            header = open_fasta.readline().rstrip()
            if header == '':
                break
            sequence = open_fasta.readline().rstrip()
            gene = Gene(header,sequence)
            self.add_new_gene(gene)
            for motif in motifs:
                motif_object = Motif(motif)
                gene.find_matches(motif_object, sequence)
        open_fasta.close()
        return None


class PyCairoDraw(cairo.Context):
    '''
    The PyCairoDraw class is responsible for drawing gene structures, exons, introns, and motif matches using the PyCairo library.
    '''
    def __init__(self, surface):
        self.motif_color_dict = {}
        self.surface = surface
        self.set_source_rgb(1, 1, 1)
        self.paint()


    def pyc_legend(self,x:float,y:int) -> None:
        # draw Exon key
        self.rectangle(
            x+580,
            y+10,
            20,
            20)
        self.set_source_rgb(0.08, 0.95, 0.95)
        self.set_line_width(6)
        self.fill()

        self.move_to(x + 607, y + 25)
        self.set_source_rgb(1,1,1)
        self.set_font_size(18)
        self.select_font_face('Arial', cairo.FONT_SLANT_NORMAL, cairo.FONT_WEIGHT_NORMAL)
        self.show_text('Exon')

        # draw gene key for Intron
        self.move_to(x + 570, y + 40)
        self.set_source_rgb(1,1,1)
        self.line_to(x + 604, y + 40)
        self.stroke()

        self.move_to(x + 607, y + 43)
        self.set_source_rgb(1,1,1)
        self.set_font_size(18)
        self.select_font_face('Arial', cairo.FONT_SLANT_NORMAL, cairo.FONT_WEIGHT_NORMAL)
        self.show_text('Intron')

        # Set the initial y-coordinate to the bottom of the "Intron" and "Exon" keys
        y = y + 60  # Adjust the value as needed for vertical spacing

        for k in self.motif_color_dict:
            self.set_line_width(6)
            self.set_source_rgb(
                self.motif_color_dict[k][0],
                self.motif_color_dict[k][1],
                self.motif_color_dict[k][2])

            # Move to the correct position under "Intron" and "Exon" keys
            self.move_to(x+570, y-3)
            self.line_to(x+604, y-3)
            self.stroke()

            self.move_to(x+607, y)
            self.set_source_rgb(1,1,1)
            self.set_font_size(18)
            self.select_font_face('Arial', cairo.FONT_SLANT_NORMAL, cairo.FONT_WEIGHT_NORMAL)
            self.show_text(k)

            # Increase the vertical spacing between motif keys
            y += 20  

    def motif_colors(self) -> tuple[float, float, float]:
        red = random.uniform(0, 1)
        green = min(random.uniform(0, 1) - 0.2, 1.0)
        blue = random.uniform(0, 1)
        self.set_source_rgb(red, green, blue)
        return red, green, blue

    
    def draw_motif(self,motif_object:Motif) -> None:
            if motif_object.motif not in self.motif_color_dict:
                self.motif_color_dict[motif_object.motif] = self.motif_colors()
            self.set_source_rgb(
                self.motif_color_dict[motif_object.motif][0],
                self.motif_color_dict[motif_object.motif][1],
                self.motif_color_dict[motif_object.motif][2]

            )
            self.set_line_width(3)
            return None

    def draw_gene(self, start: float, gene_object: Gene, start_next_gene: int) -> None:
        self.set_line_width(3)  # Thick black line
        self.set_source_rgb(1,1,1)  # Black color
        self.move_to(start, start_next_gene)
        self.line_to(gene_object.length + start, start_next_gene)
        self.stroke()
        return None

    def draw_exons(self, gene_object: Gene, start: float, init_height: int) -> None:
        for exon in gene_object.exons:
            self.set_source_rgb(0.08, 0.95, 0.95)
            self.set_line_width(4)
            self.rectangle(
                exon[0] + start,
                init_height+30,
                exon[1] - exon[0],
                10)
            self.fill()
        return None

    def current_gene_title(self,surface_width:float,height:int,gene:str,chr_location:str) -> None:
        self.set_source_rgb(1,1,1)
        self.set_font_size(14)
        self.select_font_face('Arial',cairo.FONT_SLANT_NORMAL,cairo.FONT_WEIGHT_NORMAL)
        self.move_to(surface_width*.03,height-7.7)
        self.show_text(gene)

        self.set_source_rgb(1,1,1)
        self.set_font_size(12)
        self.select_font_face('Arial',cairo.FONT_SLANT_NORMAL,cairo.FONT_WEIGHT_NORMAL)
        self.move_to(surface_width*.063,height-7.6)
        self.show_text(chr_location)

    def draw_pc_image(
        self,
        gene_collection:DrawGene,
        output_file:str,
        surface_width:float,
        surface_height:int,
        init_height:int=100) -> None:
        start = surface_width * .03
        
        #background color
        self.save()
        self.set_source_rgb(0.3,0.3,0.3)
        self.paint()
        self.restore()

        for gene_object in gene_collection.genes:
            motif_position_dict = {}
            start_next_gene = init_height + 35

            self.current_gene_title(surface_width,init_height,gene_object.gene,gene_object.chrom)
            self.draw_gene(start,gene_object,start_next_gene)
            self.draw_exons(gene_object,start,init_height)

            for motif_object in gene_object.matches:

                self.draw_motif(motif_object)

                for motif_matches in gene_object.matches[motif_object]:
                    adjust_motif_height = 4
                    while True:
                        motif_overlap = set(range(motif_matches,motif_matches+len(motif_object.motif)+ 8)) 
                        if adjust_motif_height in motif_position_dict.keys():
                            if motif_overlap.isdisjoint(motif_position_dict[adjust_motif_height]):
                                motif_position_dict[adjust_motif_height].update(motif_overlap)   
                                break                          
                            else:
                                adjust_motif_height += 4
                        else:
                            motif_position_dict[adjust_motif_height] = motif_overlap
                            break
                        
                    self.move_to(motif_matches+start,start_next_gene-adjust_motif_height)
                    self.line_to(motif_matches+start+len(motif_object.motif),start_next_gene-adjust_motif_height)
                    self.stroke()

            # adjust height for next gene
            init_height += 100

        # generate legend
        self.pyc_legend(x=surface_width*0.1,y=50)
        # save as png
        self.surface.write_to_png(output_file)
        self.surface.finish()
        return None
    
def motif_list(motif_file:str) -> list:
    f = open(motif_file,'r')
    motifs = []
    while True:
        motif = f.readline().rstrip()
        if motif == '':
            break
        motifs.append(motif)
    f.close()
    return motifs

def define_dims(gene_collection:DrawGene) -> tuple[float,int]:
    height = 100
    surface_width = gene_collection.longest_gene + gene_collection.longest_gene * .2
    surface_height = len(gene_collection.genes) * 100 + height

    return surface_width, surface_height

def png_outfile(fasta_input: str, out_directory) -> str:
    png_ext = '.png'
    # out_directory = "/Users/leycufur/bioinfo/Bi625"  # Specify the desired output directory
    prefix = os.path.splitext(os.path.basename(fasta_input))[0]
    png_output_file = os.path.join(out_directory, f"{prefix}{png_ext}")
    return png_output_file


def main():
    random.seed(77)
    args = get_args()
    
    oneline_file = oneline_fasta(args.fasta_input)
    output_file = png_outfile(args.fasta_input,args.out_directory)
    
    output_directory = os.path.dirname(output_file)
    os.makedirs(output_directory, exist_ok=True)

    motifs = motif_list(args.motif_input)

    gene_collection = DrawGene()
    gene_collection.make_oneline_fasta(oneline_file,motifs)

    # Adjust the scaling factor to make the image bigger
    scale_factor = 2 
    surface_width, surface_height = define_dims(gene_collection)
    surface = cairo.SVGSurface('plot.svg', surface_width * scale_factor, surface_height * scale_factor)

    context = PyCairoDraw(surface)
    context.scale(scale_factor, scale_factor)  # Scale the entire content

    context.draw_pc_image(
        gene_collection,
        output_file=output_file,
        surface_width=surface_width*scale_factor,
        surface_height=surface_height*scale_factor)

    # remove temporary fasta and plot
    os.remove('plot.svg')
    os.remove(f'singleline_{fasta_file}')
    

if __name__ == '__main__':
    main()
