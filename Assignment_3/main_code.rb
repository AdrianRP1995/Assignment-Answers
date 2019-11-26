
require 'net/http'
require 'bio'

#We import BioRuby and the module to use HTTP

#===========================================================================================================================
#FUNCTIONS
#===========================================================================================================================


def fetch_uri(str)
    # This function retrieves a response from a web
    #Input: a string with the URI

    add = URI(str)
    res = Net::HTTP.get_response(add)

    return res  # return the response object to be used later

end


#===========================================================================================================================


def load_file(file)
  # Function to read the genes from a file
  #Input: a string with the name of the file

  f = File.open(file, "r")
  genes = Array.new #We create an array to contain each line (gene ID) of the file

  f.each_line do |line|
    genes << line.delete("\n") # We remove \n from each line
  end

  f.close
  return genes.uniq # uniq eliminates every possible duplication we could have of each gene id

end


#===========================================================================================================================


def gene_info (geneid)
  # Function that access the ensebl database and creates a Bio:EMBL object of a given gene, and turns it into a Bioseq object
  #Input: a gene ID

  add = "http://www.ebi.ac.uk/Tools/dbfetch/dbfetch?db=ensemblgenomesgene&format=embl&id=#{geneid}"
  res = fetch_uri(add) #We retrieve the info from the gene entry in EBI

  return nil unless res #We check if we retrieved the information correctly
  record = res.body
  entry = Bio::EMBL.new(record) #We make the information intro a BIO::EMBL object

  bioseq = entry.to_biosequence #We use to_biosequence method to turn the information into a bioseq object

  return bioseq

end


#===========================================================================================================================


def find_seq_exon (exonid, target_seq_match, lenght, exon_pos, strand)
   # Function that determines our target is in a given exon, taking into acount if the exon is in the forward or reverse strand
   # Input: and exon ID, the target we are checking, lenght of the target, the exon positions where the target is, and the orientation of the strand

   exon_target = Hash.new
   # A hash containing matched pos in exon [beginning, end] and [exonid, strand]

   case strand #A case in which we check if we are working with a forward (+) or (-) strand

   when '+'

     target_seq_match.each do |match_init|

        match_end = match_init + $len_target - 1 #We iterate over each position of the target match

        if (match_init >= exon_pos[0].to_i) && (match_init <= exon_pos[1].to_i) && (match_end >= exon_pos[0].to_i) && (match_end <= exon_pos[1].to_i)
           # We check whether the target is inside the exon, checking each position
           exon_target[[match_init, match_end]] = [exonid, '+'] #We introduce the data in the hash
        end

     end


   when '-'

      target_seq_match.each do |match_init|

        match_end = match_init + $len_target - 1 #Again, we iterate over each position of the target match

        if (match_init >= exon_pos[0].to_i) && (match_init <= exon_pos[1].to_i) && (match_end >= exon_pos[0].to_i) && (match_end <= exon_pos[1].to_i)
          # We check whether the target is inside the exon, checking each position

           # We need to work with the position that would correspond to the foward strand, so we convert the position:
           m_end = lenght - match_end - 1
           m_init = lenght - match_init - 1
          #When we reverse a reverse sequence, the position moves 1 position forward.
          #To correct this, we substract 1 from the positions.

           exon_target[[m_end, m_init]] = [exonid, '-'] #We introduce the converted data in the hash

        end

     end


   end

   if not exon_target.empty? # If there isn't a target in our exon, the function returns nothing
      return exon_target
   end

end


#===========================================================================================================================


def get_exon_tar (seq_obj)
  # Function that given the Bio::Sequence object returns a hash in which the keys are the coordinates of the target's matches inside exons.
  #Input: a Bio::Sequence object which contains our sequence

  len_seq = seq_obj.length() # Length of the sequence
  target_pos_exon = Hash.new # Hash that will contain the pos targeted inside exons as keys and the strand as values

  # We get the target's matches in both foward and reverse strand
  match_seq_for = seq_obj.gsub(/#{$target}/).map{Regexp.last_match.begin(0)}
  match_seq_rev = seq_obj.reverse.tr('atgc', 'tacg').gsub(/#{$target}/).map{Regexp.last_match.begin(0)}

  seq_obj.features.each do |feature|

    pos_for = feature.position

    next unless (feature.feature == 'exon' && (not pos_for =~ /[A-Z]/))
    # We look for the exon type feature and we ommit trans-splicing for simplicity purposes

    exonid = feature.qualifiers[0].value.gsub('exonid=', '') # We format the string

    if pos_for =~ /complement/ #We have to reverse the "-" strand to find the targets in it

        pos_for = pos_for.tr('complement()', '').split('..') #Getting a 2 elements array containg initial and end position
        pos_rev = [] #We convert it to the reverse strand

        pos_for.each do |pos|
          pos_rev.insert(0, len_seq - pos.to_i) # We use insert to give the correct order of the coordinates
        end
        target_exon_response = find_seq_exon(exonid, match_seq_rev, len_seq, pos_rev, '-')
        # We call "find_seq_exon" get the matches that are inside of the exon.
        # We send the function the matches and positions of the exon in the reverse strand

        if not target_exon_response.nil? #If a response is retrieved, we add the targets to the hash. If not, we do nothing
         target_pos_exon = target_pos_exon.merge(target_exon_response)
        end


    else #We look for targets in the "+" strand

        pos_for = pos_for.split('..') # Getting a 2 elements array containg initial and end position

        target_exon_response= find_seq_exon(exonid, match_seq_for, len_seq, pos_for, '+')
        # We call "find_seq_exon" get the matches that are inside of the exon.
        # We send the function the matches and positions of the exon in the forward strand

        if not target_exon_response.nil? #If a response is retrieved, we add the targets to the hash. If not, we do nothing

         target_pos_exon = target_pos_exon.merge(target_exon_response)
        end

    end


  end

  return target_pos_exon #We return the hash with all the exon positions

end


#===========================================================================================================================


def feats(geneid, targets, bioseq)
  #Function that iterates over the hash with the target in exons and add them as features to the Bio:EMBL object.

  features = Array.new

  targets.each do |target, ex_strand|

     feat = Bio::Feature.new("#{$target.upcase}_in_exon", "#{target[0]}..#{target[1]}")

     feat.append(Bio::Feature::Qualifier.new('nucleotide_motif', "#{$target.upcase}_in_#{ex_strand[0]}"))
     #Feature qualifier needed to writer a recognizable GFF3 file     # This format will be needed for the GFF3

     feat.append(Bio::Feature::Qualifier.new('strand', ex_strand[1]))

     $gff_genes.puts "#{geneid}\t.\t#{feat.feature}\t#{target[0]}\t#{target[1]}\t.\t#{ex_strand[1]}\t.\tID=#{ex_strand[0]}"
     #We print the feature in the GFF3 file

     features << feat
  end

  bioseq.features.concat(features) #We add the new features created to the previous ones

end


#===========================================================================================================================


def get_chr (geneid, seq_obj)
  # Function that takes a Bio:Sequence object and returns its chromosome and position
  #Input: a gene ID, and its sequence object

  bs_pa = seq_obj.primary_accession

  return false unless bs_pa

  chr_array = bs_pa.split(":")

  $gff_chr.puts "#{chr_array[2]}\t.\tgene\t#{chr_array[3]}\t#{chr_array[4]}\t.\t+\t.\tID=#{geneid}"
  # This is the line that will be prnted to the GFF

  # We return:
  #   - Chromosome number ---> [2]
  #   - Chromosome gene start position ---> [3]
  #   - Chromosome gene end position ---> [4]
  return chr_array[2], chr_array[3], chr_array[4]

end


#===========================================================================================================================


def create_file(file)
  #Function that creates a file
  #Input: file name

  if File.exists?(file)
    File.delete(file) # We remove the file in case it exits to update it
  end

  return File.open(file, "a+")

end


#===========================================================================================================================


def conv_to_chr(gene, targets, chr)
  # Given the gene ID, the hash containing the targets, and the information
  # about the chromosome, this method translates the coordinates to the ones
  # refering to the chromosome. It prints them on the GFF3 chromosome file


  targets.each do |pos, ex_strand|
    initial_chr = chr[1].to_i + pos[0].to_i
    end_chr = chr[1].to_i + pos[1].to_i

    $gff_chr.puts "#{chr[0]}\t.\tnucleotide_motif\t#{initial_chr}\t#{end_chr}\t.\t#{ex_strand[1]}\t.\tID=#{ex_strand[0]};parent=#{gene}"
  end


end


#===========================================================================================================================
#END OF FUNCTIONS
#===========================================================================================================================




#===========================================================================================================================

#MAIN PROGRAM

#===========================================================================================================================


$target = "cttctt" #This global variable contains the sequence we will search in the exons
$len_target = $target.length() #This global variable contains the lenght of the target

puts "Starting. The answer should be ready in about 5 minutes!\n"

#We create global variables with the creating files objects:
$gff_genes = create_file("genes_with_target.gff3") #GFF3 with gene osition
$gff_chr = create_file("chromosomes_with_target.gff3") #GFF3 with chromosome position
$no_targets = create_file("genes_without_target.txt") #TXT with a list of genes without the target


#===========================================================================================================================

# We add the headers to each file
$gff_genes.puts "##gff-version 3"
$gff_chr.puts "##gff-version 3"
$no_targets.puts "Report of genes without CTTCTT in their exons\n\n"


gene_file = load_file("ArabidopsisSubNetwork_GeneList.txt") #We read the file with the gene list


#===========================================================================================================================
#Main part of the program to search the targets in the exons
#===========================================================================================================================


gene_file.each do |gene|

  seq_obj = gene_info(gene) #This line creates the Bio:EMBL object from each gene

  unless seq_obj.nil?
    target_hash = get_exon_tar(seq_obj) # We search for the targets inside exons of the gene sequence

    if target_hash.empty?
      $no_targets.puts gene
      # If the target is not in the exons, we add it to the file genes_without_target.txt

    else
      feats(gene, target_hash, seq_obj) #We create new features and add them to each sequence object
      chr = get_chr(gene, seq_obj) #We return the chromosome information
      conv_to_chr(gene, target_hash, chr) #We convert all the info to chromosome positions

    end

  end

end


#===========================================================================================================================


puts "\nDONE! :D\n\n"
puts "The output is contained in the next files: "
puts "\t- genes_with_target.gff3"
puts "\t- chromosomes_with_target.gff3"
puts "\t- genes_without_target.txt"
puts "\n\nGoodbye!!"


#===========================================================================================================================
