#==============================================================#
#Bioinformatic programming challenges: Assignment 4            #
#Title: Searching for orthologues                              #
#Author: Adrián López Beltrán                                  #
#==============================================================#


require 'bio'
require 'stringio'
require 'io/console'
#We require the modules we need


#===============================================================================================


#===============================================================================================
#FUNCTION DECLARATION
#===============================================================================================


def create_file(file)
#Function that deletes a file with a given name and creates it again everytime it is called.
#Input: name of the file

  if File.exists? (file)
    File.delete(file)
  end

  return File.open(file, 'w')
end


#===============================================================================================
#MAIN CODE
#===============================================================================================

#To make the BLAST, we will use one of the files as a "search file", and other as a "target file". We will create BLAST database objects with them, and then will run the BLAST.
#Whit the BLAST finished,


#To continue or orthology study, one option would be to make a phylogenetic tree to infer the relationships beetween these genes.
#We could also try to constrast the information in some orthology database, such as https://www.genome.jp/kegg/ko.html


#WARNING: The target file 'arabidopsis.fa' (original name 'TAIR10_seq_20110103_representative_gene_model_updated' WILL NOT BE INCLUDED in the Github repository.

#===============================================================================================
#Global variables. References for the parameters found on:
#  1) https://www.ncbi.nlm.nih.gov/pubmed/18042555
#  2) https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4094424/


$eval_threshold = 10**-6
$overlap = 50


#===============================================================================================


#We create variables to introduce the FASTA file names and their types, in order to change them faster if we need to.
search = 'pep.fa'
search_type = 'prot'
target = 'arabidopsis.fa'
target_type = 'nucl'


#We create the file which will contain our list of orthologues
out_file = create_file('orthologues.txt')
out_file.puts "ORTHOLOGUES FOUND in files #{search} and #{target}\n"


#We change the name of the file variables so we can create the databases used in BLAST.
db_search = search.to_s + '_db'
db_target = target.to_s + '_db'


#We create the databases for BLAST, in a previoulsy created folder called "BLAST_db"
system("makeblastdb -in '#{search}' -dbtype #{search_type} -out ./BLAST_db/#{db_search}")
system("makeblastdb -in '#{target}' -dbtype #{target_type} -out ./BLAST_db/#{db_target}")


#We create the Bio::BLAST factory objects
factory_search = Bio::Blast.local('blastx', "./BLAST_db/#{db_search}")
factory_target = Bio::Blast.local('tblastn', "./BLAST_db/#{db_target}")


#We create a Bio::Fasta Format object with each of the original files
ff_search = Bio::FastaFormat.open(search)
ff_target = Bio::FastaFormat.open(target)


# We create a hash to store all sequences in the target file, and make them accesible by their ID.
target_hash = Hash.new
ff_target.each do |seq_target|
  target_hash[(seq_target.entry_id).to_s] = (seq_target.seq).to_s
end


#We will search for the best reciprocal hits using the results of the BLAST.

count = 1
ff_search.each do |seq_search| # We iterate over each sequence in the search

  search_id = (seq_search.entry_id).to_s # We store the ID in search to later know if it is a reciprocal best hit
  report_target = factory_target.query(seq_search)

  if report_target.hits[0] # Only if there have been hits continue.

    target_id = (report_target.hits[0].definition.match(/(\w+\.\w+)|/)).to_s # We get ID that will correspond to target ID

    if (report_target.hits[0].evalue <= $eval_threshold) and (report_target.hits[0].overlap >= $overlap) # We check the stablished parameters

      report_search = factory_search.query(">#{target_id}\n#{target_hash[target_id]}")
      # We look in the hash with the previous ID to get the sequence and query the factory

      if report_search.hits[0] # Again, only continue if there have been hits

        match = (report_search.hits[0].definition.match(/(\w+\.\w+)|/)).to_s # We get the ID that will match with the ID in the search

        if (report_search.hits[0].evalue <= $eval_threshold) and (report_search.hits[0].overlap >= $overlap) # Check parameters

          if search_id == match # If the match and the search ID match, it means that this is a reciprocal best hit

            out_file.puts "#{search_id}\t\t#{target_id}" # We write it in the output file

            count += 1

          end

        end

      end

    end

  end


end

#===============================================================================================


# We print the count of orthologues to our orthologues file.
out_file.puts "\n\nNumber of orthologues found: #{count}"


puts "DONE!!!!\n\n"
puts "You can browse the output in the file output_ortologues.txt"


out_file.close
